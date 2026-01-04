//! Simple error handling without anyhow/thiserror

use std::error::Error;
use std::fmt;
use std::io;

#[derive(Debug)]
pub struct AppError {
    message: String,
    source: Option<Box<dyn Error + Send + Sync>>,
}

impl AppError {
    pub fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
            source: None,
        }
    }

    pub fn with_source(message: impl Into<String>, source: impl Error + Send + Sync + 'static) -> Self {
        Self {
            message: message.into(),
            source: Some(Box::new(source)),
        }
    }
}

impl fmt::Display for AppError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.message)?;
        if let Some(ref source) = self.source {
            write!(f, ": {}", source)?;
        }
        Ok(())
    }
}

impl Error for AppError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        self.source.as_ref().map(|e| e.as_ref() as &(dyn Error + 'static))
    }
}

impl From<io::Error> for AppError {
    fn from(err: io::Error) -> Self {
        AppError::with_source("IO error", err)
    }
}

impl From<String> for AppError {
    fn from(s: String) -> Self {
        AppError::new(s)
    }
}

impl From<&str> for AppError {
    fn from(s: &str) -> Self {
        AppError::new(s)
    }
}

pub type Result<T> = std::result::Result<T, AppError>;

/// Context trait for adding context to errors
pub trait Context<T> {
    fn context(self, msg: impl Into<String>) -> Result<T>;
}

impl<T, E: Error + Send + Sync + 'static> Context<T> for std::result::Result<T, E> {
    fn context(self, msg: impl Into<String>) -> Result<T> {
        self.map_err(|e| AppError::with_source(msg, e))
    }
}

impl<T> Context<T> for Option<T> {
    fn context(self, msg: impl Into<String>) -> Result<T> {
        self.ok_or_else(|| AppError::new(msg))
    }
}
