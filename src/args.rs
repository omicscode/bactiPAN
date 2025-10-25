use clap::{Parser, Subcommand};
#[derive(Debug, Parser)]
#[command(
    name = "bactiPAN",
    version = "1.0",
    about = "bacterial Short and Long Read Pangenome.
       ************************************************
       Gaurav Sablok,
       Email: codeprog@icloud.com
      ************************************************"
)]
pub struct CommandParse {
    /// subcommands for the specific actions
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// form input
    Selectoption {
        /// input string
        optionvalue: String,
    },
}
