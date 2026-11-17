//! Shim for Solidcore's `surface`: the extension refusal reasons NLib reports.

pub mod extend {
    /// Why an extension was refused.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    #[non_exhaustive]
    pub enum ExtendRefusal {
        /// A weight of the extension would not be positive within the distance asked.
        WeightVanishes,
        /// A sphere's latitude would run past its pole.
        PastPole,
        /// A cone would run past its apex.
        PastApex,
        /// An angle would run over a full turn.
        OverATurn,
        /// The length across the side could not be measured to the stop.
        NotExact,
        /// A curve made longer keeping its control polygon's proportions whose growth in
        /// arc length is not certified within its tolerance of the length asked.
        LengthNotMet,
    }
}
