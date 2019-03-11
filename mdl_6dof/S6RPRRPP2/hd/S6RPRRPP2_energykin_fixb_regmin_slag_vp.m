% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:33:20
% EndTime: 2019-03-09 04:33:20
% DurationCPUTime: 0.12s
% Computational Cost: add. (249->42), mult. (512->90), div. (0->0), fcn. (285->6), ass. (0->34)
t165 = pkin(4) + pkin(5);
t150 = qJD(1) ^ 2;
t164 = t150 / 0.2e1;
t163 = cos(qJ(4));
t145 = sin(pkin(9));
t137 = (pkin(1) * t145 + pkin(7)) * qJD(1);
t148 = sin(qJ(3));
t149 = cos(qJ(3));
t160 = t148 * qJD(2) + t149 * t137;
t131 = qJD(3) * pkin(8) + t160;
t146 = cos(pkin(9));
t156 = -pkin(1) * t146 - pkin(2);
t132 = (-pkin(3) * t149 - pkin(8) * t148 + t156) * qJD(1);
t147 = sin(qJ(4));
t162 = t163 * t131 + t147 * t132;
t161 = t149 * qJD(2) - t148 * t137;
t159 = qJD(1) * t148;
t158 = t149 * qJD(1);
t157 = qJD(1) * qJD(3);
t140 = -qJD(4) + t158;
t125 = -t140 * qJ(5) + t162;
t155 = qJD(3) * pkin(3) + t161;
t154 = -t147 * t131 + t163 * t132;
t153 = qJD(5) - t154;
t136 = t147 * qJD(3) + t163 * t159;
t152 = t136 * qJ(5) + t155;
t138 = t156 * qJD(1);
t135 = -t163 * qJD(3) + t147 * t159;
t126 = t135 * pkin(4) - t152;
t124 = t140 * pkin(4) + t153;
t123 = -t165 * t135 + qJD(6) + t152;
t122 = t135 * qJ(6) + t125;
t121 = -t136 * qJ(6) + t165 * t140 + t153;
t1 = [t164, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t145 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t150, t148 ^ 2 * t164, t148 * t150 * t149, t148 * t157, t149 * t157, qJD(3) ^ 2 / 0.2e1, t161 * qJD(3) - t138 * t158, -t160 * qJD(3) + t138 * t159, t136 ^ 2 / 0.2e1, -t136 * t135, -t136 * t140, t135 * t140, t140 ^ 2 / 0.2e1, -t135 * t155 - t154 * t140, -t136 * t155 + t162 * t140, t124 * t140 + t126 * t135, t124 * t136 - t125 * t135, -t125 * t140 - t126 * t136, t125 ^ 2 / 0.2e1 + t126 ^ 2 / 0.2e1 + t124 ^ 2 / 0.2e1, t121 * t140 - t123 * t135, -t122 * t140 + t123 * t136, -t121 * t136 + t122 * t135, t122 ^ 2 / 0.2e1 + t121 ^ 2 / 0.2e1 + t123 ^ 2 / 0.2e1;];
T_reg  = t1;
