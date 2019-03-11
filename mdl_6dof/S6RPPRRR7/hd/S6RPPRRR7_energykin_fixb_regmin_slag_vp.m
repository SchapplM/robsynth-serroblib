% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% T_reg [1x31]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:47
% EndTime: 2019-03-09 02:33:47
% DurationCPUTime: 0.12s
% Computational Cost: add. (291->44), mult. (645->96), div. (0->0), fcn. (448->8), ass. (0->39)
t159 = qJD(1) ^ 2;
t166 = t159 / 0.2e1;
t151 = sin(pkin(10));
t152 = cos(pkin(10));
t155 = sin(qJ(4));
t158 = cos(qJ(4));
t139 = (-t151 * t155 + t152 * t158) * qJD(1);
t142 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t162 = -pkin(7) * qJD(1) + t142;
t136 = t162 * t151;
t137 = t162 * t152;
t161 = -t155 * t136 + t158 * t137;
t124 = qJD(4) * pkin(4) - t139 * pkin(8) + t161;
t138 = (t151 * t158 + t152 * t155) * qJD(1);
t164 = t158 * t136 + t155 * t137;
t125 = -t138 * pkin(8) + t164;
t154 = sin(qJ(5));
t157 = cos(qJ(5));
t165 = t154 * t124 + t157 * t125;
t163 = qJD(1) * t151;
t145 = qJD(1) * qJ(2) + qJD(3);
t140 = pkin(3) * t163 + t145;
t129 = t157 * t138 + t154 * t139;
t160 = t157 * t124 - t154 * t125;
t131 = t138 * pkin(4) + t140;
t156 = cos(qJ(6));
t153 = sin(qJ(6));
t149 = qJD(4) + qJD(5);
t148 = t152 ^ 2;
t147 = t151 ^ 2;
t146 = -pkin(1) * qJD(1) + qJD(2);
t130 = -t154 * t138 + t157 * t139;
t128 = qJD(6) + t129;
t127 = t156 * t130 + t153 * t149;
t126 = t153 * t130 - t156 * t149;
t121 = t129 * pkin(5) - t130 * pkin(9) + t131;
t120 = t149 * pkin(9) + t165;
t119 = -t149 * pkin(5) - t160;
t1 = [t166, 0, 0, t146 * qJD(1), t159 * qJ(2), qJ(2) ^ 2 * t166 + t146 ^ 2 / 0.2e1, t145 * t163, t145 * t152 * qJD(1) (-t147 - t148) * t142 * qJD(1), t145 ^ 2 / 0.2e1 + (t147 / 0.2e1 + t148 / 0.2e1) * t142 ^ 2, t139 ^ 2 / 0.2e1, -t139 * t138, t139 * qJD(4), -t138 * qJD(4), qJD(4) ^ 2 / 0.2e1, t161 * qJD(4) + t140 * t138, -t164 * qJD(4) + t140 * t139, t130 ^ 2 / 0.2e1, -t130 * t129, t130 * t149, -t129 * t149, t149 ^ 2 / 0.2e1, t131 * t129 + t160 * t149, t131 * t130 - t165 * t149, t127 ^ 2 / 0.2e1, -t127 * t126, t127 * t128, -t126 * t128, t128 ^ 2 / 0.2e1 (-t153 * t120 + t156 * t121) * t128 + t119 * t126 -(t156 * t120 + t153 * t121) * t128 + t119 * t127;];
T_reg  = t1;
