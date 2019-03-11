% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRRR8
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
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:36:19
% EndTime: 2019-03-09 02:36:19
% DurationCPUTime: 0.09s
% Computational Cost: add. (270->44), mult. (591->96), div. (0->0), fcn. (402->8), ass. (0->39)
t164 = qJD(1) ^ 2;
t172 = t164 / 0.2e1;
t171 = cos(qJ(5));
t157 = sin(pkin(10));
t149 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t167 = -pkin(7) * qJD(1) + t149;
t142 = t167 * t157;
t158 = cos(pkin(10));
t143 = t167 * t158;
t161 = sin(qJ(4));
t163 = cos(qJ(4));
t169 = t163 * t142 + t161 * t143;
t133 = qJD(4) * pkin(8) + t169;
t145 = (t157 * t163 + t158 * t161) * qJD(1);
t146 = (-t157 * t161 + t158 * t163) * qJD(1);
t151 = qJD(1) * qJ(2) + qJD(3);
t168 = qJD(1) * t157;
t147 = pkin(3) * t168 + t151;
t134 = t145 * pkin(4) - t146 * pkin(8) + t147;
t160 = sin(qJ(5));
t170 = t171 * t133 + t160 * t134;
t166 = -t160 * t133 + t171 * t134;
t165 = -t161 * t142 + t163 * t143;
t132 = -qJD(4) * pkin(4) - t165;
t144 = qJD(5) + t145;
t162 = cos(qJ(6));
t159 = sin(qJ(6));
t155 = t158 ^ 2;
t154 = t157 ^ 2;
t152 = -qJD(1) * pkin(1) + qJD(2);
t141 = qJD(6) + t144;
t137 = t160 * qJD(4) + t171 * t146;
t136 = -t171 * qJD(4) + t160 * t146;
t128 = -t159 * t136 + t162 * t137;
t127 = t162 * t136 + t159 * t137;
t126 = t136 * pkin(5) + t132;
t125 = -t136 * pkin(9) + t170;
t124 = t144 * pkin(5) - t137 * pkin(9) + t166;
t1 = [t172, 0, 0, t152 * qJD(1), t164 * qJ(2), qJ(2) ^ 2 * t172 + t152 ^ 2 / 0.2e1, t151 * t168, t151 * t158 * qJD(1) (-t154 - t155) * t149 * qJD(1), t151 ^ 2 / 0.2e1 + (t154 / 0.2e1 + t155 / 0.2e1) * t149 ^ 2, t146 ^ 2 / 0.2e1, -t146 * t145, t146 * qJD(4), -t145 * qJD(4), qJD(4) ^ 2 / 0.2e1, t165 * qJD(4) + t147 * t145, -t169 * qJD(4) + t147 * t146, t137 ^ 2 / 0.2e1, -t137 * t136, t137 * t144, -t136 * t144, t144 ^ 2 / 0.2e1, t132 * t136 + t166 * t144, t132 * t137 - t170 * t144, t128 ^ 2 / 0.2e1, -t128 * t127, t128 * t141, -t127 * t141, t141 ^ 2 / 0.2e1 (t162 * t124 - t159 * t125) * t141 + t126 * t127 -(t159 * t124 + t162 * t125) * t141 + t126 * t128;];
T_reg  = t1;
