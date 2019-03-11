% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:21:21
% EndTime: 2019-03-09 02:21:21
% DurationCPUTime: 0.12s
% Computational Cost: add. (281->50), mult. (670->103), div. (0->0), fcn. (473->10), ass. (0->41)
t186 = cos(qJ(5));
t169 = sin(pkin(10));
t163 = (pkin(1) * t169 + qJ(3)) * qJD(1);
t170 = cos(pkin(11));
t166 = t170 * qJD(2);
t168 = sin(pkin(11));
t150 = t166 + (-pkin(7) * qJD(1) - t163) * t168;
t155 = qJD(2) * t168 + t163 * t170;
t182 = qJD(1) * t170;
t151 = pkin(7) * t182 + t155;
t174 = sin(qJ(4));
t176 = cos(qJ(4));
t184 = t150 * t174 + t151 * t176;
t140 = qJD(4) * pkin(8) + t184;
t171 = cos(pkin(10));
t181 = -pkin(1) * t171 - pkin(2);
t158 = qJD(3) + (-pkin(3) * t170 + t181) * qJD(1);
t183 = qJD(1) * t168;
t159 = t174 * t183 - t176 * t182;
t160 = (t168 * t176 + t170 * t174) * qJD(1);
t145 = t159 * pkin(4) - t160 * pkin(8) + t158;
t173 = sin(qJ(5));
t185 = t140 * t186 + t145 * t173;
t180 = -t140 * t173 + t145 * t186;
t179 = t150 * t176 - t151 * t174;
t157 = qJD(5) + t159;
t139 = -qJD(4) * pkin(4) - t179;
t177 = qJD(1) ^ 2;
t175 = cos(qJ(6));
t172 = sin(qJ(6));
t162 = qJD(1) * t181 + qJD(3);
t156 = qJD(6) + t157;
t154 = -t163 * t168 + t166;
t153 = qJD(4) * t173 + t160 * t186;
t152 = -qJD(4) * t186 + t160 * t173;
t142 = -t152 * t172 + t153 * t175;
t141 = t152 * t175 + t153 * t172;
t137 = pkin(5) * t152 + t139;
t136 = -pkin(9) * t152 + t185;
t135 = pkin(5) * t157 - pkin(9) * t153 + t180;
t1 = [t177 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t169 ^ 2 / 0.2e1 + t171 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t177, -t162 * t182, t162 * t183 (-t154 * t168 + t155 * t170) * qJD(1), t155 ^ 2 / 0.2e1 + t154 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1, t160 ^ 2 / 0.2e1, -t159 * t160, qJD(4) * t160, -t159 * qJD(4), qJD(4) ^ 2 / 0.2e1, qJD(4) * t179 + t158 * t159, -qJD(4) * t184 + t158 * t160, t153 ^ 2 / 0.2e1, -t153 * t152, t153 * t157, -t152 * t157, t157 ^ 2 / 0.2e1, t139 * t152 + t157 * t180, t139 * t153 - t157 * t185, t142 ^ 2 / 0.2e1, -t142 * t141, t142 * t156, -t141 * t156, t156 ^ 2 / 0.2e1 (t135 * t175 - t136 * t172) * t156 + t137 * t141 -(t135 * t172 + t136 * t175) * t156 + t137 * t142;];
T_reg  = t1;
