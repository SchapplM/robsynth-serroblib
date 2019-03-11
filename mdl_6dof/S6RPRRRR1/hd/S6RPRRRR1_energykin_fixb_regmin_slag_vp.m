% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:55:16
% EndTime: 2019-03-09 06:55:16
% DurationCPUTime: 0.13s
% Computational Cost: add. (325->48), mult. (734->106), div. (0->0), fcn. (521->10), ass. (0->43)
t170 = qJD(1) ^ 2;
t182 = t170 / 0.2e1;
t181 = cos(qJ(4));
t165 = sin(qJ(4));
t166 = sin(qJ(3));
t169 = cos(qJ(3));
t150 = (t165 * t169 + t181 * t166) * qJD(1);
t160 = qJD(3) + qJD(4);
t161 = sin(pkin(11));
t153 = (pkin(1) * t161 + pkin(7)) * qJD(1);
t159 = t169 * qJD(2);
t146 = qJD(3) * pkin(3) + t159 + (-pkin(8) * qJD(1) - t153) * t166;
t176 = qJD(1) * t169;
t178 = t166 * qJD(2) + t169 * t153;
t147 = pkin(8) * t176 + t178;
t173 = t181 * t146 - t165 * t147;
t134 = t160 * pkin(4) - t150 * pkin(9) + t173;
t177 = qJD(1) * t166;
t149 = t165 * t177 - t181 * t176;
t179 = t165 * t146 + t181 * t147;
t136 = -t149 * pkin(9) + t179;
t164 = sin(qJ(5));
t168 = cos(qJ(5));
t180 = t164 * t134 + t168 * t136;
t175 = qJD(1) * qJD(3);
t162 = cos(pkin(11));
t174 = -pkin(1) * t162 - pkin(2);
t140 = t168 * t149 + t164 * t150;
t172 = t168 * t134 - t164 * t136;
t151 = (-pkin(3) * t169 + t174) * qJD(1);
t142 = t149 * pkin(4) + t151;
t167 = cos(qJ(6));
t163 = sin(qJ(6));
t157 = qJD(5) + t160;
t154 = t174 * qJD(1);
t141 = -t164 * t149 + t168 * t150;
t139 = qJD(6) + t140;
t138 = t167 * t141 + t163 * t157;
t137 = t163 * t141 - t167 * t157;
t132 = t140 * pkin(5) - t141 * pkin(10) + t142;
t131 = t157 * pkin(10) + t180;
t130 = -t157 * pkin(5) - t172;
t1 = [t182, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t161 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t170, t166 ^ 2 * t182, t166 * t170 * t169, t166 * t175, t169 * t175, qJD(3) ^ 2 / 0.2e1, -t154 * t176 + (-t166 * t153 + t159) * qJD(3), -t178 * qJD(3) + t154 * t177, t150 ^ 2 / 0.2e1, -t150 * t149, t150 * t160, -t149 * t160, t160 ^ 2 / 0.2e1, t151 * t149 + t173 * t160, t151 * t150 - t179 * t160, t141 ^ 2 / 0.2e1, -t141 * t140, t141 * t157, -t140 * t157, t157 ^ 2 / 0.2e1, t142 * t140 + t172 * t157, t142 * t141 - t180 * t157, t138 ^ 2 / 0.2e1, -t138 * t137, t138 * t139, -t137 * t139, t139 ^ 2 / 0.2e1 (-t163 * t131 + t167 * t132) * t139 + t130 * t137 -(t167 * t131 + t163 * t132) * t139 + t130 * t138;];
T_reg  = t1;
