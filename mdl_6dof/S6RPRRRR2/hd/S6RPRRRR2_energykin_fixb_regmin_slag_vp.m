% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRR2
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
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:58:29
% EndTime: 2019-03-09 06:58:29
% DurationCPUTime: 0.12s
% Computational Cost: add. (317->48), mult. (674->106), div. (0->0), fcn. (469->10), ass. (0->43)
t176 = qJD(1) ^ 2;
t188 = t176 / 0.2e1;
t187 = cos(qJ(5));
t166 = qJD(3) + qJD(4);
t167 = sin(pkin(11));
t160 = (pkin(1) * t167 + pkin(7)) * qJD(1);
t175 = cos(qJ(3));
t165 = t175 * qJD(2);
t172 = sin(qJ(3));
t150 = qJD(3) * pkin(3) + t165 + (-pkin(8) * qJD(1) - t160) * t172;
t182 = qJD(1) * t175;
t184 = t172 * qJD(2) + t175 * t160;
t153 = pkin(8) * t182 + t184;
t171 = sin(qJ(4));
t174 = cos(qJ(4));
t185 = t171 * t150 + t174 * t153;
t140 = t166 * pkin(9) + t185;
t183 = qJD(1) * t172;
t156 = t171 * t183 - t174 * t182;
t157 = (t171 * t175 + t172 * t174) * qJD(1);
t168 = cos(pkin(11));
t180 = -pkin(1) * t168 - pkin(2);
t158 = (-pkin(3) * t175 + t180) * qJD(1);
t145 = t156 * pkin(4) - t157 * pkin(9) + t158;
t170 = sin(qJ(5));
t186 = t187 * t140 + t170 * t145;
t181 = qJD(1) * qJD(3);
t179 = -t170 * t140 + t187 * t145;
t178 = t174 * t150 - t171 * t153;
t139 = -t166 * pkin(4) - t178;
t155 = qJD(5) + t156;
t173 = cos(qJ(6));
t169 = sin(qJ(6));
t161 = t180 * qJD(1);
t154 = qJD(6) + t155;
t152 = t187 * t157 + t170 * t166;
t151 = t170 * t157 - t187 * t166;
t142 = -t169 * t151 + t173 * t152;
t141 = t173 * t151 + t169 * t152;
t137 = t151 * pkin(5) + t139;
t136 = -t151 * pkin(10) + t186;
t135 = t155 * pkin(5) - t152 * pkin(10) + t179;
t1 = [t188, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t167 ^ 2 / 0.2e1 + t168 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t176, t172 ^ 2 * t188, t172 * t176 * t175, t172 * t181, t175 * t181, qJD(3) ^ 2 / 0.2e1, -t161 * t182 + (-t172 * t160 + t165) * qJD(3), -t184 * qJD(3) + t161 * t183, t157 ^ 2 / 0.2e1, -t157 * t156, t166 * t157, -t166 * t156, t166 ^ 2 / 0.2e1, t158 * t156 + t178 * t166, t158 * t157 - t185 * t166, t152 ^ 2 / 0.2e1, -t152 * t151, t152 * t155, -t151 * t155, t155 ^ 2 / 0.2e1, t139 * t151 + t179 * t155, t139 * t152 - t186 * t155, t142 ^ 2 / 0.2e1, -t142 * t141, t142 * t154, -t141 * t154, t154 ^ 2 / 0.2e1 (t173 * t135 - t169 * t136) * t154 + t137 * t141 -(t169 * t135 + t173 * t136) * t154 + t137 * t142;];
T_reg  = t1;
