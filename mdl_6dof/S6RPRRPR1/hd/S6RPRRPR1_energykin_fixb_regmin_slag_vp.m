% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:58:44
% EndTime: 2019-03-09 04:58:44
% DurationCPUTime: 0.13s
% Computational Cost: add. (335->47), mult. (765->103), div. (0->0), fcn. (525->10), ass. (0->42)
t176 = qJD(1) ^ 2;
t186 = t176 / 0.2e1;
t185 = cos(qJ(4));
t172 = sin(qJ(4));
t173 = sin(qJ(3));
t175 = cos(qJ(3));
t157 = (t172 * t175 + t173 * t185) * qJD(1);
t166 = qJD(3) + qJD(4);
t168 = sin(pkin(10));
t160 = (pkin(1) * t168 + pkin(7)) * qJD(1);
t165 = t175 * qJD(2);
t153 = qJD(3) * pkin(3) + t165 + (-pkin(8) * qJD(1) - t160) * t173;
t181 = qJD(1) * t175;
t183 = qJD(2) * t173 + t160 * t175;
t154 = pkin(8) * t181 + t183;
t178 = t153 * t185 - t154 * t172;
t141 = pkin(4) * t166 - qJ(5) * t157 + t178;
t182 = qJD(1) * t173;
t156 = t172 * t182 - t181 * t185;
t184 = t153 * t172 + t154 * t185;
t143 = -qJ(5) * t156 + t184;
t167 = sin(pkin(11));
t169 = cos(pkin(11));
t138 = t141 * t167 + t143 * t169;
t180 = qJD(1) * qJD(3);
t170 = cos(pkin(10));
t179 = -pkin(1) * t170 - pkin(2);
t147 = -t156 * t169 - t157 * t167;
t137 = t141 * t169 - t143 * t167;
t158 = (-pkin(3) * t175 + t179) * qJD(1);
t149 = t156 * pkin(4) + qJD(5) + t158;
t174 = cos(qJ(6));
t171 = sin(qJ(6));
t161 = t179 * qJD(1);
t148 = -t156 * t167 + t157 * t169;
t146 = qJD(6) - t147;
t145 = t148 * t174 + t166 * t171;
t144 = t148 * t171 - t166 * t174;
t139 = -t147 * pkin(5) - t148 * pkin(9) + t149;
t136 = pkin(9) * t166 + t138;
t135 = -pkin(5) * t166 - t137;
t1 = [t186, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t168 ^ 2 / 0.2e1 + t170 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t176, t173 ^ 2 * t186, t173 * t176 * t175, t173 * t180, t175 * t180, qJD(3) ^ 2 / 0.2e1, -t161 * t181 + (-t160 * t173 + t165) * qJD(3), -qJD(3) * t183 + t161 * t182, t157 ^ 2 / 0.2e1, -t157 * t156, t157 * t166, -t156 * t166, t166 ^ 2 / 0.2e1, t156 * t158 + t166 * t178, t157 * t158 - t166 * t184, -t137 * t148 + t138 * t147, t138 ^ 2 / 0.2e1 + t137 ^ 2 / 0.2e1 + t149 ^ 2 / 0.2e1, t145 ^ 2 / 0.2e1, -t145 * t144, t145 * t146, -t144 * t146, t146 ^ 2 / 0.2e1 (-t136 * t171 + t139 * t174) * t146 + t135 * t144 -(t136 * t174 + t139 * t171) * t146 + t135 * t145;];
T_reg  = t1;
