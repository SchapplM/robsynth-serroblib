% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:26:42
% EndTime: 2019-03-09 15:26:42
% DurationCPUTime: 0.14s
% Computational Cost: add. (440->50), mult. (1049->102), div. (0->0), fcn. (747->8), ass. (0->45)
t194 = pkin(4) + pkin(9);
t193 = -pkin(8) - pkin(7);
t180 = qJD(1) ^ 2;
t192 = t180 / 0.2e1;
t191 = cos(qJ(3));
t179 = cos(qJ(2));
t190 = t179 * t180;
t176 = sin(qJ(3));
t177 = sin(qJ(2));
t166 = (t176 * t179 + t191 * t177) * qJD(1);
t172 = qJD(2) + qJD(3);
t188 = qJD(1) * t177;
t168 = qJD(2) * pkin(2) + t193 * t188;
t187 = qJD(1) * t179;
t169 = t193 * t187;
t183 = t191 * t168 + t176 * t169;
t151 = t172 * pkin(3) - t166 * qJ(4) + t183;
t165 = t176 * t188 - t191 * t187;
t189 = t176 * t168 - t191 * t169;
t154 = -t165 * qJ(4) + t189;
t173 = sin(pkin(10));
t174 = cos(pkin(10));
t148 = t173 * t151 + t174 * t154;
t186 = qJD(1) * qJD(2);
t185 = t177 * t186;
t184 = t179 * t186;
t147 = t174 * t151 - t173 * t154;
t146 = -t172 * qJ(5) - t148;
t170 = (-pkin(2) * t179 - pkin(1)) * qJD(1);
t182 = qJD(5) - t147;
t160 = -t173 * t165 + t174 * t166;
t161 = t165 * pkin(3) + qJD(4) + t170;
t181 = -t160 * qJ(5) + t161;
t178 = cos(qJ(6));
t175 = sin(qJ(6));
t159 = t174 * t165 + t173 * t166;
t158 = qJD(6) + t160;
t156 = t175 * t159 + t178 * t172;
t155 = -t178 * t159 + t175 * t172;
t149 = t159 * pkin(4) + t181;
t145 = -t172 * pkin(4) + t182;
t144 = t194 * t159 + t181;
t143 = -t159 * pkin(5) - t146;
t142 = t160 * pkin(5) - t194 * t172 + t182;
t1 = [t192, 0, 0, t177 ^ 2 * t192, t177 * t190, t185, t184, qJD(2) ^ 2 / 0.2e1, pkin(1) * t190 - pkin(7) * t185, -t180 * pkin(1) * t177 - pkin(7) * t184, t166 ^ 2 / 0.2e1, -t166 * t165, t166 * t172, -t165 * t172, t172 ^ 2 / 0.2e1, t170 * t165 + t183 * t172, t170 * t166 - t189 * t172, -t147 * t160 - t148 * t159, t148 ^ 2 / 0.2e1 + t147 ^ 2 / 0.2e1 + t161 ^ 2 / 0.2e1, t145 * t160 + t146 * t159, t145 * t172 - t149 * t159, -t146 * t172 - t149 * t160, t149 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1 + t145 ^ 2 / 0.2e1, t156 ^ 2 / 0.2e1, -t156 * t155, t156 * t158, -t155 * t158, t158 ^ 2 / 0.2e1 (t178 * t142 - t175 * t144) * t158 + t143 * t155 -(t175 * t142 + t178 * t144) * t158 + t143 * t156;];
T_reg  = t1;
