% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:01:04
% EndTime: 2019-03-09 16:01:04
% DurationCPUTime: 0.13s
% Computational Cost: add. (441->53), mult. (906->107), div. (0->0), fcn. (597->8), ass. (0->42)
t198 = qJD(1) ^ 2;
t209 = t198 / 0.2e1;
t208 = cos(pkin(10));
t197 = cos(qJ(2));
t207 = t197 * t198;
t193 = sin(qJ(3));
t196 = cos(qJ(3));
t194 = sin(qJ(2));
t205 = qJD(1) * t194;
t179 = qJD(2) * t193 + t196 * t205;
t204 = t197 * qJD(1);
t186 = -qJD(3) + t204;
t175 = (-pkin(2) * t197 - pkin(8) * t194 - pkin(1)) * qJD(1);
t183 = pkin(7) * t204 + qJD(2) * pkin(8);
t200 = t196 * t175 - t183 * t193;
t199 = qJD(4) - t200;
t163 = -t179 * qJ(5) + (pkin(3) + pkin(4)) * t186 + t199;
t206 = t175 * t193 + t183 * t196;
t168 = -t186 * qJ(4) + t206;
t178 = -qJD(2) * t196 + t193 * t205;
t165 = qJ(5) * t178 + t168;
t190 = sin(pkin(10));
t157 = t163 * t190 + t165 * t208;
t182 = -qJD(2) * pkin(2) + pkin(7) * t205;
t203 = qJD(1) * qJD(2);
t202 = t194 * t203;
t201 = t197 * t203;
t156 = t163 * t208 - t165 * t190;
t169 = pkin(3) * t178 - qJ(4) * t179 + t182;
t166 = -pkin(4) * t178 + qJD(5) - t169;
t195 = cos(qJ(6));
t192 = sin(qJ(6));
t185 = qJD(6) + t186;
t172 = t178 * t190 + t179 * t208;
t171 = -t178 * t208 + t179 * t190;
t167 = pkin(3) * t186 + t199;
t160 = -t171 * t192 + t172 * t195;
t159 = t171 * t195 + t172 * t192;
t158 = pkin(5) * t171 + t166;
t155 = -pkin(9) * t171 + t157;
t154 = pkin(5) * t186 - pkin(9) * t172 + t156;
t1 = [t209, 0, 0, t194 ^ 2 * t209, t194 * t207, t202, t201, qJD(2) ^ 2 / 0.2e1, pkin(1) * t207 - pkin(7) * t202, -pkin(1) * t194 * t198 - pkin(7) * t201, t179 ^ 2 / 0.2e1, -t179 * t178, -t179 * t186, t178 * t186, t186 ^ 2 / 0.2e1, t178 * t182 - t186 * t200, t179 * t182 + t186 * t206, t167 * t186 + t169 * t178, t167 * t179 - t168 * t178, -t168 * t186 - t169 * t179, t168 ^ 2 / 0.2e1 + t169 ^ 2 / 0.2e1 + t167 ^ 2 / 0.2e1, t156 * t186 + t166 * t171, -t157 * t186 + t166 * t172, -t156 * t172 - t157 * t171, t157 ^ 2 / 0.2e1 + t156 ^ 2 / 0.2e1 + t166 ^ 2 / 0.2e1, t160 ^ 2 / 0.2e1, -t160 * t159, t160 * t185, -t159 * t185, t185 ^ 2 / 0.2e1 (t154 * t195 - t155 * t192) * t185 + t158 * t159 -(t154 * t192 + t155 * t195) * t185 + t158 * t160;];
T_reg  = t1;
