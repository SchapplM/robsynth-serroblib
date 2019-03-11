% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:26:07
% EndTime: 2019-03-09 09:26:07
% DurationCPUTime: 0.12s
% Computational Cost: add. (354->52), mult. (805->107), div. (0->0), fcn. (539->8), ass. (0->41)
t197 = qJD(1) ^ 2;
t207 = t197 / 0.2e1;
t206 = cos(qJ(5));
t196 = cos(qJ(2));
t205 = t196 * t197;
t194 = sin(qJ(2));
t176 = (-pkin(2) * t196 - qJ(3) * t194 - pkin(1)) * qJD(1);
t202 = t196 * qJD(1);
t183 = pkin(7) * t202 + qJD(2) * qJ(3);
t189 = sin(pkin(10));
t190 = cos(pkin(10));
t171 = t176 * t190 - t183 * t189;
t166 = pkin(3) * t202 + qJD(4) - t171;
t203 = qJD(1) * t194;
t179 = qJD(2) * t189 + t190 * t203;
t161 = pkin(4) * t202 - pkin(8) * t179 + t166;
t172 = t176 * t189 + t183 * t190;
t167 = -qJ(4) * t202 + t172;
t178 = -qJD(2) * t190 + t189 * t203;
t164 = pkin(8) * t178 + t167;
t193 = sin(qJ(5));
t204 = t161 * t193 + t164 * t206;
t201 = qJD(1) * qJD(2);
t182 = -qJD(2) * pkin(2) + pkin(7) * t203 + qJD(3);
t200 = t194 * t201;
t199 = t196 * t201;
t198 = t161 * t206 - t164 * t193;
t185 = qJD(5) + t202;
t165 = pkin(3) * t178 - qJ(4) * t179 + t182;
t163 = -pkin(4) * t178 - t165;
t195 = cos(qJ(6));
t192 = sin(qJ(6));
t184 = qJD(6) + t185;
t170 = t178 * t193 + t179 * t206;
t169 = -t178 * t206 + t179 * t193;
t158 = -t169 * t192 + t170 * t195;
t157 = t169 * t195 + t170 * t192;
t156 = pkin(5) * t169 + t163;
t155 = -pkin(9) * t169 + t204;
t154 = pkin(5) * t185 - pkin(9) * t170 + t198;
t1 = [t207, 0, 0, t194 ^ 2 * t207, t194 * t205, t200, t199, qJD(2) ^ 2 / 0.2e1, pkin(1) * t205 - pkin(7) * t200, -pkin(1) * t194 * t197 - pkin(7) * t199, -t171 * t202 + t178 * t182, t172 * t202 + t179 * t182, -t171 * t179 - t172 * t178, t172 ^ 2 / 0.2e1 + t171 ^ 2 / 0.2e1 + t182 ^ 2 / 0.2e1, t165 * t178 + t166 * t202, t166 * t179 - t167 * t178, -t165 * t179 - t167 * t202, t167 ^ 2 / 0.2e1 + t165 ^ 2 / 0.2e1 + t166 ^ 2 / 0.2e1, t170 ^ 2 / 0.2e1, -t170 * t169, t170 * t185, -t169 * t185, t185 ^ 2 / 0.2e1, t163 * t169 + t185 * t198, t163 * t170 - t185 * t204, t158 ^ 2 / 0.2e1, -t158 * t157, t158 * t184, -t157 * t184, t184 ^ 2 / 0.2e1 (t154 * t195 - t155 * t192) * t184 + t156 * t157 -(t154 * t192 + t155 * t195) * t184 + t156 * t158;];
T_reg  = t1;
