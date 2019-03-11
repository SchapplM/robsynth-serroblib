% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:59:45
% EndTime: 2019-03-09 13:59:45
% DurationCPUTime: 0.13s
% Computational Cost: add. (325->52), mult. (667->109), div. (0->0), fcn. (437->8), ass. (0->43)
t196 = qJD(1) ^ 2;
t209 = t196 / 0.2e1;
t208 = cos(qJ(5));
t195 = cos(qJ(2));
t207 = t195 * t196;
t203 = qJD(1) * t195;
t192 = sin(qJ(2));
t204 = qJD(1) * t192;
t176 = -qJD(1) * pkin(1) - pkin(2) * t203 - qJ(3) * t204;
t167 = pkin(3) * t203 - t176;
t191 = sin(qJ(4));
t194 = cos(qJ(4));
t173 = (t191 * t192 + t194 * t195) * qJD(1);
t174 = (-t191 * t195 + t192 * t194) * qJD(1);
t159 = t173 * pkin(4) - t174 * pkin(9) + t167;
t185 = qJD(2) - qJD(4);
t202 = pkin(7) * t204 + qJD(3);
t168 = -pkin(8) * t204 + (-pkin(2) - pkin(3)) * qJD(2) + t202;
t178 = pkin(7) * t203 + qJD(2) * qJ(3);
t175 = -pkin(8) * t203 + t178;
t205 = t191 * t168 + t194 * t175;
t162 = -t185 * pkin(9) + t205;
t190 = sin(qJ(5));
t206 = t190 * t159 + t208 * t162;
t201 = qJD(1) * qJD(2);
t200 = t192 * t201;
t199 = t195 * t201;
t198 = t208 * t159 - t190 * t162;
t197 = t194 * t168 - t191 * t175;
t161 = t185 * pkin(4) - t197;
t172 = qJD(5) + t173;
t193 = cos(qJ(6));
t189 = sin(qJ(6));
t177 = -qJD(2) * pkin(2) + t202;
t171 = qJD(6) + t172;
t165 = t208 * t174 - t190 * t185;
t164 = t190 * t174 + t208 * t185;
t156 = -t189 * t164 + t193 * t165;
t155 = t193 * t164 + t189 * t165;
t154 = t164 * pkin(5) + t161;
t153 = -t164 * pkin(10) + t206;
t152 = t172 * pkin(5) - t165 * pkin(10) + t198;
t1 = [t209, 0, 0, t192 ^ 2 * t209, t192 * t207, t200, t199, qJD(2) ^ 2 / 0.2e1, pkin(1) * t207 - pkin(7) * t200, -t196 * pkin(1) * t192 - pkin(7) * t199, -t177 * qJD(2) - t176 * t203 (t177 * t192 + t178 * t195) * qJD(1), t178 * qJD(2) - t176 * t204, t178 ^ 2 / 0.2e1 + t176 ^ 2 / 0.2e1 + t177 ^ 2 / 0.2e1, t174 ^ 2 / 0.2e1, -t174 * t173, -t174 * t185, t173 * t185, t185 ^ 2 / 0.2e1, t167 * t173 - t197 * t185, t167 * t174 + t205 * t185, t165 ^ 2 / 0.2e1, -t165 * t164, t165 * t172, -t164 * t172, t172 ^ 2 / 0.2e1, t161 * t164 + t198 * t172, t161 * t165 - t206 * t172, t156 ^ 2 / 0.2e1, -t156 * t155, t156 * t171, -t155 * t171, t171 ^ 2 / 0.2e1 (t193 * t152 - t189 * t153) * t171 + t154 * t155 -(t189 * t152 + t193 * t153) * t171 + t154 * t156;];
T_reg  = t1;
