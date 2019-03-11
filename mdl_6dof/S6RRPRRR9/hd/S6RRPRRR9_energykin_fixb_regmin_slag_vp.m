% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:13:58
% EndTime: 2019-03-09 14:13:58
% DurationCPUTime: 0.18s
% Computational Cost: add. (859->58), mult. (2065->120), div. (0->0), fcn. (1687->12), ass. (0->52)
t253 = cos(qJ(4));
t229 = sin(pkin(6));
t239 = qJD(1) ^ 2;
t252 = t229 ^ 2 * t239;
t247 = cos(pkin(6)) * qJD(1);
t226 = qJD(2) + t247;
t228 = sin(pkin(12));
t230 = cos(pkin(12));
t235 = sin(qJ(2));
t248 = qJD(1) * t229;
t244 = t235 * t248;
t215 = -t230 * t226 + t228 * t244;
t216 = t228 * t226 + t230 * t244;
t234 = sin(qJ(4));
t207 = -t234 * t215 + t216 * t253;
t238 = cos(qJ(2));
t243 = t238 * t248;
t221 = -qJD(4) + t243;
t246 = pkin(1) * t247;
t249 = pkin(8) * t243 + t235 * t246;
t212 = t226 * qJ(3) + t249;
t214 = (-pkin(2) * t238 - qJ(3) * t235 - pkin(1)) * t248;
t202 = -t228 * t212 + t230 * t214;
t199 = -pkin(3) * t243 - t216 * pkin(9) + t202;
t203 = t230 * t212 + t228 * t214;
t201 = -t215 * pkin(9) + t203;
t242 = t199 * t253 - t234 * t201;
t188 = -t221 * pkin(4) - t207 * pkin(10) + t242;
t206 = t215 * t253 + t234 * t216;
t250 = t234 * t199 + t201 * t253;
t190 = -t206 * pkin(10) + t250;
t233 = sin(qJ(5));
t237 = cos(qJ(5));
t251 = t233 * t188 + t237 * t190;
t245 = t238 * t252;
t194 = t237 * t206 + t233 * t207;
t241 = t237 * t188 - t233 * t190;
t240 = -pkin(8) * t244 + t238 * t246;
t209 = -t226 * pkin(2) + qJD(3) - t240;
t205 = t215 * pkin(3) + t209;
t196 = t206 * pkin(4) + t205;
t236 = cos(qJ(6));
t232 = sin(qJ(6));
t218 = -qJD(5) + t221;
t195 = -t233 * t206 + t237 * t207;
t193 = qJD(6) + t194;
t192 = t236 * t195 - t232 * t218;
t191 = t232 * t195 + t236 * t218;
t186 = t194 * pkin(5) - t195 * pkin(11) + t196;
t185 = -t218 * pkin(11) + t251;
t184 = t218 * pkin(5) - t241;
t1 = [t239 / 0.2e1, 0, 0, t235 ^ 2 * t252 / 0.2e1, t235 * t245, t226 * t244, t226 * t243, t226 ^ 2 / 0.2e1, pkin(1) * t245 + t226 * t240, -pkin(1) * t235 * t252 - t226 * t249, -t202 * t243 + t209 * t215, t203 * t243 + t209 * t216, -t202 * t216 - t203 * t215, t203 ^ 2 / 0.2e1 + t202 ^ 2 / 0.2e1 + t209 ^ 2 / 0.2e1, t207 ^ 2 / 0.2e1, -t207 * t206, -t207 * t221, t206 * t221, t221 ^ 2 / 0.2e1, t205 * t206 - t221 * t242, t205 * t207 + t221 * t250, t195 ^ 2 / 0.2e1, -t195 * t194, -t195 * t218, t194 * t218, t218 ^ 2 / 0.2e1, t196 * t194 - t218 * t241, t196 * t195 + t218 * t251, t192 ^ 2 / 0.2e1, -t192 * t191, t192 * t193, -t191 * t193, t193 ^ 2 / 0.2e1 (-t232 * t185 + t236 * t186) * t193 + t184 * t191 -(t236 * t185 + t232 * t186) * t193 + t184 * t192;];
T_reg  = t1;
