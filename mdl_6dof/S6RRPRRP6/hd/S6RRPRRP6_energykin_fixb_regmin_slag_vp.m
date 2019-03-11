% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:14:09
% EndTime: 2019-03-09 12:14:09
% DurationCPUTime: 0.17s
% Computational Cost: add. (694->54), mult. (1890->110), div. (0->0), fcn. (1493->10), ass. (0->47)
t226 = sin(pkin(6));
t235 = qJD(1) ^ 2;
t247 = t226 ^ 2 * t235;
t225 = sin(pkin(11));
t227 = cos(pkin(11));
t234 = cos(qJ(2));
t243 = qJD(1) * t226;
t238 = t234 * t243;
t231 = sin(qJ(2));
t239 = t231 * t243;
t215 = -t225 * t239 + t227 * t238;
t214 = qJD(4) - t215;
t242 = cos(pkin(6)) * qJD(1);
t241 = pkin(1) * t242;
t222 = t234 * t241;
t223 = qJD(2) + t242;
t209 = t223 * pkin(2) + t222 + (-pkin(8) - qJ(3)) * t239;
t244 = pkin(8) * t238 + t231 * t241;
t212 = qJ(3) * t238 + t244;
t200 = t225 * t209 + t227 * t212;
t198 = pkin(9) * t223 + t200;
t216 = (t225 * t234 + t227 * t231) * t243;
t217 = qJD(3) + (-pkin(2) * t234 - pkin(1)) * t243;
t204 = -pkin(3) * t215 - pkin(9) * t216 + t217;
t230 = sin(qJ(4));
t233 = cos(qJ(4));
t245 = t233 * t198 + t230 * t204;
t192 = pkin(10) * t214 + t245;
t199 = t209 * t227 - t225 * t212;
t197 = -pkin(3) * t223 - t199;
t207 = t216 * t230 - t233 * t223;
t208 = t216 * t233 + t223 * t230;
t194 = pkin(4) * t207 - pkin(10) * t208 + t197;
t229 = sin(qJ(5));
t232 = cos(qJ(5));
t246 = t232 * t192 + t229 * t194;
t240 = t234 * t247;
t237 = -t230 * t198 + t204 * t233;
t236 = -t229 * t192 + t194 * t232;
t191 = -pkin(4) * t214 - t237;
t206 = qJD(5) + t207;
t202 = t208 * t232 + t214 * t229;
t201 = t208 * t229 - t232 * t214;
t189 = pkin(5) * t201 - qJ(6) * t202 + t191;
t188 = qJ(6) * t206 + t246;
t187 = -pkin(5) * t206 + qJD(6) - t236;
t1 = [t235 / 0.2e1, 0, 0, t231 ^ 2 * t247 / 0.2e1, t231 * t240, t223 * t239, t223 * t238, t223 ^ 2 / 0.2e1, pkin(1) * t240 + (-pkin(8) * t239 + t222) * t223, -pkin(1) * t231 * t247 - t244 * t223, -t199 * t216 + t200 * t215, t200 ^ 2 / 0.2e1 + t199 ^ 2 / 0.2e1 + t217 ^ 2 / 0.2e1, t208 ^ 2 / 0.2e1, -t208 * t207, t208 * t214, -t207 * t214, t214 ^ 2 / 0.2e1, t197 * t207 + t237 * t214, t197 * t208 - t245 * t214, t202 ^ 2 / 0.2e1, -t202 * t201, t202 * t206, -t201 * t206, t206 ^ 2 / 0.2e1, t191 * t201 + t236 * t206, t191 * t202 - t246 * t206, -t187 * t206 + t189 * t201, t187 * t202 - t188 * t201, t188 * t206 - t189 * t202, t188 ^ 2 / 0.2e1 + t189 ^ 2 / 0.2e1 + t187 ^ 2 / 0.2e1;];
T_reg  = t1;
