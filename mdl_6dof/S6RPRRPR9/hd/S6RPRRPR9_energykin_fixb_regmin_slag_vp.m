% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRPR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:31:33
% EndTime: 2019-03-09 05:31:33
% DurationCPUTime: 0.27s
% Computational Cost: add. (1117->69), mult. (3606->140), div. (0->0), fcn. (3027->14), ass. (0->58)
t288 = sin(pkin(12));
t292 = cos(pkin(12));
t290 = sin(pkin(6));
t307 = qJD(1) * t290;
t302 = t292 * t307;
t294 = cos(pkin(6));
t306 = qJD(1) * t294;
t305 = pkin(1) * t306;
t279 = qJ(2) * t302 + t288 * t305;
t289 = sin(pkin(7));
t293 = cos(pkin(7));
t310 = t290 * t292;
t267 = (t289 * t294 + t293 * t310) * qJD(1) * pkin(9) + t279;
t284 = t292 * t305;
t312 = t288 * t290;
t269 = t284 + (pkin(2) * t294 + (-pkin(9) * t293 - qJ(2)) * t312) * qJD(1);
t274 = qJD(2) + (-pkin(9) * t288 * t289 - pkin(2) * t292 - pkin(1)) * t307;
t297 = sin(qJ(3));
t299 = cos(qJ(3));
t314 = -t297 * t267 + (t269 * t293 + t274 * t289) * t299;
t313 = cos(qJ(4));
t311 = t289 * t297;
t309 = t293 * t297;
t271 = (t294 * t311 + (t288 * t299 + t292 * t309) * t290) * qJD(1);
t276 = t289 * t302 - t293 * t306 - qJD(3);
t296 = sin(qJ(4));
t262 = t313 * t271 - t296 * t276;
t303 = t288 * t307;
t270 = t297 * t303 + (-t289 * t306 - t293 * t302) * t299;
t268 = qJD(4) + t270;
t259 = -t289 * t269 + t293 * t274;
t252 = t270 * pkin(3) - t271 * pkin(10) + t259;
t304 = t299 * t267 + t269 * t309 + t274 * t311;
t255 = -t276 * pkin(10) + t304;
t301 = t313 * t252 - t296 * t255;
t244 = t268 * pkin(4) - t262 * qJ(5) + t301;
t261 = t296 * t271 + t313 * t276;
t308 = t296 * t252 + t313 * t255;
t246 = -t261 * qJ(5) + t308;
t287 = sin(pkin(13));
t291 = cos(pkin(13));
t241 = t287 * t244 + t291 * t246;
t257 = -t291 * t261 - t287 * t262;
t240 = t291 * t244 - t287 * t246;
t254 = t276 * pkin(3) - t314;
t247 = t261 * pkin(4) + qJD(5) + t254;
t298 = cos(qJ(6));
t295 = sin(qJ(6));
t285 = -pkin(1) * t307 + qJD(2);
t278 = -qJ(2) * t303 + t284;
t258 = -t287 * t261 + t291 * t262;
t256 = qJD(6) - t257;
t249 = t298 * t258 + t295 * t268;
t248 = t295 * t258 - t298 * t268;
t242 = -t257 * pkin(5) - t258 * pkin(11) + t247;
t239 = t268 * pkin(11) + t241;
t238 = -t268 * pkin(5) - t240;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0 (t278 * t294 - t285 * t310) * qJD(1) (-t279 * t294 + t285 * t312) * qJD(1) (-t278 * t288 + t279 * t292) * t307, t279 ^ 2 / 0.2e1 + t278 ^ 2 / 0.2e1 + t285 ^ 2 / 0.2e1, t271 ^ 2 / 0.2e1, -t271 * t270, -t271 * t276, t270 * t276, t276 ^ 2 / 0.2e1, t259 * t270 - t276 * t314, t259 * t271 + t304 * t276, t262 ^ 2 / 0.2e1, -t262 * t261, t262 * t268, -t261 * t268, t268 ^ 2 / 0.2e1, t254 * t261 + t301 * t268, t254 * t262 - t308 * t268, -t240 * t258 + t241 * t257, t241 ^ 2 / 0.2e1 + t240 ^ 2 / 0.2e1 + t247 ^ 2 / 0.2e1, t249 ^ 2 / 0.2e1, -t249 * t248, t249 * t256, -t248 * t256, t256 ^ 2 / 0.2e1 (-t295 * t239 + t298 * t242) * t256 + t238 * t248 -(t298 * t239 + t295 * t242) * t256 + t238 * t249;];
T_reg  = t1;
