% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_energykin_fixb_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:25:03
% EndTime: 2019-03-09 01:25:03
% DurationCPUTime: 0.23s
% Computational Cost: add. (739->58), mult. (1945->127), div. (0->0), fcn. (1668->16), ass. (0->57)
t284 = cos(pkin(7));
t278 = t284 * qJD(2) + qJD(3);
t280 = sin(pkin(8));
t283 = cos(pkin(8));
t294 = cos(qJ(3));
t281 = sin(pkin(7));
t306 = qJD(2) * t281;
t301 = t294 * t306;
t318 = t278 * t280 + t283 * t301;
t295 = cos(qJ(2));
t308 = qJD(1) * sin(pkin(6));
t273 = qJD(2) * pkin(2) + t295 * t308;
t307 = qJD(1) * cos(pkin(6));
t317 = t273 * t284 + t281 * t307;
t290 = sin(qJ(2));
t270 = pkin(10) * t306 + t290 * t308;
t289 = sin(qJ(3));
t304 = t294 * t270 + t317 * t289;
t252 = t318 * pkin(11) + t304;
t309 = t317 * t294;
t255 = t278 * pkin(3) + (-pkin(11) * t283 * t306 - t270) * t289 + t309;
t277 = t284 * t307;
t258 = t277 + (-t273 + (-pkin(11) * t280 * t289 - pkin(3) * t294) * qJD(2)) * t281;
t288 = sin(qJ(4));
t293 = cos(qJ(4));
t316 = -t288 * t252 + (t255 * t283 + t258 * t280) * t293;
t296 = qJD(2) ^ 2;
t313 = t281 ^ 2 * t296;
t312 = t280 * t288;
t311 = t283 * t288;
t265 = -t283 * t278 + t280 * t301 - qJD(4);
t305 = t293 * t252 + t255 * t311 + t258 * t312;
t242 = -t265 * pkin(12) + t305;
t245 = -t280 * t255 + t283 * t258;
t302 = t289 * t306;
t261 = t288 * t302 - t318 * t293;
t262 = t278 * t312 + (t289 * t293 + t294 * t311) * t306;
t244 = t261 * pkin(4) - t262 * pkin(12) + t245;
t287 = sin(qJ(5));
t292 = cos(qJ(5));
t310 = t292 * t242 + t287 * t244;
t300 = qJD(2) * t308;
t253 = t287 * t262 + t292 * t265;
t298 = -t287 * t242 + t292 * t244;
t241 = t265 * pkin(4) - t316;
t291 = cos(qJ(6));
t286 = sin(qJ(6));
t264 = -t281 * t273 + t277;
t260 = qJD(5) + t261;
t254 = t292 * t262 - t287 * t265;
t251 = qJD(6) + t253;
t247 = t291 * t254 + t286 * t260;
t246 = t286 * t254 - t291 * t260;
t239 = t253 * pkin(5) - t254 * pkin(13) + t241;
t238 = t260 * pkin(13) + t310;
t237 = -t260 * pkin(5) - t298;
t1 = [qJD(1) ^ 2 / 0.2e1, t296 / 0.2e1, t295 * t300, -t290 * t300, t289 ^ 2 * t313 / 0.2e1, t289 * t294 * t313, t278 * t302, t278 * t301, t278 ^ 2 / 0.2e1, -t264 * t301 + (-t289 * t270 + t309) * t278, t264 * t302 - t304 * t278, t262 ^ 2 / 0.2e1, -t261 * t262, -t265 * t262, t265 * t261, t265 ^ 2 / 0.2e1, t245 * t261 - t265 * t316, t245 * t262 + t305 * t265, t254 ^ 2 / 0.2e1, -t254 * t253, t254 * t260, -t253 * t260, t260 ^ 2 / 0.2e1, t241 * t253 + t298 * t260, t241 * t254 - t310 * t260, t247 ^ 2 / 0.2e1, -t247 * t246, t247 * t251, -t246 * t251, t251 ^ 2 / 0.2e1 (-t286 * t238 + t291 * t239) * t251 + t237 * t246 -(t291 * t238 + t286 * t239) * t251 + t237 * t247;];
T_reg  = t1;
