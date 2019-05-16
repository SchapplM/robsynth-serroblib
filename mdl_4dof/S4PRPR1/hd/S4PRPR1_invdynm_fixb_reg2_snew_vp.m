% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4PRPR1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:57
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4PRPR1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR1_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_invdynm_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:57:15
% EndTime: 2019-05-04 18:57:17
% DurationCPUTime: 1.54s
% Computational Cost: add. (3642->165), mult. (5803->164), div. (0->0), fcn. (3576->6), ass. (0->99)
t311 = -qJD(2) + qJD(4);
t308 = t311 ^ 2;
t309 = qJDD(2) - qJDD(4);
t317 = sin(qJ(4));
t319 = cos(qJ(4));
t286 = t319 * t308 - t317 * t309;
t318 = sin(qJ(2));
t320 = cos(qJ(2));
t327 = t317 * t308 + t319 * t309;
t259 = t318 * t286 - t320 * t327;
t313 = g(3) - qJDD(1);
t273 = pkin(5) * t286 + t319 * t313;
t367 = pkin(5) * t327 + t317 * t313;
t226 = pkin(4) * t259 - t318 * t273 + t320 * t367;
t315 = sin(pkin(6));
t316 = cos(pkin(6));
t257 = t320 * t286 + t318 * t327;
t378 = t315 * t257 + t316 * t259;
t390 = -pkin(4) * t257 + t320 * t273 + t318 * t367;
t398 = qJ(1) * t378 + t316 * t226 - t315 * t390;
t391 = t316 * t257 - t315 * t259;
t397 = -qJ(1) * t391 + t315 * t226 + t316 * t390;
t321 = qJD(2) ^ 2;
t298 = t318 * qJDD(2) + t320 * t321;
t299 = t320 * qJDD(2) - t318 * t321;
t265 = t316 * t298 + t315 * t299;
t275 = pkin(4) * t299 + t318 * t313;
t340 = -pkin(4) * t298 + t320 * t313;
t236 = -qJ(1) * t265 - t315 * t275 + t316 * t340;
t307 = 2 * qJD(3) * qJD(2);
t300 = t315 * g(1) - t316 * g(2);
t301 = t316 * g(1) + t315 * g(2);
t343 = -t318 * t300 + t320 * t301;
t341 = t307 - t343;
t342 = qJDD(2) * qJ(3);
t325 = t341 + t342;
t364 = pkin(3) + pkin(2);
t250 = -t321 * t364 + t325;
t314 = qJDD(2) * pkin(2);
t333 = t320 * t300 + t318 * t301;
t330 = -qJDD(3) + t333;
t263 = -t321 * qJ(3) - t314 - t330;
t254 = -qJDD(2) * pkin(3) + t263;
t228 = t317 * t250 - t319 * t254;
t229 = t319 * t250 + t317 * t254;
t217 = t319 * t228 - t317 * t229;
t329 = t317 * t228 + t319 * t229;
t207 = t318 * t217 - t320 * t329;
t381 = t320 * t217 + t318 * t329;
t393 = -t315 * t207 + t316 * t381;
t392 = t316 * t207 + t315 * t381;
t338 = -t318 * t333 - t320 * t343;
t242 = t318 * t343 - t320 * t333;
t351 = t316 * t242;
t387 = -t315 * t338 + t351;
t356 = t315 * t242;
t386 = t316 * t338 + t356;
t336 = -t315 * t298 + t316 * t299;
t365 = qJ(1) * t336 + t316 * t275 + t315 * t340;
t255 = -t321 * pkin(2) + t325;
t233 = t318 * t255 - t320 * t263;
t339 = t320 * t255 + t318 * t263;
t380 = -t315 * t233 + t316 * t339;
t379 = t316 * t233 + t315 * t339;
t363 = pkin(1) * t313;
t362 = pkin(3) * t217;
t360 = pkin(5) * t217;
t359 = pkin(5) * t329;
t357 = qJ(3) * t313;
t352 = t315 * t313;
t347 = t316 * t313;
t344 = -pkin(2) * t263 + qJ(3) * t255;
t334 = -t315 * t300 - t316 * t301;
t252 = -pkin(1) * t298 + t343;
t332 = pkin(2) * t217 + qJ(3) * t329 + t362;
t331 = pkin(3) * t286 + t229;
t328 = t316 * t300 - t315 * t301;
t326 = pkin(3) * t327 + t228;
t324 = 0.2e1 * t314 + t330;
t323 = pkin(2) * t286 + qJ(3) * t327 + t331;
t322 = pkin(2) * t327 - qJ(3) * t286 + t326;
t306 = 0.2e1 * t342;
t297 = pkin(1) * t299;
t251 = t297 + t333;
t245 = t297 + t324;
t244 = -t252 + t306 + t307;
t239 = pkin(1) * t242;
t238 = pkin(4) * t338 + t363;
t223 = -pkin(4) * t233 + (-pkin(2) * t318 + qJ(3) * t320) * t313;
t222 = pkin(4) * t339 + (pkin(2) * t320 + qJ(3) * t318 + pkin(1)) * t313;
t221 = pkin(1) * t233 + t344;
t220 = -pkin(1) * t259 + t322;
t219 = pkin(1) * t257 + t323;
t212 = t357 + t360;
t211 = t313 * t364 - t359;
t206 = -pkin(4) * t381 - t318 * t211 + t320 * t212;
t205 = -pkin(4) * t207 + t320 * t211 + t318 * t212 + t363;
t204 = pkin(1) * t381 + t332;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t352, -t347, -t328, -qJ(1) * t328, 0, 0, t336, 0, -t265, 0, -t365, -t236, t387, pkin(4) * t351 + qJ(1) * t387 - t315 * t238, 0, t336, 0, 0, t265, 0, -t365, -t379, t236, -qJ(1) * t379 - t315 * t222 + t316 * t223, 0, 0, t378, 0, -t391, 0, t398, t397, t393, -qJ(1) * t393 - t315 * t205 + t316 * t206; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t347, -t352, t334, qJ(1) * t334, 0, 0, t265, 0, t336, 0, t236, -t365, t386, pkin(4) * t356 + qJ(1) * t386 + t316 * t238, 0, t265, 0, 0, -t336, 0, t236, t380, t365, qJ(1) * t380 + t316 * t222 + t315 * t223, 0, 0, -t391, 0, -t378, 0, t397, -t398, t392, -qJ(1) * t392 + t316 * t205 + t315 * t206; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t300, t301, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t251, t252, 0, -t239, 0, 0, 0, qJDD(2), 0, 0, t245, 0, t244, t221, 0, 0, 0, 0, 0, t309, t220, t219, 0, t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t313, -t300, 0, 0, 0, t299, 0, -t298, 0, -t275, -t340, t242, pkin(4) * t242, 0, t299, 0, 0, t298, 0, -t275, -t233, t340, t223, 0, 0, t259, 0, -t257, 0, t226, t390, t381, t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t313, 0, -t301, 0, 0, 0, t298, 0, t299, 0, t340, -t275, t338, t238, 0, t298, 0, 0, -t299, 0, t340, t339, t275, t222, 0, 0, -t257, 0, -t259, 0, t390, -t226, t207, t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t300, t301, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t251, t252, 0, -t239, 0, 0, 0, qJDD(2), 0, 0, t245, 0, t244, t221, 0, 0, 0, 0, 0, t309, t220, t219, 0, t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t321, 0, 0, -t313, -t333, 0, 0, qJDD(2), 0, 0, t321, 0, 0, t263, t313, t357, 0, 0, -t327, 0, -t286, 0, t367, t273, t217, t212; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t321, 0, qJDD(2), 0, t313, 0, -t343, 0, 0, t321, 0, 0, -qJDD(2), 0, t313, t255, 0, pkin(2) * t313, 0, 0, -t286, 0, t327, 0, t273, -t367, -t329, t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t333, t343, 0, 0, 0, 0, 0, qJDD(2), 0, 0, t324, 0, t306 + t341, t344, 0, 0, 0, 0, 0, t309, t322, t323, 0, t332; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, 0, t321, 0, 0, t263, t313, 0, 0, 0, -t327, 0, -t286, 0, t367, t273, t217, t360; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, 0, -t263, 0, t255, 0, 0, 0, 0, 0, 0, t309, t326, t331, 0, t362; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t321, 0, 0, qJDD(2), 0, -t313, -t255, 0, 0, 0, 0, t286, 0, -t327, 0, -t273, t367, t329, -pkin(3) * t313 + t359; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t309, 0, -t308, 0, 0, t313, t228, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t308, 0, -t309, 0, -t313, 0, t229, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t309, -t228, -t229, 0, 0;];
m_new_reg  = t1;
