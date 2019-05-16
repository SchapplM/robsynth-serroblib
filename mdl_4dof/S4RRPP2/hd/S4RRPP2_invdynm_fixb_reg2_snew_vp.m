% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4RRPP2
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:20
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4RRPP2_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP2_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP2_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_invdynm_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:20:32
% EndTime: 2019-05-04 19:20:33
% DurationCPUTime: 1.47s
% Computational Cost: add. (2531->165), mult. (3158->142), div. (0->0), fcn. (1616->4), ass. (0->96)
t326 = sin(qJ(2));
t322 = qJDD(1) + qJDD(2);
t328 = cos(qJ(2));
t357 = t328 * t322;
t323 = (qJD(1) + qJD(2));
t370 = t323 ^ 2;
t392 = -t326 * t370 + t357;
t401 = pkin(5) * t392;
t276 = t326 * g(3) + t401;
t360 = t326 * t322;
t373 = t328 * t370 + t360;
t395 = pkin(5) * t373;
t280 = -t328 * g(3) + t395;
t327 = sin(qJ(1));
t329 = cos(qJ(1));
t397 = -t327 * t373 + t329 * t392;
t405 = pkin(4) * t397;
t409 = -t329 * t276 + t327 * t280 - t405;
t398 = t327 * t392 + t329 * t373;
t404 = pkin(4) * t398;
t408 = t327 * t276 + t329 * t280 + t404;
t285 = pkin(1) * t392;
t303 = t329 * g(1) + t327 * g(2);
t331 = qJD(1) ^ 2;
t298 = -t331 * pkin(1) - t303;
t302 = t327 * g(1) - t329 * g(2);
t338 = qJDD(1) * pkin(1) + t302;
t274 = t326 * t298 - t328 * t338;
t275 = t328 * t298 + t326 * t338;
t349 = t326 * t274 + t328 * t275;
t245 = t328 * t274 - t326 * t275;
t356 = t329 * t245;
t400 = -t327 * t349 + t356;
t359 = t327 * t245;
t399 = t329 * t349 + t359;
t396 = pkin(1) * t373;
t353 = (2 * qJD(3) * t323) + t275;
t363 = t322 * qJ(3);
t345 = 0.2e1 * t363 + t353;
t391 = t345 + t396;
t258 = -t370 * pkin(2) + t353 + t363;
t384 = t370 * pkin(3);
t252 = t258 - t384;
t312 = t322 * pkin(2);
t342 = -qJDD(3) - t274;
t383 = qJ(3) * t370;
t336 = -t342 - t383;
t262 = -t312 + t336;
t311 = t322 * pkin(3);
t257 = t262 - t311;
t232 = t328 * t252 + t326 * t257;
t351 = t326 * t252 - t328 * t257;
t390 = t327 * t232 + t329 * t351;
t389 = t329 * t232 - t327 * t351;
t234 = t326 * t258 - t328 * t262;
t350 = t328 * t258 + t326 * t262;
t388 = -t327 * t234 + t329 * t350;
t387 = t329 * t234 + t327 * t350;
t325 = g(3) + qJDD(4);
t304 = qJ(4) * t370 + t325;
t264 = qJ(4) * t357 - t326 * t304 - t401;
t335 = qJ(4) * t360 + t328 * t304 - t395;
t386 = t329 * t264 - t327 * t335 - t405;
t385 = t327 * t264 + t329 * t335 - t404;
t369 = pkin(3) * t257;
t366 = qJ(4) * t252;
t365 = qJ(4) * t257;
t364 = qJ(4) * t322;
t354 = -pkin(2) * t262 + qJ(3) * t258;
t346 = -t327 * t302 - t329 * t303;
t344 = -pkin(2) * t257 + qJ(3) * t252 - t369;
t301 = t329 * qJDD(1) - t327 * t331;
t343 = -pkin(4) * t301 - t327 * g(3);
t339 = t329 * t302 - t327 * t303;
t310 = 0.2e1 * t312;
t337 = t310 + t342;
t333 = 0.2e1 * t311 - t336;
t332 = t310 + t333 - t383;
t330 = pkin(1) * g(3);
t300 = t327 * qJDD(1) + t329 * t331;
t284 = -pkin(4) * t300 + t329 * g(3);
t260 = -t274 + t285;
t259 = -t275 - t396;
t251 = t285 + t337;
t249 = qJ(3) * t325 - t365;
t247 = -t366 + (pkin(2) + pkin(3)) * t325;
t242 = pkin(1) * t245;
t241 = pkin(5) * t349 + t330;
t240 = t332 + t285;
t228 = -pkin(5) * t234 + (-pkin(2) * t326 + qJ(3) * t328) * g(3);
t227 = pkin(5) * t350 + t330 + (pkin(2) * t328 + qJ(3) * t326) * g(3);
t226 = pkin(1) * t234 + t354;
t225 = -pkin(5) * t351 - t326 * t247 + t328 * t249;
t224 = pkin(1) * t325 + pkin(5) * t232 + t328 * t247 + t326 * t249;
t223 = pkin(1) * t351 + t344;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t301, 0, -t300, 0, t343, -t284, -t339, -pkin(4) * t339, 0, 0, t397, 0, -t398, 0, t409, t408, t400, pkin(4) * t400 + pkin(5) * t356 - t327 * t241, 0, t397, 0, 0, t398, 0, t409, -t387, -t408, -pkin(4) * t387 - t327 * t227 + t329 * t228, 0, 0, -t397, 0, -t398, 0, t386, t385, t390, -pkin(4) * t390 - t327 * t224 + t329 * t225; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t300, 0, t301, 0, t284, t343, t346, pkin(4) * t346, 0, 0, t398, 0, t397, 0, -t408, t409, t399, pkin(4) * t399 + pkin(5) * t359 + t329 * t241, 0, t398, 0, 0, -t397, 0, -t408, t388, -t409, pkin(4) * t388 + t329 * t227 + t327 * t228, 0, 0, -t398, 0, t397, 0, t385, -t386, -t389, pkin(4) * t389 + t329 * t224 + t327 * t225; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t302, t303, 0, 0, 0, 0, 0, 0, 0, t322, t260, t259, 0, -t242, 0, 0, 0, t322, 0, 0, t251, 0, t391, t226, 0, 0, 0, 0, 0, t322, t240, t391, 0, t223; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t331, 0, 0, -g(3), -t302, 0, 0, 0, t392, 0, -t373, 0, -t276, t280, t245, pkin(5) * t245, 0, t392, 0, 0, t373, 0, -t276, -t234, -t280, t228, 0, 0, -t392, 0, -t373, 0, t264, t335, t351, t225; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t331, 0, qJDD(1), 0, g(3), 0, -t303, 0, 0, 0, t373, 0, t392, 0, -t280, -t276, t349, t241, 0, t373, 0, 0, -t392, 0, -t280, t350, t276, t227, 0, 0, -t373, 0, t392, 0, t335, -t264, -t232, t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t302, t303, 0, 0, 0, 0, 0, 0, 0, t322, t260, t259, 0, -t242, 0, 0, 0, t322, 0, 0, t251, 0, t391, t226, 0, 0, 0, 0, 0, t322, t240, t391, 0, t223; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t322, 0, -t370, 0, 0, -g(3), t274, 0, 0, t322, 0, 0, t370, 0, 0, t262, g(3), qJ(3) * g(3), 0, 0, -t322, 0, -t370, 0, t364, t304, -t257, t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t370, 0, t322, 0, g(3), 0, t275, 0, 0, t370, 0, 0, -t322, 0, g(3), t258, 0, pkin(2) * g(3), 0, 0, -t370, 0, t322, 0, t304, -t364, -t252, t247; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t322, -t274, -t275, 0, 0, 0, 0, 0, t322, 0, 0, t337, 0, t345, t354, 0, 0, 0, 0, 0, t322, t332, t345, 0, t344; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t322, 0, 0, t370, 0, 0, t262, g(3), 0, 0, 0, -t322, 0, -t370, 0, t364, t304, -t257, -t365; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t322, 0, 0, -t262, 0, t258, 0, 0, 0, 0, 0, 0, t322, t312 + t333, t252 + t384, 0, -t369; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t370, 0, 0, t322, 0, -g(3), -t258, 0, 0, 0, 0, t370, 0, -t322, 0, -t304, t364, t252, -pkin(3) * t325 + t366; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t322, 0, -t370, 0, 0, t325, -t257, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t370, 0, -t322, 0, -t325, 0, t252, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t322, t257, -t252, 0, 0;];
m_new_reg  = t1;
