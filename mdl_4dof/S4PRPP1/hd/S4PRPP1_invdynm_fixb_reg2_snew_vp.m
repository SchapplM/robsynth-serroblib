% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4PRPP1
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
%   pkin=[a2,a3,a4,d2,theta1]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4PRPP1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP1_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP1_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_invdynm_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:53:29
% EndTime: 2019-05-04 18:53:31
% DurationCPUTime: 1.12s
% Computational Cost: add. (1581->156), mult. (2642->134), div. (0->0), fcn. (1600->4), ass. (0->90)
t278 = g(3) - qJDD(1);
t335 = qJD(2) ^ 2;
t264 = pkin(3) * t335 - t278;
t284 = sin(qJ(2));
t285 = cos(qJ(2));
t306 = t285 * qJDD(2);
t256 = -t284 * t335 + t306;
t331 = pkin(4) * t256;
t235 = -pkin(3) * t306 + t284 * t264 - t331;
t281 = sin(pkin(5));
t282 = cos(pkin(5));
t307 = t284 * qJDD(2);
t257 = t285 * t335 + t307;
t253 = pkin(4) * t257;
t339 = -pkin(3) * t307 - t285 * t264 - t253;
t299 = t281 * t256 + t282 * t257;
t362 = qJ(1) * t299;
t365 = t281 * t235 + t282 * t339 - t362;
t239 = t284 * t278 + t331;
t337 = t285 * t278 - t253;
t364 = -t281 * t239 + t282 * t337 - t362;
t301 = t282 * t256 - t281 * t257;
t346 = qJ(1) * t301;
t363 = t282 * t239 + t281 * t337 + t346;
t258 = t281 * g(1) - t282 * g(2);
t259 = t282 * g(1) + t281 * g(2);
t296 = t285 * t258 + t284 * t259;
t310 = -t284 * t258 + t285 * t259;
t302 = -t284 * t296 - t285 * t310;
t203 = t284 * t310 - t285 * t296;
t319 = t282 * t203;
t359 = -t281 * t302 + t319;
t326 = t281 * t203;
t358 = t282 * t302 + t326;
t275 = qJDD(2) * qJ(4);
t279 = qJDD(2) * pkin(2);
t292 = -qJDD(3) + t296;
t288 = t335 * qJ(3) + t279 + t292;
t309 = (qJD(4) * qJD(2));
t214 = -t275 - t288 - (2 * t309);
t305 = (2 * qJD(3) * qJD(2)) - t310;
t308 = qJDD(2) * qJ(3);
t220 = -t335 * pkin(2) + t305 + t308;
t286 = qJDD(4) + t220;
t215 = -t335 * qJ(4) + t286;
t193 = t285 * t214 - t284 * t215;
t304 = t284 * t214 + t285 * t215;
t354 = t281 * t193 + t282 * t304;
t353 = t282 * t193 - t281 * t304;
t196 = t284 * t220 + t285 * t288;
t303 = t285 * t220 - t284 * t288;
t352 = -t281 * t196 + t282 * t303;
t351 = t282 * t196 + t281 * t303;
t347 = -t282 * t235 + t281 * t339 + t346;
t338 = 0.2e1 * t275 + (2 * t309);
t289 = -0.2e1 * t279 - t292;
t254 = pkin(1) * t256;
t255 = pkin(1) * t257;
t334 = pkin(1) * t278;
t333 = pkin(3) * t214;
t332 = pkin(3) * t215;
t328 = qJ(3) * t278;
t327 = qJ(4) * t214;
t320 = t281 * t278;
t313 = t282 * t278;
t311 = pkin(2) * t288 + qJ(3) * t220;
t297 = -t281 * t258 - t282 * t259;
t295 = 0.2e1 * t308 + t305;
t293 = -pkin(2) * t214 + qJ(3) * t215 - t327;
t291 = t282 * t258 - t281 * t259;
t290 = qJDD(4) + t295;
t287 = -t289 + t338;
t280 = pkin(3) * qJDD(2);
t218 = t310 - t255;
t217 = t254 + t296;
t212 = t289 - t254;
t209 = t255 + t295;
t208 = t255 + t290;
t207 = t328 + t333;
t206 = t254 + t287;
t205 = t332 + (pkin(2) + qJ(4)) * t278;
t200 = pkin(1) * t203;
t199 = pkin(4) * t302 + t334;
t190 = -pkin(4) * t196 + (-pkin(2) * t284 + qJ(3) * t285) * t278;
t189 = pkin(4) * t303 + (pkin(2) * t285 + qJ(3) * t284 + pkin(1)) * t278;
t188 = pkin(1) * t196 + t311;
t187 = pkin(4) * t193 - t284 * t205 + t285 * t207;
t186 = pkin(4) * t304 + t285 * t205 + t284 * t207 + t334;
t185 = -pkin(1) * t193 + t293;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t320, -t313, -t291, -qJ(1) * t291, 0, 0, t301, 0, -t299, 0, -t363, -t364, t359, pkin(4) * t319 + qJ(1) * t359 - t281 * t199, 0, -t301, t299, 0, 0, 0, -t351, t363, t364, -qJ(1) * t351 - t281 * t189 + t282 * t190, 0, t299, t301, 0, 0, 0, t353, t365, -t347, qJ(1) * t353 - t281 * t186 + t282 * t187; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t313, -t320, t297, qJ(1) * t297, 0, 0, t299, 0, t301, 0, t364, -t363, t358, pkin(4) * t326 + qJ(1) * t358 + t282 * t199, 0, -t299, -t301, 0, 0, 0, t352, -t364, t363, qJ(1) * t352 + t282 * t189 + t281 * t190, 0, -t301, t299, 0, 0, 0, t354, t347, t365, qJ(1) * t354 + t282 * t186 + t281 * t187; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t258, t259, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t217, t218, 0, -t200, qJDD(2), 0, 0, 0, 0, 0, 0, t212, t209, t188, qJDD(2), 0, 0, 0, 0, 0, 0, t208, t206, t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t278, -t258, 0, 0, 0, t256, 0, -t257, 0, -t239, -t337, t203, pkin(4) * t203, 0, -t256, t257, 0, 0, 0, -t196, t239, t337, t190, 0, t257, t256, 0, 0, 0, t193, t339, t235, t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t278, 0, -t259, 0, 0, 0, t257, 0, t256, 0, t337, -t239, t302, t199, 0, -t257, -t256, 0, 0, 0, t303, -t337, t239, t189, 0, -t256, t257, 0, 0, 0, t304, -t235, t339, t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t258, t259, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t217, t218, 0, -t200, qJDD(2), 0, 0, 0, 0, 0, 0, t212, t209, t188, qJDD(2), 0, 0, 0, 0, 0, 0, t208, t206, t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t335, 0, 0, -t278, -t296, 0, 0, -qJDD(2), t335, 0, 0, 0, -t288, 0, t278, t328, 0, t335, qJDD(2), 0, 0, 0, t214, -t264, -t280, t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t335, 0, qJDD(2), 0, t278, 0, -t310, 0, 0, -t335, -qJDD(2), 0, 0, 0, t220, -t278, 0, pkin(2) * t278, 0, -qJDD(2), t335, 0, 0, 0, t215, t280, -t264, t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t296, t310, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, t289, t295, t311, qJDD(2), 0, 0, 0, 0, 0, 0, t290, t287, t293; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, -t288, t220, 0, qJDD(2), 0, 0, 0, 0, 0, 0, t286, t288 + t338, -t327; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t335, 0, 0, 0, t288, 0, -t278, 0, 0, -t335, -qJDD(2), 0, 0, 0, -t214, t264, t280, -t333; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t335, qJDD(2), 0, 0, 0, -t220, t278, 0, 0, 0, qJDD(2), -t335, 0, 0, 0, -t215, -t280, t264, -qJ(4) * t278 - t332; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, t215, -t214, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t335, 0, 0, 0, -t215, 0, -t278, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t335, qJDD(2), 0, 0, 0, t214, t278, 0, 0;];
m_new_reg  = t1;
