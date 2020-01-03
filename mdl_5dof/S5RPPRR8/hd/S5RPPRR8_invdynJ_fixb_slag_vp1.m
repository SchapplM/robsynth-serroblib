% Calculate vector of inverse dynamics joint torques for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR8_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR8_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR8_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR8_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR8_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR8_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:04
% EndTime: 2019-12-31 18:01:13
% DurationCPUTime: 6.51s
% Computational Cost: add. (9698->409), mult. (11355->503), div. (0->0), fcn. (11488->8), ass. (0->218)
t302 = pkin(8) + qJ(4);
t276 = sin(t302);
t277 = cos(t302);
t357 = sin(qJ(1));
t358 = cos(qJ(1));
t141 = t358 * t276 - t357 * t277;
t205 = qJDD(1) - qJDD(4);
t206 = qJD(1) - qJD(4);
t338 = cos(pkin(8));
t190 = pkin(3) * t338 + pkin(2);
t207 = sin(pkin(8));
t289 = t357 * t207;
t226 = pkin(3) * t289 + t358 * t190;
t299 = t358 * pkin(2);
t120 = -t299 + t226;
t178 = t357 * t190;
t283 = qJD(1) * t358;
t274 = t207 * t283;
t181 = pkin(3) * t274;
t210 = qJD(1) ^ 2;
t197 = qJD(2) * t357;
t282 = qJD(1) * t357;
t310 = t358 * pkin(1) + t357 * qJ(2);
t313 = qJ(2) * t283 + t197;
t233 = qJDD(1) * t310 - qJDD(2) * t358 + (-pkin(1) * t282 + t197 + t313) * qJD(1);
t297 = t357 * pkin(2);
t216 = qJDD(1) * t299 - t210 * t297 + t233;
t214 = qJD(1) * (t181 + (t297 - t178) * qJD(1)) + qJDD(1) * t120 + t216;
t208 = sin(qJ(5));
t209 = cos(qJ(5));
t267 = rSges(6,1) * t208 + rSges(6,2) * t209;
t174 = -rSges(6,1) * t209 + rSges(6,2) * t208;
t155 = t174 * qJD(5);
t307 = qJD(5) * t155;
t140 = -t357 * t276 - t358 * t277;
t325 = t140 * t209;
t326 = t140 * t208;
t78 = -rSges(6,1) * t325 + rSges(6,2) * t326 + t141 * rSges(6,3);
t339 = -t140 * pkin(4) + pkin(7) * t141 + t78;
t101 = t206 * t140;
t102 = t206 * t141;
t305 = qJD(5) * t208;
t239 = t102 * t209 + t140 * t305;
t304 = qJD(5) * t209;
t240 = -t102 * t208 + t140 * t304;
t32 = t239 * rSges(6,1) + t240 * rSges(6,2) + t101 * rSges(6,3);
t368 = t102 * pkin(4) + t101 * pkin(7) + t32;
t80 = qJD(5) * t101 + qJDD(5) * t141;
t9 = -t141 * t307 + t339 * t205 + t368 * t206 + t80 * t267 + t214;
t401 = t9 - g(2);
t198 = qJD(2) * t358;
t142 = qJD(1) * t310 - t198;
t314 = qJD(1) * t198 + qJDD(2) * t357;
t236 = -t210 * t299 + t314;
t290 = t358 * t207;
t250 = pkin(3) * t290 - t178;
t119 = t297 + t250;
t200 = t358 * qJ(2);
t298 = t357 * pkin(1);
t171 = t298 - t200;
t273 = -t171 - t297;
t248 = t119 + t273;
t212 = (-qJD(1) * t120 - t142) * qJD(1) + t248 * qJDD(1) + t236;
t241 = -t101 * t209 + t141 * t305;
t242 = t101 * t208 + t141 * t304;
t31 = t241 * rSges(6,1) + t242 * rSges(6,2) + t102 * rSges(6,3);
t380 = t101 * pkin(4) - t102 * pkin(7) - t31;
t322 = t141 * t209;
t323 = t141 * t208;
t77 = -rSges(6,1) * t322 + rSges(6,2) * t323 - t140 * rSges(6,3);
t41 = -t141 * pkin(4) - t140 * pkin(7) + t77;
t81 = qJD(5) * t102 - qJDD(5) * t140;
t8 = -t140 * t307 - t41 * t205 + t380 * t206 - t81 * t267 + t212;
t400 = t8 - g(1);
t318 = t141 * rSges(5,1) - t140 * rSges(5,2);
t54 = -t101 * rSges(5,1) - t102 * rSges(5,2);
t19 = t205 * t318 - t206 * t54 + t212;
t399 = t19 - g(1);
t398 = t206 * t41;
t336 = Icges(6,4) * t209;
t255 = -Icges(6,2) * t208 + t336;
t67 = -Icges(6,6) * t141 + t140 * t255;
t340 = t208 * t67;
t337 = Icges(6,4) * t208;
t257 = Icges(6,1) * t209 - t337;
t71 = -Icges(6,5) * t141 + t140 * t257;
t258 = t209 * t71 - t340;
t68 = Icges(6,6) * t140 + t141 * t255;
t341 = t208 * t68;
t72 = Icges(6,5) * t140 + t141 * t257;
t396 = -t209 * t72 + t341;
t253 = Icges(6,5) * t209 - Icges(6,6) * t208;
t64 = Icges(6,3) * t140 + t141 * t253;
t351 = -t140 * t64 - t322 * t72;
t63 = -Icges(6,3) * t141 + t140 * t253;
t395 = -(t63 - t341) * t141 + t351;
t306 = qJD(5) * t267;
t394 = t141 * t306 + t206 * t339;
t259 = -t208 * t71 - t209 * t67;
t261 = -t208 * t72 - t209 * t68;
t308 = qJD(5) * t141;
t309 = qJD(5) * t140;
t359 = -t206 / 0.2e1;
t393 = ((-t259 * t140 + t261 * t141) * qJD(5) + t259 * t309 - t261 * t308) * t359;
t103 = rSges(5,1) * t140 + rSges(5,2) * t141;
t55 = t102 * rSges(5,1) - t101 * rSges(5,2);
t20 = -t103 * t205 + t206 * t55 + t214;
t392 = t20 - g(2);
t391 = t141 * t63;
t388 = t64 * t141;
t386 = t64 + t340;
t385 = t206 * t318;
t252 = Icges(6,5) * t208 + Icges(6,6) * t209;
t87 = t252 * t140;
t86 = t252 * t141;
t370 = qJD(1) * t171 - t197 + t313;
t254 = Icges(6,2) * t209 + t337;
t256 = Icges(6,1) * t208 + t336;
t251 = -t208 * t254 + t209 * t256;
t48 = t141 * t251 + t87;
t44 = t48 * t206;
t49 = t140 * t251 - t86;
t45 = t49 * t206;
t150 = -t338 * t358 - t289;
t271 = t357 * t338;
t151 = -t271 + t290;
t111 = -t150 * rSges(4,1) - t151 * rSges(4,2);
t272 = t310 + t299;
t376 = t111 + t272;
t375 = t151 * rSges(4,1) - t150 * rSges(4,2);
t275 = t310 + t226;
t311 = t358 * rSges(3,1) + t357 * rSges(3,3);
t371 = t310 + t311;
t366 = t140 * t306 - t398;
t365 = -qJD(1) * t275 + t198;
t364 = pkin(2) * t282 + t181 + t370 + (-t119 - t298 - t178) * qJD(1);
t345 = t254 * t141 - t72;
t347 = -t256 * t141 - t68;
t363 = t208 * t345 + t209 * t347;
t360 = m(4) + m(5);
t350 = t140 * t63 + t322 * t71;
t349 = t325 * t72 - t388;
t348 = t325 * t71 - t391;
t346 = -t256 * t140 - t67;
t344 = t254 * t140 - t71;
t321 = t253 * t206;
t143 = t150 * qJD(1);
t144 = -qJD(1) * t271 + t274;
t317 = t144 * rSges(4,1) - t143 * rSges(4,2);
t316 = -t254 + t257;
t315 = -t255 - t256;
t303 = m(6) + t360;
t301 = -t358 / 0.2e1;
t300 = t357 / 0.2e1;
t292 = t357 * rSges(3,1);
t281 = t309 / 0.2e1;
t280 = -t308 / 0.2e1;
t279 = t308 / 0.2e1;
t28 = Icges(6,4) * t239 + Icges(6,2) * t240 + Icges(6,6) * t101;
t30 = Icges(6,1) * t239 + Icges(6,4) * t240 + Icges(6,5) * t101;
t218 = qJD(5) * t259 + t208 * t28 - t209 * t30;
t27 = Icges(6,4) * t241 + Icges(6,2) * t242 + Icges(6,6) * t102;
t29 = Icges(6,1) * t241 + Icges(6,4) * t242 + Icges(6,5) * t102;
t219 = qJD(5) * t261 + t208 * t27 - t209 * t29;
t25 = Icges(6,5) * t241 + Icges(6,6) * t242 + Icges(6,3) * t102;
t26 = Icges(6,5) * t239 + Icges(6,6) * t240 + Icges(6,3) * t101;
t270 = -(-t101 * t396 - t102 * t64 - t140 * t25 + t141 * t219) * t140 + t141 * (t101 * t258 - t102 * t63 - t140 * t26 + t141 * t218);
t269 = -t140 * (-t101 * t64 + t102 * t396 + t140 * t219 + t141 * t25) + t141 * (-t101 * t63 - t102 * t258 + t140 * t218 + t141 * t26);
t268 = t143 * rSges(4,1) + t144 * rSges(4,2);
t15 = -t323 * t68 - t351;
t16 = -t323 * t67 + t350;
t266 = -t140 * t15 + t141 * t16;
t17 = -t326 * t68 + t349;
t18 = -t326 * t67 + t348;
t265 = -t140 * t17 + t141 * t18;
t225 = qJD(1) * t248 + t197;
t23 = t225 + t366;
t224 = (t120 + t272) * qJD(1) - t198;
t24 = t224 + t394;
t264 = -t140 * t23 - t141 * t24;
t263 = -t140 * t32 - t141 * t31;
t262 = t140 * t78 + t141 * t77;
t115 = -t208 * t256 - t209 * t254;
t249 = t375 + t273;
t235 = -t298 - t297;
t234 = -t292 - t298;
t177 = rSges(2,1) * t358 - rSges(2,2) * t357;
t173 = rSges(2,1) * t357 + rSges(2,2) * t358;
t228 = t208 * t344 + t209 * t346;
t221 = (t208 * t315 + t209 * t316) * t206;
t220 = t250 - t171;
t153 = t255 * qJD(5);
t154 = t257 * qJD(5);
t215 = qJD(5) * t115 - t153 * t208 + t154 * t209;
t11 = -t396 * qJD(5) - t208 * t29 - t209 * t27;
t12 = qJD(5) * t258 - t208 * t30 - t209 * t28;
t152 = t253 * qJD(5);
t13 = t101 * t251 - t102 * t252 + t140 * t152 + t141 * t215;
t14 = -t101 * t252 - t102 * t251 + t140 * t215 - t141 * t152;
t6 = qJD(5) * t266 + t44;
t7 = qJD(5) * t265 + t45;
t211 = t206 * (-t154 * t208 + t254 * t305 + (-qJD(5) * t256 - t153) * t209) + t279 * t6 - (-t259 + t49) * t80 / 0.2e1 - (-t261 + t48) * t81 / 0.2e1 + (t12 + t14) * t280 + (-Icges(5,3) + t115) * t205 + (t11 + t13 + t7) * t281;
t202 = t358 * rSges(3,3);
t195 = rSges(3,3) * t283;
t172 = t292 - t202;
t94 = t267 * t140;
t93 = t267 * t141;
t84 = qJD(1) * t249 + t197;
t61 = qJDD(1) * t311 + qJD(1) * (-rSges(3,1) * t282 + t195) + t233;
t60 = -qJD(1) * t142 - t210 * t311 + t314 + (-t171 - t172) * qJDD(1);
t47 = -t103 * t206 + t224;
t46 = t225 + t385;
t40 = qJD(1) * t317 + qJDD(1) * t111 + t216;
t39 = (-t142 + t268) * qJD(1) + t249 * qJDD(1) + t236;
t33 = qJD(5) * t262 - qJD(3);
t10 = (-t140 * t68 - t141 * t67) * t208 + t349 + t350;
t5 = qJD(5) * t263 - t77 * t80 + t78 * t81 + qJDD(3);
t1 = [-t211 + (-t44 + ((t17 + (-t258 + t64) * t141) * t141 + (t386 * t140 + t18 - t348 - t391) * t140) * qJD(5)) * t280 + (t45 + ((t396 * t141 + t15 + t348) * t141 + (-t141 * t386 - t10 + t16) * t140) * qJD(5)) * t281 - m(2) * (-g(1) * t173 + g(2) * t177) + t393 + (t401 * (t275 + t339) + (t365 + t380) * t23 + (-t267 * t309 + t23 + t364 + t368 + t398) * t24 + t400 * (t220 - t41)) * m(6) + (t392 * (-t103 + t275) + (t365 - t54) * t46 + (t364 + t46 + t55 - t385) * t47 + t399 * (t220 + t318)) * m(5) + ((t40 - g(2)) * t376 + (-qJD(1) * t272 + t198 + t268) * t84 + (t84 + t317 + (t297 + t235 - t375) * qJD(1) + t370) * (qJD(1) * t376 - t198) + (t39 - g(1)) * (t200 + t235 + t375)) * m(4) + ((-g(2) + t61) * t371 + (-g(1) + t60) * (t200 + t202 + t234) + (t195 + (t172 + t234) * qJD(1) + t370) * (qJD(1) * t371 - t198)) * m(3) + (m(2) * (t173 ^ 2 + t177 ^ 2) + Icges(2,3) + Icges(3,2) + Icges(4,3)) * qJDD(1); (-m(3) - t303) * (g(1) * t357 - g(2) * t358) + 0.2e1 * (t300 * t8 + t301 * t9) * m(6) + 0.2e1 * (t19 * t300 + t20 * t301) * m(5) + 0.2e1 * (t300 * t39 + t301 * t40) * m(4) + 0.2e1 * (t300 * t60 + t301 * t61) * m(3); m(6) * t5 + g(3) * t303 + qJDD(3) * t360; t211 + (-t45 + ((-t15 - t395) * t141 + (-t16 + (-t396 + t63) * t140 - t388) * t140) * qJD(5)) * t281 + (t44 + ((t10 - t17) * t141 + (t258 * t140 - t18 + t395) * t140) * qJD(5)) * t280 + t393 + (t400 * t41 + (t366 - t368) * t24 + (-t380 - t394) * t23 - t401 * t339) * m(6) + (t46 * t54 - t47 * t55 - (-t206 * t47 + t399) * t318 + (t206 * t46 + t392) * t103) * m(5); t101 * t7 / 0.2e1 + t141 * (qJD(5) * t269 + t14 * t206 + t17 * t81 + t18 * t80 + t205 * t49) / 0.2e1 + t80 * t265 / 0.2e1 + (t101 * t18 + t102 * t17 + t269) * t279 + t102 * t6 / 0.2e1 - t140 * (qJD(5) * t270 + t13 * t206 + t15 * t81 + t16 * t80 + t205 * t48) / 0.2e1 + t81 * t266 / 0.2e1 - (t101 * t16 + t102 * t15 + t270) * t309 / 0.2e1 + t205 * (t140 * t261 - t141 * t259) / 0.2e1 + t206 * (-t101 * t259 - t102 * t261 - t11 * t140 + t12 * t141) / 0.2e1 + ((t87 * t308 - t321) * t141 + (t221 + (-t363 * t140 + (-t86 + t228) * t141) * qJD(5)) * t140) * t280 + ((t86 * t309 + t321) * t140 + (t221 + (t228 * t141 + (-t363 - t87) * t140) * qJD(5)) * t141) * t281 + ((t208 * t316 - t209 * t315) * t206 + ((t140 * t345 - t141 * t344) * t209 + (-t140 * t347 + t141 * t346) * t208) * qJD(5)) * t359 + (-t5 * t262 + t33 * (t101 * t77 - t102 * t78 - t263) + t264 * t155 - (-t24 * t101 + t23 * t102 - t8 * t140 - t141 * t9) * t267 - (-t23 * t93 + t24 * t94) * t206 - (t33 * (t140 * t94 + t141 * t93) + t264 * t174) * qJD(5) - g(1) * t94 - g(2) * t93 - g(3) * t174) * m(6);];
tau = t1;
