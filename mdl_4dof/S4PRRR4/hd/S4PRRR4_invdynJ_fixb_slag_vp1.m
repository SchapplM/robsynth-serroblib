% Calculate vector of inverse dynamics joint torques for
% S4PRRR4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR4_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR4_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR4_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR4_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:30
% EndTime: 2019-12-31 16:32:38
% DurationCPUTime: 6.68s
% Computational Cost: add. (9868->517), mult. (8958->701), div. (0->0), fcn. (7099->6), ass. (0->291)
t222 = pkin(7) + qJ(2);
t216 = sin(t222);
t224 = qJ(3) + qJ(4);
t219 = cos(t224);
t341 = t216 * t219;
t218 = sin(t224);
t342 = t216 * t218;
t217 = cos(t222);
t354 = Icges(5,6) * t217;
t108 = Icges(5,4) * t341 - Icges(5,2) * t342 - t354;
t214 = Icges(5,4) * t219;
t162 = Icges(5,1) * t218 + t214;
t403 = -t162 * t216 - t108;
t260 = -Icges(5,2) * t218 + t214;
t109 = Icges(5,6) * t216 + t217 * t260;
t402 = -t162 * t217 - t109;
t401 = t162 + t260;
t173 = Icges(5,4) * t342;
t358 = Icges(5,5) * t217;
t110 = Icges(5,1) * t341 - t173 - t358;
t351 = t108 * t218;
t259 = -t110 * t219 + t351;
t248 = t259 * t216;
t159 = Icges(5,5) * t219 - Icges(5,6) * t218;
t107 = Icges(5,3) * t216 + t159 * t217;
t360 = Icges(5,4) * t218;
t163 = Icges(5,1) * t219 - t360;
t111 = Icges(5,5) * t216 + t163 * t217;
t336 = t217 * t219;
t374 = t216 * t107 + t111 * t336;
t400 = -t248 - t374;
t176 = rSges(5,2) * t342;
t112 = rSges(5,1) * t341 - t217 * rSges(5,3) - t176;
t207 = t216 * rSges(5,3);
t337 = t217 * t218;
t113 = rSges(5,1) * t336 - rSges(5,2) * t337 + t207;
t372 = rSges(5,2) * t219;
t164 = rSges(5,1) * t218 + t372;
t130 = t164 * t216;
t131 = t164 * t217;
t223 = qJD(3) + qJD(4);
t152 = t216 * t223;
t153 = t217 * t223;
t212 = t217 * pkin(5);
t156 = pkin(2) * t216 - t212;
t227 = -pkin(6) - pkin(5);
t198 = t217 * t227;
t226 = cos(qJ(3));
t221 = t226 * pkin(3);
t376 = t221 + pkin(2);
t317 = -t216 * t376 - t198;
t100 = t156 + t317;
t157 = t217 * pkin(2) + t216 * pkin(5);
t276 = -t216 * t227 + t217 * t376;
t101 = t276 - t157;
t35 = t112 * t152 + t113 * t153 + qJD(1) + (-t100 * t216 + t101 * t217) * qJD(3);
t215 = t219 * rSges(5,1);
t398 = -t218 * rSges(5,2) + t215;
t225 = sin(qJ(3));
t305 = qJD(3) * t225;
t293 = t217 * t305;
t274 = pkin(3) * t293;
t244 = -t153 * t164 - t274;
t329 = t100 - t112;
t296 = -t156 + t329;
t41 = qJD(2) * t296 + t244;
t298 = pkin(3) * t305;
t179 = t216 * t298;
t328 = t101 + t113;
t42 = -t152 * t164 - t179 + (t157 + t328) * qJD(2);
t302 = qJD(2) * qJD(3);
t148 = qJDD(3) * t216 + t217 * t302;
t301 = qJD(2) * qJD(4);
t104 = qJDD(4) * t216 + t217 * t301 + t148;
t194 = t216 * t302;
t105 = t216 * t301 + t194 + (-qJDD(3) - qJDD(4)) * t217;
t149 = -qJDD(3) * t217 + t194;
t299 = t223 * t372;
t310 = qJD(2) * t216;
t309 = qJD(2) * t217;
t318 = rSges(5,3) * t309 + qJD(2) * t176;
t333 = t218 * t223;
t65 = -t217 * t299 + (-t217 * t333 - t219 * t310) * rSges(5,1) + t318;
t251 = t164 * t223;
t66 = -t216 * t251 + (t217 * t398 + t207) * qJD(2);
t199 = pkin(5) * t309;
t82 = -t274 - t199 + (-t216 * t221 - t198) * qJD(2);
t83 = -t179 + (t217 * t221 + (-pkin(5) - t227) * t216) * qJD(2);
t8 = -t100 * t148 - t101 * t149 + t104 * t112 - t105 * t113 + t152 * t66 + t153 * t65 + qJDD(1) + (t216 * t83 + t217 * t82) * qJD(3);
t399 = -t41 * (qJD(2) * t130 - t153 * t398) - t35 * (-t152 * t130 - t131 * t153) - t42 * (-qJD(2) * t131 - t152 * t398) + t8 * (t216 * t112 + t217 * t113);
t220 = Icges(4,4) * t226;
t261 = -Icges(4,2) * t225 + t220;
t185 = Icges(4,1) * t225 + t220;
t182 = Icges(4,5) * t226 - Icges(4,6) * t225;
t181 = Icges(4,5) * t225 + Icges(4,6) * t226;
t245 = qJD(3) * t181;
t361 = Icges(4,4) * t225;
t186 = Icges(4,1) * t226 - t361;
t120 = Icges(4,5) * t216 + t186 * t217;
t118 = Icges(4,6) * t216 + t217 * t261;
t348 = t118 * t225;
t256 = -t120 * t226 + t348;
t353 = Icges(4,3) * t217;
t397 = -t217 * t245 + (-t182 * t216 + t256 + t353) * qJD(2);
t340 = t216 * t225;
t191 = Icges(4,4) * t340;
t339 = t216 * t226;
t359 = Icges(4,5) * t217;
t119 = Icges(4,1) * t339 - t191 - t359;
t355 = Icges(4,6) * t217;
t117 = Icges(4,4) * t339 - Icges(4,2) * t340 - t355;
t349 = t117 * t225;
t257 = -t119 * t226 + t349;
t116 = Icges(4,3) * t216 + t182 * t217;
t311 = qJD(2) * t116;
t396 = qJD(2) * t257 - t216 * t245 + t311;
t158 = Icges(5,5) * t218 + Icges(5,6) * t219;
t249 = t158 * t223;
t350 = t109 * t218;
t352 = Icges(5,3) * t217;
t395 = -t217 * t249 + (-t111 * t219 - t159 * t216 + t350 + t352) * qJD(2);
t312 = qJD(2) * t107;
t394 = qJD(2) * t259 - t216 * t249 + t312;
t115 = Icges(4,5) * t339 - Icges(4,6) * t340 - t353;
t47 = -t217 * t115 - t216 * t257;
t160 = Icges(5,2) * t219 + t360;
t254 = t160 * t218 - t162 * t219;
t393 = qJD(2) * t254 + t159 * t223;
t183 = Icges(4,2) * t226 + t361;
t252 = t225 * t183 - t185 * t226;
t392 = qJD(2) * t252 + t182 * qJD(3);
t323 = -Icges(4,2) * t339 + t119 - t191;
t325 = t185 * t216 + t117;
t391 = -t225 * t323 - t226 * t325;
t390 = qJD(2) * t401 + t152 * (-t160 * t217 + t111) - t153 * (-Icges(5,2) * t341 + t110 - t173);
t389 = m(2) + m(3);
t388 = t104 / 0.2e1;
t387 = t105 / 0.2e1;
t386 = t148 / 0.2e1;
t385 = t149 / 0.2e1;
t384 = -t152 / 0.2e1;
t383 = t152 / 0.2e1;
t382 = -t153 / 0.2e1;
t381 = t153 / 0.2e1;
t380 = t216 / 0.2e1;
t379 = -t217 / 0.2e1;
t378 = -qJD(2) / 0.2e1;
t377 = qJD(2) / 0.2e1;
t106 = Icges(5,5) * t341 - Icges(5,6) * t342 - t352;
t375 = -t216 * t106 - t110 * t336;
t373 = rSges(4,1) * t226;
t187 = rSges(4,1) * t225 + rSges(4,2) * t226;
t143 = t187 * t217;
t307 = qJD(3) * t216;
t208 = t216 * rSges(4,3);
t334 = t217 * t226;
t335 = t217 * t225;
t122 = rSges(4,1) * t334 - rSges(4,2) * t335 + t208;
t90 = t122 + t157;
t68 = qJD(2) * t90 - t187 * t307;
t371 = t143 * t68;
t306 = qJD(3) * t217;
t294 = t187 * t306;
t313 = rSges(4,2) * t340 + t217 * rSges(4,3);
t121 = rSges(4,1) * t339 - t313;
t321 = -t121 - t156;
t67 = qJD(2) * t321 - t294;
t370 = t216 * t67;
t368 = t42 * t217;
t367 = t67 * t217;
t366 = qJDD(2) / 0.2e1;
t365 = -t216 * t115 - t119 * t334;
t364 = t216 * t116 + t120 * t334;
t347 = t158 * t216;
t346 = t158 * t217;
t345 = t160 * t223;
t344 = t181 * t216;
t343 = t181 * t217;
t332 = t226 * qJD(3) ^ 2;
t57 = -t216 * t254 - t346;
t331 = t57 * qJD(2);
t78 = -t216 * t252 - t343;
t330 = t78 * qJD(2);
t324 = -t185 * t217 - t118;
t322 = -t183 * t217 + t120;
t320 = qJD(2) * (-pkin(2) * t310 + t199) + qJDD(2) * t157;
t308 = qJD(2) * t225;
t316 = rSges(4,2) * t216 * t308 + rSges(4,3) * t309;
t315 = -t183 + t186;
t314 = t185 + t261;
t304 = qJD(3) * t226;
t303 = t182 * qJD(2);
t300 = t112 * t309 + t216 * t66 + t217 * t65;
t292 = -pkin(2) - t373;
t291 = t310 / 0.2e1;
t290 = t309 / 0.2e1;
t289 = -t307 / 0.2e1;
t288 = t307 / 0.2e1;
t287 = -t306 / 0.2e1;
t286 = t306 / 0.2e1;
t243 = -pkin(3) * t225 - t164;
t285 = qJD(2) * t111 + t403 * t223;
t284 = (-t163 * t216 + t358) * qJD(2) + t402 * t223;
t283 = qJD(2) * t109 + t110 * t223 - t216 * t345;
t282 = t111 * t223 - t217 * t345 + (-t216 * t260 + t354) * qJD(2);
t97 = t120 * t339;
t281 = t217 * t116 - t97;
t280 = -t106 + t350;
t279 = -t115 + t348;
t278 = t401 * t223;
t277 = t163 * t223 - t345;
t135 = t398 * t223;
t273 = -pkin(3) * t304 - t135;
t86 = t111 * t341;
t272 = t109 * t342 - t86;
t155 = rSges(3,1) * t217 - rSges(3,2) * t216;
t154 = rSges(3,1) * t216 + rSges(3,2) * t217;
t188 = -rSges(4,2) * t225 + t373;
t70 = t118 * t226 + t120 * t225;
t246 = qJD(3) * t183;
t74 = -t217 * t246 + (-t216 * t261 + t355) * qJD(2);
t247 = qJD(3) * t185;
t76 = -t217 * t247 + (-t186 * t216 + t359) * qJD(2);
t232 = -qJD(3) * t70 - t225 * t74 + t226 * t76 + t311;
t69 = t117 * t226 + t119 * t225;
t75 = qJD(2) * t118 - t216 * t246;
t77 = qJD(2) * t120 - t216 * t247;
t233 = qJD(2) * t115 - qJD(3) * t69 - t225 * t75 + t226 * t77;
t270 = -(t216 * t396 + t233 * t217) * t217 + (t216 * t397 + t232 * t217) * t216;
t269 = -(t233 * t216 - t217 * t396) * t217 + (t232 * t216 - t217 * t397) * t216;
t268 = -t216 * t42 - t217 * t41;
t48 = -t118 * t340 - t281;
t267 = t216 * t48 - t217 * t47;
t49 = -t117 * t335 - t365;
t50 = -t118 * t335 + t364;
t266 = t216 * t50 - t217 * t49;
t265 = -t68 * t216 - t367;
t80 = -rSges(4,2) * t217 * t304 + (-t226 * t310 - t293) * rSges(4,1) + t316;
t142 = t187 * t216;
t81 = -qJD(3) * t142 + (t188 * t217 + t208) * qJD(2);
t264 = t216 * t81 + t217 * t80;
t51 = t108 * t219 + t110 * t218;
t255 = t121 * t216 + t122 * t217;
t253 = t183 * t226 + t185 * t225;
t242 = qJD(2) * t159 - t152 * t346 + t153 * t347;
t241 = -t225 * t322 + t226 * t324;
t240 = (-t225 * t314 + t226 * t315) * qJD(2);
t235 = -t218 * t282 + t219 * t284 + t312;
t10 = t216 * t395 + t235 * t217;
t236 = qJD(2) * t106 - t218 * t283 + t219 * t285;
t11 = t236 * t216 - t217 * t394;
t12 = t235 * t216 - t217 * t395;
t43 = -t106 * t217 - t248;
t44 = -t107 * t217 - t272;
t20 = t152 * t44 - t153 * t43 + t331;
t45 = -t108 * t337 - t375;
t46 = -t109 * t337 + t374;
t58 = -t217 * t254 + t347;
t53 = t58 * qJD(2);
t21 = t152 * t46 - t153 * t45 + t53;
t237 = t402 * t152 - t403 * t153 + (-t160 + t163) * qJD(2);
t230 = -t218 * t390 + t237 * t219;
t28 = t218 * t285 + t219 * t283;
t29 = t218 * t284 + t219 * t282;
t234 = qJD(2) * t158 - t218 * t278 + t219 * t277;
t30 = t216 * t393 + t234 * t217;
t31 = t234 * t216 - t217 * t393;
t52 = t109 * t219 + t111 * t218;
t9 = t216 * t394 + t236 * t217;
t239 = (qJD(2) * t30 + qJDD(2) * t58 + t10 * t152 + t104 * t46 + t105 * t45 - t153 * t9) * t380 + (t237 * t218 + t219 * t390) * t378 + t20 * t291 + t21 * t290 + (qJD(2) * t31 + qJDD(2) * t57 + t104 * t44 + t105 * t43 - t11 * t153 + t12 * t152) * t379 + (t216 * t46 - t217 * t45) * t388 + (t216 * t44 - t217 * t43) * t387 + (t10 * t216 - t217 * t9 + (t216 * t45 + t217 * t46) * qJD(2)) * t383 + (t216 * t52 - t217 * t51) * t366 + (-t11 * t217 + t12 * t216 + (t216 * t43 + t217 * t44) * qJD(2)) * t382 + (t216 * t29 - t217 * t28 + (t216 * t51 + t217 * t52) * qJD(2)) * t377 + (t216 * t242 + t217 * t230) * t384 + (t216 * t230 - t217 * t242) * t381;
t167 = t261 * qJD(3);
t168 = t186 * qJD(3);
t231 = qJD(2) * t181 - qJD(3) * t253 - t167 * t225 + t168 * t226;
t170 = t188 * qJD(3);
t151 = qJD(2) * t156;
t147 = t157 * qJD(2);
t79 = -t217 * t252 + t344;
t71 = t79 * qJD(2);
t56 = qJD(3) * t255 + qJD(1);
t40 = t231 * t216 - t217 * t392;
t39 = t216 * t392 + t231 * t217;
t37 = qJD(2) * t80 + qJDD(2) * t122 - t148 * t187 - t170 * t307 + t320;
t36 = -t170 * t306 + t149 * t187 + t321 * qJDD(2) + (-t147 - t81) * qJD(2);
t34 = -qJD(3) * t256 + t225 * t76 + t226 * t74;
t33 = -t257 * qJD(3) + t225 * t77 + t226 * t75;
t32 = qJD(3) * t264 + t121 * t148 - t122 * t149 + qJDD(1);
t27 = qJD(3) * t266 + t71;
t26 = qJD(3) * t267 + t330;
t23 = -t104 * t164 - t135 * t152 + t328 * qJDD(2) + (t65 + t82) * qJD(2) + (-t148 * t225 - t216 * t332) * pkin(3) + t320;
t22 = t105 * t164 - t135 * t153 + (t149 * t225 - t217 * t332) * pkin(3) + t296 * qJDD(2) + (-t147 - t66 - t83) * qJD(2);
t1 = [t389 * qJDD(1) + m(4) * t32 + m(5) * t8 + (-m(4) - m(5) - t389) * g(3); (t71 + ((t48 - t97 + (t116 + t349) * t217 + t365) * t217 + t364 * t216) * qJD(3)) * t286 + (t53 + (t44 + (t107 + t351) * t217 + t272 + t375) * t153 + (-t217 * t280 - t400 + t43) * t152) * t381 - m(3) * (-g(1) * t154 + g(2) * t155) + (t52 + t58) * t388 + (t51 + t57) * t387 + (t70 + t79) * t386 + (t69 + t78) * t385 + (-t331 + (t46 + t400) * t153 + (t280 * t216 + t45 - t86) * t152 + ((t107 + t259) * t152 + t280 * t153) * t217 + t20) * t384 + (t29 + t30) * t383 + (-t330 + ((t217 * t279 - t364 + t50) * t217 + (t216 * t279 + t281 + t49) * t216) * qJD(3) + t26) * t289 + (t34 + t39) * t288 + (-qJD(3) * t252 + t167 * t226 + t168 * t225 + t277 * t218 + t278 * t219) * qJD(2) + (t21 + t28 + t31) * t382 + (t33 + t40 + t27) * t287 + (t41 * (t179 + (rSges(5,1) * t333 + t299) * t216) + t42 * t318 + (-t251 - t298) * t368 + ((t41 * (-t398 - t376) - t42 * t227) * t217 + (t41 * (t227 - rSges(5,3)) + t42 * (-t215 - t376)) * t216) * qJD(2) - (qJD(2) * t329 - t151 + t244 - t41) * t42 + (t23 - g(2)) * (t113 + t276) + (t22 - g(1)) * (-t112 + t317)) * m(5) + (t68 * (t199 + t316) + (t187 * t370 - t371) * qJD(3) + ((-pkin(2) - t188) * t367 + (t67 * (-rSges(4,3) - pkin(5)) + t68 * t292) * t216) * qJD(2) - (-qJD(2) * t121 - t151 - t294 - t67) * t68 + (t37 - g(2)) * t90 + (t36 - g(1)) * (t292 * t216 + t212 + t313)) * m(4) + (t253 + m(3) * (t154 ^ 2 + t155 ^ 2) + Icges(3,3) + t160 * t219 + t162 * t218) * qJDD(2); ((t49 * t216 + t50 * t217) * qJD(2) + t270) * t288 + ((t47 * t216 + t48 * t217) * qJD(2) + t269) * t287 + ((t225 * t315 + t226 * t314) * qJD(2) + ((t216 * t322 - t217 * t323) * t226 + (t216 * t324 + t217 * t325) * t225) * qJD(3)) * t378 + ((-t306 * t344 - t303) * t217 + (t240 + (t241 * t216 + (t343 - t391) * t217) * qJD(3)) * t216) * t286 + ((-t307 * t343 + t303) * t216 + (t240 + (-t391 * t217 + (t344 + t241) * t216) * qJD(3)) * t217) * t289 + t239 + (t216 * t34 - t217 * t33 + (t69 * t216 + t217 * t70) * qJD(2)) * t377 + (qJD(2) * t40 + qJD(3) * t269 + qJDD(2) * t78 + t148 * t48 + t149 * t47) * t379 + (qJD(2) * t39 + qJD(3) * t270 + qJDD(2) * t79 + t148 * t50 + t149 * t49) * t380 + t26 * t291 + t27 * t290 + (t216 * t70 - t217 * t69) * t366 + t266 * t386 + t267 * t385 + (-(-t308 * t368 + (t268 * t226 + t35 * (-t216 ^ 2 - t217 ^ 2) * t225) * qJD(3)) * pkin(3) + t35 * t300 + (t22 * t243 + t41 * t273 + t8 * t101 + t35 * t82 + (-t35 * t100 + t243 * t42) * qJD(2)) * t217 + (t23 * t243 + t42 * t273 - t8 * t100 + t35 * t83 + (t41 * t164 - t328 * t35) * qJD(2)) * t216 - g(3) * (t398 + t221) - (g(1) * t217 + g(2) * t216) * t243 + t399) * m(5) + (t32 * t255 + t56 * ((t121 * t217 - t122 * t216) * qJD(2) + t264) + t265 * t170 + (-t37 * t216 - t36 * t217 + (-t217 * t68 + t370) * qJD(2)) * t187 - (t142 * t67 - t371) * qJD(2) - (t56 * (-t142 * t216 - t143 * t217) + t265 * t188) * qJD(3) + g(1) * t143 + g(2) * t142 - g(3) * t188) * m(4); t239 + (t35 * (-t113 * t310 + t300) + t268 * t135 + (-t23 * t216 - t22 * t217 + (t41 * t216 - t368) * qJD(2)) * t164 + g(1) * t131 + g(2) * t130 - g(3) * t398 + t399) * m(5);];
tau = t1;
