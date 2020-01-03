% Calculate vector of inverse dynamics joint torques for
% S4PRPR5
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
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR5_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR5_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR5_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR5_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:58
% EndTime: 2019-12-31 16:23:16
% DurationCPUTime: 15.22s
% Computational Cost: add. (7388->514), mult. (11340->779), div. (0->0), fcn. (10736->8), ass. (0->269)
t211 = qJ(2) + pkin(7);
t206 = sin(t211);
t207 = cos(t211);
t216 = sin(qJ(2));
t218 = cos(qJ(2));
t278 = -Icges(3,5) * t216 - Icges(3,6) * t218;
t422 = -Icges(4,5) * t206 - Icges(4,6) * t207 + t278;
t426 = Icges(3,3) + Icges(4,3);
t425 = t422 * qJD(2);
t424 = Icges(3,5) * t218 + Icges(4,5) * t207 - Icges(3,6) * t216 - Icges(4,6) * t206;
t212 = sin(pkin(6));
t213 = cos(pkin(6));
t423 = t212 * t213;
t210 = t213 ^ 2;
t193 = rSges(3,1) * t216 + rSges(3,2) * t218;
t209 = t212 ^ 2;
t392 = t209 + t210;
t421 = t193 * t392;
t420 = t424 * t212 - t426 * t213;
t419 = t426 * t212 + t424 * t213;
t418 = t425 * t212;
t417 = t425 * t213;
t215 = sin(qJ(4));
t368 = rSges(5,2) * t215;
t217 = cos(qJ(4));
t370 = rSges(5,1) * t217;
t301 = -t368 + t370;
t416 = t206 * rSges(5,3) + t207 * t301;
t361 = Icges(3,4) * t218;
t284 = -Icges(3,2) * t216 + t361;
t132 = -Icges(3,6) * t213 + t212 * t284;
t362 = Icges(3,4) * t216;
t289 = Icges(3,1) * t218 - t362;
t134 = -Icges(3,5) * t213 + t212 * t289;
t283 = -Icges(3,2) * t218 - t362;
t288 = -Icges(3,1) * t216 - t361;
t334 = qJD(2) * t207;
t335 = qJD(2) * t206;
t359 = Icges(4,4) * t207;
t360 = Icges(4,4) * t206;
t282 = -Icges(4,2) * t206 + t359;
t114 = -Icges(4,6) * t213 + t212 * t282;
t287 = Icges(4,1) * t207 - t360;
t116 = -Icges(4,5) * t213 + t212 * t287;
t403 = t114 * t207 + t116 * t206;
t414 = ((-Icges(4,2) * t207 - t360) * t335 - (-Icges(4,1) * t206 - t359) * t334 - (-t216 * t283 + t218 * t288) * qJD(2)) * t212 + (t132 * t218 + t134 * t216 + t403) * qJD(2);
t412 = t114 * t206 - t116 * t207 + t132 * t216 - t134 * t218;
t411 = 0.2e1 * qJD(2);
t410 = 2 * qJDD(2);
t376 = t216 * pkin(2);
t408 = t392 * t376;
t407 = t422 * t212;
t406 = t422 * t213;
t208 = t218 * pkin(2);
t394 = t207 * rSges(4,1) - rSges(4,2) * t206;
t396 = t394 + t208;
t393 = t207 * pkin(3) + t206 * pkin(5);
t391 = g(1) * t213 + g(2) * t212;
t275 = Icges(5,5) * t217 - Icges(5,6) * t215;
t102 = -Icges(5,3) * t207 + t206 * t275;
t355 = Icges(5,4) * t217;
t280 = -Icges(5,2) * t215 + t355;
t104 = -Icges(5,6) * t207 + t206 * t280;
t356 = Icges(5,4) * t215;
t285 = Icges(5,1) * t217 - t356;
t106 = -Icges(5,5) * t207 + t206 * t285;
t329 = qJD(4) * t206;
t333 = qJD(2) * t212;
t177 = t213 * t329 + t333;
t331 = qJD(2) * t213;
t178 = t212 * t329 - t331;
t344 = t213 * t215;
t345 = t212 * t217;
t164 = -t207 * t344 + t345;
t349 = t206 * t213;
t343 = t213 * t217;
t346 = t212 * t215;
t165 = t207 * t343 + t346;
t357 = Icges(5,4) * t165;
t67 = Icges(5,2) * t164 + Icges(5,6) * t349 + t357;
t145 = Icges(5,4) * t164;
t69 = Icges(5,1) * t165 + Icges(5,5) * t349 + t145;
t293 = -t215 * t67 + t217 * t69;
t162 = -t207 * t346 - t343;
t350 = t206 * t212;
t163 = t207 * t345 - t344;
t358 = Icges(5,4) * t163;
t66 = Icges(5,2) * t162 + Icges(5,6) * t350 + t358;
t144 = Icges(5,4) * t162;
t68 = Icges(5,1) * t163 + Icges(5,5) * t350 + t144;
t294 = -t215 * t66 + t217 * t68;
t390 = -(-t102 * t213 - t293) * t177 - (-t102 * t212 - t294) * t178;
t389 = -t216 * (t283 * t212 + t134) - t218 * (-t288 * t212 + t132);
t149 = (-Icges(5,2) * t217 - t356) * t206;
t328 = qJD(4) * t207;
t222 = t177 * (-Icges(5,2) * t165 + t145 + t69) + t178 * (-Icges(5,2) * t163 + t144 + t68) - t328 * (t106 + t149);
t219 = qJD(2) ^ 2;
t388 = -m(4) - m(5);
t326 = qJD(2) * qJD(4);
t243 = qJDD(4) * t206 + t207 * t326;
t325 = qJDD(2) * t212;
t100 = t213 * t243 + t325;
t387 = t100 / 0.2e1;
t324 = qJDD(2) * t213;
t101 = t212 * t243 - t324;
t386 = t101 / 0.2e1;
t170 = -qJDD(4) * t207 + t206 * t326;
t385 = t170 / 0.2e1;
t384 = -t177 / 0.2e1;
t383 = t177 / 0.2e1;
t382 = -t178 / 0.2e1;
t381 = t178 / 0.2e1;
t110 = -qJ(3) * t213 + t208 * t212;
t111 = qJ(3) * t212 + t208 * t213;
t371 = t212 * t110 + t213 * t111;
t367 = pkin(2) * qJD(2);
t65 = Icges(5,5) * t165 + Icges(5,6) * t164 + Icges(5,3) * t349;
t21 = t162 * t67 + t163 * t69 + t350 * t65;
t366 = t21 * t213;
t64 = Icges(5,5) * t163 + Icges(5,6) * t162 + Icges(5,3) * t350;
t22 = t164 * t66 + t165 * t68 + t349 * t64;
t365 = t212 * t22;
t36 = t102 * t350 + t104 * t162 + t106 * t163;
t364 = t36 * t206;
t37 = t102 * t349 + t104 * t164 + t106 * t165;
t363 = t37 * t206;
t351 = t102 * t207;
t348 = t207 * t212;
t347 = t207 * t213;
t320 = t216 * t367;
t330 = qJD(3) * t213;
t181 = -t212 * t320 - t330;
t205 = qJD(3) * t212;
t182 = -t213 * t320 + t205;
t339 = t212 * t181 + t213 * t182;
t321 = t206 * t368;
t338 = rSges(5,3) * t348 + t212 * t321;
t337 = rSges(5,3) * t347 + t213 * t321;
t323 = qJDD(3) * t213;
t322 = t206 * t370;
t319 = t218 * t367;
t318 = t110 * t333 + t111 * t331 + qJD(1);
t317 = t215 * t335;
t316 = t217 * t335;
t315 = t207 * t333;
t314 = t207 * t331;
t310 = t331 / 0.2e1;
t309 = -t328 / 0.2e1;
t308 = t328 / 0.2e1;
t183 = rSges(4,1) * t206 + rSges(4,2) * t207;
t236 = -t183 - t376;
t185 = pkin(3) * t206 - pkin(5) * t207;
t306 = -t185 - t376;
t108 = -rSges(5,3) * t207 + t206 * t301;
t305 = -t108 + t306;
t168 = t394 * qJD(2);
t304 = -t168 - t319;
t194 = rSges(3,1) * t218 - rSges(3,2) * t216;
t169 = t393 * qJD(2);
t248 = (-qJDD(2) * t216 - t218 * t219) * pkin(2);
t224 = -qJD(2) * t169 - qJDD(2) * t185 + t248;
t90 = -qJD(4) * t165 + t213 * t317;
t91 = qJD(4) * t164 - t213 * t316;
t52 = rSges(5,1) * t91 + rSges(5,2) * t90 + rSges(5,3) * t314;
t153 = (-rSges(5,1) * t215 - rSges(5,2) * t217) * t206;
t63 = qJD(2) * t416 + qJD(4) * t153;
t71 = rSges(5,1) * t165 + rSges(5,2) * t164 + rSges(5,3) * t349;
t16 = -t100 * t108 + t170 * t71 - t177 * t63 + t212 * t224 - t328 * t52 - t323;
t204 = qJDD(3) * t212;
t88 = -qJD(4) * t163 + t212 * t317;
t89 = qJD(4) * t162 - t212 * t316;
t51 = rSges(5,1) * t89 + rSges(5,2) * t88 + rSges(5,3) * t315;
t70 = rSges(5,1) * t163 + rSges(5,2) * t162 + rSges(5,3) * t350;
t17 = t101 * t108 - t170 * t70 + t178 * t63 + t213 * t224 + t328 * t51 + t204;
t300 = -t16 * t213 + t17 * t212;
t299 = t65 * t177 + t64 * t178;
t20 = t162 * t66 + t163 * t68 + t350 * t64;
t298 = t20 * t212 + t366;
t23 = t164 * t67 + t165 * t69 + t349 * t65;
t297 = t213 * t23 + t365;
t24 = t206 * t294 - t207 * t64;
t25 = t206 * t293 - t207 * t65;
t296 = t24 * t212 + t25 * t213;
t295 = -t212 * t71 + t213 * t70;
t292 = qJD(2) * t236;
t291 = qJD(2) * t306;
t290 = t110 * t325 + t111 * t324 + t181 * t333 + t182 * t331 + qJDD(1);
t274 = -t104 * t215 + t106 * t217;
t118 = -rSges(4,3) * t213 + t212 * t394;
t119 = rSges(4,3) * t212 + t213 * t394;
t271 = t118 * t212 + t119 * t213;
t268 = t392 * t194;
t146 = t393 * t212;
t147 = t393 * t213;
t267 = t146 * t212 + t147 * t213;
t266 = qJD(2) * t421;
t265 = -t169 - t63 - t319;
t45 = Icges(5,5) * t89 + Icges(5,6) * t88 + Icges(5,3) * t315;
t261 = t206 * t45 + t334 * t64;
t46 = Icges(5,5) * t91 + Icges(5,6) * t90 + Icges(5,3) * t314;
t260 = t206 * t46 + t334 * t65;
t103 = Icges(5,3) * t206 + t207 * t275;
t148 = (-Icges(5,5) * t215 - Icges(5,6) * t217) * t206;
t60 = qJD(2) * t103 + qJD(4) * t148;
t259 = t102 * t334 + t206 * t60;
t258 = t103 - t274;
t257 = qJD(2) * t185;
t18 = qJD(2) * t267 + t177 * t70 - t178 * t71 + t318;
t256 = t18 * t295;
t38 = qJD(2) * t271 + t318;
t255 = t38 * t183;
t253 = qJD(2) * t183;
t150 = (-Icges(5,1) * t215 - t355) * t206;
t235 = t278 * t212 + ((-t284 + t288) * t218 + (-t283 - t289) * t216) * t213;
t234 = t148 * t328 - t177 * (Icges(5,5) * t164 - Icges(5,6) * t165) - t178 * (Icges(5,5) * t162 - Icges(5,6) * t163);
t107 = Icges(5,5) * t206 + t207 * t285;
t105 = Icges(5,6) * t206 + t207 * t280;
t230 = t206 * t234;
t225 = -qJD(2) * t168 - qJDD(2) * t183 + t248;
t223 = (Icges(5,1) * t164 - t357 - t67) * t177 + (Icges(5,1) * t162 - t358 - t66) * t178 - (-t104 + t150) * t328;
t221 = t212 * (-(Icges(4,6) * t212 + t213 * t282) * t207 - (Icges(4,5) * t212 + t213 * t287) * t206) + t213 * t403;
t220 = (-t258 * t328 - t390) * t206;
t196 = pkin(5) * t347;
t195 = pkin(5) * t348;
t180 = t193 * t213;
t179 = t193 * t212;
t129 = t213 * t257;
t128 = t212 * t257;
t127 = t213 * t253;
t126 = t212 * t253;
t93 = t213 * t292 + t205;
t92 = t212 * t292 - t330;
t87 = -t213 * t322 + t337;
t86 = -t212 * t322 + t338;
t85 = t106 * t213;
t84 = t106 * t212;
t83 = t104 * t213;
t82 = t104 * t212;
t79 = rSges(5,1) * t164 - rSges(5,2) * t165;
t78 = rSges(5,1) * t162 - rSges(5,2) * t163;
t62 = qJD(2) * t107 + qJD(4) * t150;
t61 = qJD(2) * t105 + qJD(4) * t149;
t58 = t213 * t225 + t204;
t57 = t212 * t225 - t323;
t50 = Icges(5,1) * t91 + Icges(5,4) * t90 + Icges(5,5) * t314;
t49 = Icges(5,1) * t89 + Icges(5,4) * t88 + Icges(5,5) * t315;
t48 = Icges(5,4) * t91 + Icges(5,2) * t90 + Icges(5,6) * t314;
t47 = Icges(5,4) * t89 + Icges(5,2) * t88 + Icges(5,6) * t315;
t40 = -qJD(2) * t266 + qJDD(2) * t268 + qJDD(1);
t39 = t206 * t274 - t351;
t31 = t108 * t178 + t213 * t291 + t328 * t70 + t205;
t30 = -t108 * t177 + t212 * t291 - t328 * t71 - t330;
t19 = t271 * qJDD(2) + (-t126 * t212 - t127 * t213) * qJD(2) + t290;
t15 = (qJD(2) * t274 - t60) * t207 + (qJD(2) * t102 - t215 * t61 + t217 * t62 + (-t104 * t217 - t106 * t215) * qJD(4)) * t206;
t14 = t104 * t90 + t106 * t91 + t164 * t61 + t165 * t62 + t213 * t259;
t13 = t104 * t88 + t106 * t89 + t162 * t61 + t163 * t62 + t212 * t259;
t12 = t177 * t25 + t178 * t24 - t328 * t39;
t11 = t164 * t48 + t165 * t50 + t213 * t260 + t67 * t90 + t69 * t91;
t10 = t164 * t47 + t165 * t49 + t213 * t261 + t66 * t90 + t68 * t91;
t9 = t162 * t48 + t163 * t50 + t212 * t260 + t67 * t88 + t69 * t89;
t8 = t162 * t47 + t163 * t49 + t212 * t261 + t66 * t88 + t68 * t89;
t7 = t100 * t70 - t101 * t71 + t177 * t51 - t178 * t52 + t267 * qJDD(2) + (-t128 * t212 - t129 * t213) * qJD(2) + t290;
t6 = t177 * t23 + t178 * t22 - t328 * t37;
t5 = t177 * t21 + t178 * t20 - t328 * t36;
t4 = (qJD(2) * t293 - t46) * t207 + (qJD(2) * t65 - t215 * t48 + t217 * t50 + (-t215 * t69 - t217 * t67) * qJD(4)) * t206;
t3 = (qJD(2) * t294 - t45) * t207 + (qJD(2) * t64 - t215 * t47 + t217 * t49 + (-t215 * t68 - t217 * t66) * qJD(4)) * t206;
t2 = t10 * t178 + t100 * t23 + t101 * t22 + t11 * t177 - t14 * t328 + t170 * t37;
t1 = t100 * t21 + t101 * t20 - t13 * t328 + t170 * t36 + t177 * t9 + t178 * t8;
t26 = [m(2) * qJDD(1) + m(3) * t40 + m(4) * t19 + m(5) * t7 + (-m(2) - m(3) + t388) * g(3); -t12 * t329 / 0.2e1 + (-t10 * t213 + t11 * t212) * t383 + ((-t164 * t83 - t165 * t85) * t177 + (-t164 * t82 - t165 * t84) * t178 + (t363 + (-t105 * t164 - t107 * t165 + t365) * t207) * qJD(4) + (((t23 - t351) * qJD(4) + t299) * t207 + t220) * t213) * t384 + (t212 * t25 - t213 * t24) * t385 + (-t20 * t213 + t21 * t212) * t386 + (t212 * t23 - t213 * t22) * t387 + (t212 * t9 - t213 * t8) * t381 + ((-t162 * t83 - t163 * t85) * t177 + (-t162 * t82 - t163 * t84) * t178 + (t364 + (-t105 * t162 - t107 * t163 + t366) * t207) * qJD(4) + (((t20 - t351) * qJD(4) + t299) * t207 + t220) * t212) * t382 + (((t215 * t83 - t217 * t85 + t65) * t177 + (t215 * t82 - t217 * t84 + t64) * t178 + t39 * qJD(4)) * t206 + ((t258 * t207 + (t105 * t215 - t107 * t217 - t102) * t206 + t296) * qJD(4) + t390) * t207) * t308 - (t406 * qJD(2) * t209 + (-t389 * t213 + t221 + (t235 - t407) * t212) * t331) * t333 / 0.2e1 + ((t235 * t212 + t221 + (-t389 - t406) * t213) * t333 + t407 * t210 * qJD(2)) * t310 + ((-t3 + t6) * t213 + (t4 + t5) * t212) * t309 + (t7 * t371 + (t17 * t305 + t31 * t265 + t7 * (t147 + t71)) * t213 + (t16 * t305 + t30 * t265 + t7 * (t146 + t70)) * t212 - g(1) * (-t213 * t376 + t196 + t337) - g(2) * (-t212 * t376 + t195 + t338) - t391 * t206 * (-pkin(3) - t370) - ((t30 * t71 - t31 * t70) * t206 + (t31 * (t108 * t212 + t86) + t30 * (-t108 * t213 - t87) + t256) * t207) * qJD(4) + (t177 * t30 - t178 * t31 - g(3)) * t416 + (g(3) - (t30 * t212 + t31 * t213) * qJD(2)) * (-t393 - t208) + (t339 + (-t129 + t52) * t213 + (-t128 + t51) * t212 - t177 * t86 + t178 * t87 - (-t408 + t213 * (-pkin(3) * t349 + t196) + t212 * (-pkin(3) * t350 + t195)) * qJD(2)) * t18) * m(5) + (t19 * t371 + t38 * t339 + (t19 * t119 - t38 * t127 + t236 * t58 + t304 * t93) * t213 + (t19 * t118 - t38 * t126 + t236 * t57 + t304 * t92) * t212 - (-t38 * t408 + (-t213 * t255 - t396 * t93) * t213 + (-t212 * t255 - t396 * t92) * t212) * qJD(2) - g(3) * t396 - t391 * t236) * m(4) + (g(1) * t180 + g(2) * t179 - g(3) * t194 + t40 * t268 + (-t266 - (-t179 * t212 - t180 * t213) * qJD(2)) * (qJD(2) * t268 + qJD(1)) + (t193 * qJDD(2) + t194 * t219) * t421) * m(3) + (t2 + (t414 * t210 + (t417 * t212 - t213 * t418) * t212) * t411 + (t412 * t210 + (t419 * t212 - t213 * t420) * t212) * t410) * t212 / 0.2e1 - (t1 + (t418 * t210 + (t414 - t417) * t423) * t411 + (t420 * t210 + (t412 - t419) * t423) * t410) * t213 / 0.2e1; t388 * (g(1) * t212 - g(2) * t213) + m(4) * (t212 * t58 - t213 * t57) + m(5) * t300; t207 * t6 * t310 + t2 * t349 / 0.2e1 + (t206 * t297 - t207 * t37) * t387 + (-t14 * t207 + (t10 * t212 + t11 * t213) * t206 + (t207 * t297 + t363) * qJD(2)) * t383 + t5 * t315 / 0.2e1 + t1 * t350 / 0.2e1 + (t206 * t298 - t207 * t36) * t386 + (-t13 * t207 + (t212 * t8 + t213 * t9) * t206 + (t207 * t298 + t364) * qJD(2)) * t381 + t12 * t335 / 0.2e1 - t207 * (t100 * t25 + t101 * t24 - t15 * t328 + t170 * t39 + t177 * t4 + t178 * t3) / 0.2e1 + (t206 * t296 - t207 * t39) * t385 + (-t15 * t207 + (t212 * t3 + t213 * t4) * t206 + (t39 * t206 + t207 * t296) * qJD(2)) * t309 + (t164 * t222 + t165 * t223 - t213 * t230) * t384 + (t162 * t222 + t163 * t223 - t212 * t230) * t382 + (t234 * t207 + (-t215 * t222 + t223 * t217) * t206) * t308 + ((-t16 * t71 + t17 * t70 - t30 * t52 + t31 * t51 + (t256 + (t212 * t31 - t213 * t30) * t108) * qJD(2)) * t207 + (t31 * (-qJD(2) * t70 + t212 * t63) + t30 * (qJD(2) * t71 - t213 * t63) + t7 * t295 + t18 * (-t212 * t52 + t213 * t51) + t300 * t108) * t206 - t31 * (t153 * t178 + t328 * t78) - t30 * (-t153 * t177 - t328 * t79) - t18 * (t177 * t78 - t178 * t79) - g(1) * t79 - g(2) * t78 - g(3) * t153) * m(5);];
tau = t26;
