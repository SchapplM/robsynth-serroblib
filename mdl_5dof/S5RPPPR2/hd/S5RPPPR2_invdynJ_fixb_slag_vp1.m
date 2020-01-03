% Calculate vector of inverse dynamics joint torques for
% S5RPPPR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:20
% EndTime: 2020-01-03 11:22:42
% DurationCPUTime: 11.00s
% Computational Cost: add. (9940->526), mult. (26586->681), div. (0->0), fcn. (30448->10), ass. (0->241)
t252 = sin(pkin(9));
t255 = cos(pkin(9));
t256 = cos(pkin(7));
t254 = sin(pkin(7));
t358 = cos(pkin(8));
t308 = t254 * t358;
t204 = t252 * t308 + t256 * t255;
t170 = qJD(5) * t204 + qJD(1);
t258 = sin(qJ(1));
t307 = t258 * t358;
t253 = sin(pkin(8));
t260 = cos(qJ(1));
t346 = t260 * t253;
t207 = t256 * t307 - t346;
t350 = t254 * t258;
t160 = t207 * t255 + t252 * t350;
t306 = t260 * t358;
t347 = t258 * t253;
t206 = t256 * t347 + t306;
t257 = sin(qJ(5));
t259 = cos(qJ(5));
t121 = -t160 * t257 + t206 * t259;
t122 = t160 * t259 + t206 * t257;
t159 = t207 * t252 - t255 * t350;
t60 = Icges(6,5) * t122 + Icges(6,6) * t121 + Icges(6,3) * t159;
t357 = Icges(6,4) * t122;
t63 = Icges(6,2) * t121 + Icges(6,6) * t159 + t357;
t115 = Icges(6,4) * t121;
t66 = Icges(6,1) * t122 + Icges(6,5) * t159 + t115;
t14 = t121 * t63 + t122 * t66 + t159 * t60;
t209 = t256 * t306 + t347;
t349 = t254 * t260;
t162 = t209 * t255 + t252 * t349;
t208 = t256 * t346 - t307;
t124 = t162 * t259 + t208 * t257;
t125 = t162 * t257 - t208 * t259;
t163 = -t209 * t252 + t255 * t349;
t61 = Icges(6,5) * t124 - Icges(6,6) * t125 - Icges(6,3) * t163;
t356 = Icges(6,4) * t124;
t64 = -Icges(6,2) * t125 - Icges(6,6) * t163 + t356;
t116 = Icges(6,4) * t125;
t67 = Icges(6,1) * t124 - Icges(6,5) * t163 - t116;
t15 = -t121 * t64 - t122 * t67 - t159 * t61;
t291 = t14 * t159 + t15 * t163;
t205 = -t256 * t252 + t255 * t308;
t351 = t253 * t254;
t157 = -t205 * t257 + t259 * t351;
t158 = t205 * t259 + t257 * t351;
t94 = Icges(6,5) * t158 + Icges(6,6) * t157 + Icges(6,3) * t204;
t355 = Icges(6,4) * t158;
t95 = Icges(6,2) * t157 + Icges(6,6) * t204 + t355;
t148 = Icges(6,4) * t157;
t96 = Icges(6,1) * t158 + Icges(6,5) * t204 + t148;
t30 = t121 * t95 + t122 * t96 + t159 * t94;
t5 = qJD(5) * t291 + t30 * t170;
t333 = qJD(5) * t163;
t311 = t333 / 0.2e1;
t16 = -t124 * t66 + t125 * t63 + t163 * t60;
t398 = -t124 * t67 + t125 * t64 + t163 * t61;
t290 = t159 * t16 - t398 * t163;
t21 = -t157 * t64 - t158 * t67 - t204 * t61;
t309 = t162 * pkin(4) - pkin(6) * t163;
t305 = t209 * pkin(3) + t208 * qJ(4);
t354 = qJ(2) * t258;
t395 = t305 + t354;
t70 = t124 * rSges(6,1) - rSges(6,2) * t125 - rSges(6,3) * t163;
t411 = t70 + t309 + t395;
t369 = pkin(2) * t256;
t292 = qJ(3) * t254 + t369;
t211 = t292 * t260;
t222 = pkin(1) * t260 + t354;
t217 = qJD(1) * t222;
t335 = qJD(3) * t254;
t234 = t258 * t335;
t336 = qJD(2) * t260;
t304 = t234 - t336;
t337 = qJD(1) * t260;
t317 = t256 * t337;
t318 = t254 * t337;
t321 = pkin(2) * t317 + qJ(3) * t318 + t234;
t338 = qJD(1) * t258;
t340 = pkin(1) * t337 + qJ(2) * t338;
t408 = -qJD(1) * t211 - t217 - t304 + t321 + t340;
t97 = rSges(6,1) * t158 + rSges(6,2) * t157 + rSges(6,3) * t204;
t407 = t170 * t70 + t333 * t97;
t283 = pkin(1) + t292;
t385 = t283 * t260;
t189 = qJD(4) * t206;
t183 = -qJD(1) * t307 + t253 * t317;
t184 = t209 * qJD(1);
t276 = t184 * pkin(3) + qJ(4) * t183 + t189;
t405 = -qJD(1) * t305 - t189 + t276 + t408;
t144 = t184 * t252 - t255 * t318;
t145 = t184 * t255 + t252 * t318;
t404 = t145 * pkin(4) + pkin(6) * t144 - t336;
t326 = t162 * rSges(5,1) + rSges(5,2) * t163 + t208 * rSges(5,3);
t402 = t395 + t326;
t399 = t124 * t96 - t125 * t95 - t163 * t94;
t396 = t159 * t97;
t324 = t209 * rSges(4,1) - t208 * rSges(4,2) + rSges(4,3) * t349;
t394 = t324 + t354;
t348 = t256 * t258;
t382 = -rSges(3,1) * t348 + rSges(3,2) * t350;
t379 = rSges(3,3) * t260 + t382;
t393 = qJD(1) * t379;
t384 = t184 * rSges(4,1) - t183 * rSges(4,2) + rSges(4,3) * t318 - t336;
t383 = t145 * rSges(5,1) - t144 * rSges(5,2) + t183 * rSges(5,3) - t336;
t136 = t207 * pkin(3) + qJ(4) * t206;
t210 = pkin(2) * t348 + qJ(3) * t350;
t235 = t260 * t335;
t251 = t258 * pkin(1);
t220 = -qJ(2) * t260 + t251;
t249 = qJD(2) * t258;
t343 = qJD(1) * t220 - t249;
t301 = qJD(1) * t210 - t235 + t343;
t332 = t208 * qJD(4);
t381 = -qJD(1) * t136 - t301 + t332;
t181 = t206 * qJD(1);
t182 = t207 * qJD(1);
t114 = pkin(3) * t182 + t181 * qJ(4) - t332;
t339 = qJ(2) * t337 + t249;
t320 = t235 + t339;
t380 = -t283 * t338 - t114 + t320;
t111 = t160 * pkin(4) + pkin(6) * t159;
t69 = t122 * rSges(6,1) + t121 * rSges(6,2) + t159 * rSges(6,3);
t378 = -qJD(1) * t111 - t170 * t69 + t381;
t98 = t160 * rSges(5,1) - t159 * rSges(5,2) + t206 * rSges(5,3);
t377 = -qJD(1) * t98 + t381;
t376 = t256 ^ 2;
t375 = -m(5) - m(6);
t104 = qJD(5) * t144 + qJDD(5) * t159;
t374 = t104 / 0.2e1;
t319 = t254 * t338;
t142 = t182 * t252 - t255 * t319;
t105 = qJD(5) * t142 + qJDD(5) * t163;
t373 = t105 / 0.2e1;
t372 = -t256 / 0.2e1;
t371 = t258 / 0.2e1;
t370 = -t260 / 0.2e1;
t368 = g(2) * t260;
t169 = qJDD(5) * t204 + qJDD(1);
t146 = t157 * qJD(5);
t147 = t158 * qJD(5);
t100 = Icges(6,5) * t146 - Icges(6,6) * t147;
t101 = Icges(6,4) * t146 - Icges(6,2) * t147;
t102 = Icges(6,1) * t146 - Icges(6,4) * t147;
t19 = t100 * t204 + t101 * t157 + t102 * t158 + t146 * t96 - t147 * t95;
t34 = t157 * t95 + t158 * t96 + t204 * t94;
t366 = t34 * t169 + t19 * t170;
t365 = rSges(3,1) * t256;
t364 = rSges(3,2) * t254;
t20 = t157 * t63 + t158 * t66 + t204 * t60;
t363 = t20 * t104;
t362 = t21 * t105;
t361 = rSges(3,3) + qJ(2);
t360 = -Icges(6,2) * t158 + t148 + t96;
t359 = Icges(6,1) * t157 - t355 - t95;
t203 = pkin(1) * t338 - t339;
t345 = -t292 * t338 - t203 + t235;
t344 = t211 + t222;
t241 = rSges(3,2) * t349;
t179 = rSges(3,3) * t258 + t260 * t365 - t241;
t342 = t222 + t179;
t341 = rSges(3,1) * t317 + rSges(3,3) * t338;
t334 = qJD(5) * t159;
t331 = -m(4) + t375;
t330 = qJDD(3) * t254;
t77 = -qJD(5) * t122 - t145 * t257 + t183 * t259;
t78 = qJD(5) * t121 + t145 * t259 + t183 * t257;
t42 = t78 * rSges(6,1) + t77 * rSges(6,2) + t144 * rSges(6,3);
t328 = -t114 + t345;
t127 = t207 * rSges(4,1) - t206 * rSges(4,2) + rSges(4,3) * t350;
t323 = t305 + t344;
t322 = t324 + t344;
t316 = pkin(1) + t365;
t315 = qJD(1) * t335;
t313 = t334 / 0.2e1;
t312 = -t333 / 0.2e1;
t303 = t326 + t323;
t302 = qJD(1) * t249 - qJDD(2) * t260;
t212 = -qJDD(3) * t256 + qJDD(4) * t351;
t300 = t309 + t323;
t143 = t182 * t255 + t252 * t319;
t299 = -pkin(4) * t143 - pkin(6) * t142;
t75 = qJD(5) * t124 - t143 * t257 + t181 * t259;
t76 = qJD(5) * t125 + t143 * t259 + t181 * t257;
t35 = Icges(6,5) * t76 + Icges(6,6) * t75 + Icges(6,3) * t142;
t36 = Icges(6,5) * t78 + Icges(6,6) * t77 + Icges(6,3) * t144;
t37 = Icges(6,4) * t76 + Icges(6,2) * t75 + Icges(6,6) * t142;
t38 = Icges(6,4) * t78 + Icges(6,2) * t77 + Icges(6,6) * t144;
t39 = Icges(6,1) * t76 + Icges(6,4) * t75 + Icges(6,5) * t142;
t40 = Icges(6,1) * t78 + Icges(6,4) * t77 + Icges(6,5) * t144;
t298 = (-t124 * t40 + t125 * t38 + t142 * t60 + t163 * t36 + t63 * t75 + t66 * t76) * t159 + t163 * (-t124 * t39 + t125 * t37 - t142 * t61 + t163 * t35 - t64 * t75 - t67 * t76);
t297 = t159 * (t121 * t38 + t122 * t40 + t144 * t60 + t159 * t36 + t63 * t77 + t66 * t78) + t163 * (t121 * t37 + t122 * t39 - t144 * t61 + t159 * t35 - t64 * t77 - t67 * t78);
t7 = t146 * t66 - t147 * t63 + t157 * t38 + t158 * t40 + t204 * t36;
t8 = -t146 * t67 + t147 * t64 + t157 * t37 + t158 * t39 + t204 * t35;
t296 = t159 * t7 + t163 * t8;
t295 = t189 + t304;
t294 = -t336 + t340;
t223 = rSges(2,1) * t260 - t258 * rSges(2,2);
t221 = rSges(2,1) * t258 + rSges(2,2) * t260;
t293 = -rSges(4,1) * t182 + rSges(4,2) * t181;
t41 = rSges(6,1) * t76 + rSges(6,2) * t75 + rSges(6,3) * t142;
t289 = t159 * t41 - t163 * t42;
t288 = -t159 * t70 - t163 * t69;
t287 = t159 * (Icges(6,5) * t121 - Icges(6,6) * t122) + t163 * (Icges(6,5) * t125 + Icges(6,6) * t124);
t286 = qJD(1) * t294 + qJDD(1) * t220 - qJDD(2) * t258;
t282 = t220 + t210;
t278 = t258 * t330 + t260 * t315 + t302;
t275 = -rSges(5,1) * t143 + rSges(5,2) * t142 - rSges(5,3) * t181;
t272 = (Icges(6,1) * t121 - t357 - t63) * t159 + (Icges(6,1) * t125 + t356 + t64) * t163;
t271 = (-Icges(6,2) * t122 + t115 + t66) * t159 + (Icges(6,2) * t124 + t116 - t67) * t163;
t269 = t283 * t368;
t268 = qJD(4) * t183 + qJDD(4) * t206 + t278;
t266 = t136 + t282;
t263 = qJD(1) * t321 + qJDD(1) * t210 + t258 * t315 - t260 * t330 + t286;
t262 = qJD(1) * t276 + qJD(4) * t181 + qJDD(1) * t136 - qJDD(4) * t208 + t263;
t133 = qJD(1) * t342 - t336;
t132 = t343 - t393;
t110 = rSges(6,1) * t157 - rSges(6,2) * t158;
t107 = Icges(6,5) * t157 - Icges(6,6) * t158;
t103 = rSges(6,1) * t146 - rSges(6,2) * t147;
t93 = qJD(1) * t322 + t304;
t91 = t342 * qJDD(1) + (-t203 + t393) * qJD(1) + t302;
t90 = -qJDD(1) * t379 + ((-qJD(1) * t364 - qJD(2)) * t260 + t341) * qJD(1) + t286;
t86 = rSges(6,1) * t125 + rSges(6,2) * t124;
t85 = rSges(6,1) * t121 - rSges(6,2) * t122;
t56 = qJD(1) * t303 + t295;
t45 = t322 * qJDD(1) + (-rSges(4,3) * t319 + t293 + t345) * qJD(1) + t278;
t44 = qJD(1) * t384 + qJDD(1) * t127 + t263;
t32 = -qJD(3) * t256 + qJD(4) * t351 + qJD(5) * t288;
t28 = t303 * qJDD(1) + (t275 + t328) * qJD(1) + t268;
t27 = qJD(1) * t383 + qJDD(1) * t98 + t262;
t26 = qJD(1) * t300 + t295 + t407;
t25 = -t334 * t97 - t378;
t13 = t100 * t159 + t101 * t121 + t102 * t122 + t144 * t94 + t77 * t95 + t78 * t96;
t12 = t100 * t163 + t101 * t125 - t102 * t124 + t142 * t94 + t75 * t95 + t76 * t96;
t11 = qJD(5) * t289 - t104 * t70 - t105 * t69 + t212;
t10 = t103 * t333 + t105 * t97 + t169 * t70 - t170 * t41 + t300 * qJDD(1) + (t299 + t328) * qJD(1) + t268;
t9 = t404 * qJD(1) + qJDD(1) * t111 - t103 * t334 - t104 * t97 + t169 * t69 + t170 * t42 + t262;
t1 = [t366 + t363 / 0.2e1 + t362 / 0.2e1 - t399 * t373 + t30 * t374 + t5 * t312 - m(2) * (g(2) * t223 + g(3) * t221) + (t13 + t7) * t313 + (-g(2) * t411 - t269 + (t9 - g(3)) * (t111 + t266 + t69) + (-qJD(5) * t396 + t299 - t378 + t380 - t41) * t26 + (t385 + t411) * t10 + (-qJD(1) * t309 + t404 + t405 - t407 + t42) * t25) * m(6) + (-g(2) * t402 - t269 + (-g(3) + t27) * (t266 + t98) + (t275 + t380 - t377) * t56 + (t385 + t402) * t28 - (-qJD(1) * t326 + t383 + t405) * t377) * m(5) + (-g(2) * t394 - t269 + (t293 + t320 + (-t369 - pkin(1) + (-rSges(4,3) - qJ(3)) * t254) * t338) * t93 + (t385 + t394) * t45 + (-qJD(1) * t324 + t384 + t408 + t93) * (qJD(1) * t127 + t301) + (-g(3) + t44) * (t282 + t127)) * m(4) + (t133 * t339 + t132 * (t294 + t341) + ((rSges(3,3) * t133 - t132 * t364) * t260 + t133 * (-t316 + t364) * t258) * qJD(1) - (qJD(1) * t179 - t133 + t217 - t336) * t132 + (t91 - g(2)) * (t258 * t361 + t260 * t316 - t241) + (t90 - g(3)) * (-t260 * t361 + t251 - t382)) * m(3) + (Icges(3,2) * t376 + (Icges(5,1) * t205 + 0.2e1 * Icges(5,5) * t351) * t205 + (-0.2e1 * Icges(5,4) * t205 + Icges(5,2) * t204 - 0.2e1 * Icges(5,6) * t351) * t204 + m(2) * (t221 ^ 2 + t223 ^ 2) + Icges(2,3) + (-Icges(4,5) * t308 + Icges(4,6) * t351 + Icges(4,3) * t256) * t256 + ((Icges(4,1) * t358 - Icges(4,4) * t253) * t308 - (Icges(4,4) * t358 - Icges(4,2) * t253) * t351 + (Icges(5,3) * t253 ^ 2 + Icges(3,1)) * t254 + (-Icges(4,5) * t358 + Icges(4,6) * t253 + (2 * Icges(3,4))) * t256) * t254) * qJDD(1) + (t5 + t12 + t8) * t311; (-m(3) + t331) * (-g(3) * t258 - t368) + m(3) * (-t258 * t90 - t260 * t91) + m(4) * (-t258 * t44 - t260 * t45) + m(5) * (-t27 * t258 - t260 * t28) + m(6) * (-t10 * t260 - t258 * t9); t331 * (-g(1) * t256 + (g(2) * t258 - g(3) * t260) * t254) + 0.2e1 * (t11 * t372 + (t10 * t371 + t370 * t9) * t254) * m(6) + 0.2e1 * (t212 * t372 + (t27 * t370 + t28 * t371) * t254) * m(5) + 0.2e1 * (qJDD(3) * t376 / 0.2e1 + (t370 * t44 + t371 * t45) * t254) * m(4); t375 * (g(1) * t351 + g(2) * t206 - g(3) * t208) + m(5) * (-t181 * t377 + t183 * t56 + t206 * t28 - t208 * t27 + t212 * t351) + m(6) * (t10 * t206 + t11 * t351 + t181 * t25 + t183 * t26 - t208 * t9) + 0.2e1 * (-m(5) * (-t206 * t377 + t208 * t56) / 0.2e1 - m(6) * (t206 * t25 + t208 * t26) / 0.2e1) * qJD(1); t204 * (qJD(5) * t296 + t362 + t363 + t366) / 0.2e1 + t169 * (t159 * t20 + t163 * t21 + t204 * t34) / 0.2e1 + t170 * (t142 * t21 + t144 * t20 + t19 * t204 + t296) / 0.2e1 + t144 * t5 / 0.2e1 + t159 * (qJD(5) * t297 + t104 * t14 + t105 * t15 + t13 * t170 + t169 * t30) / 0.2e1 + (t204 * t30 + t291) * t374 + (t13 * t204 + t14 * t144 + t142 * t15 + t297) * t313 + t142 * (qJD(5) * t290 - t399 * t170) / 0.2e1 + t163 * (qJD(5) * t298 + t104 * t16 - t105 * t398 + t12 * t170 - t169 * t399) / 0.2e1 + (-t204 * t399 + t290) * t373 + (t12 * t204 - t142 * t398 + t144 * t16 + t298) * t311 - t170 * ((t204 * t107 + t157 * t360 + t158 * t359) * t170 + (t157 * t271 + t158 * t272 + t204 * t287) * qJD(5)) / 0.2e1 - ((t159 * t107 + t121 * t360 + t122 * t359) * t170 + (t121 * t271 + t122 * t272 + t159 * t287) * qJD(5)) * t334 / 0.2e1 + ((t163 * t107 - t124 * t359 + t125 * t360) * t170 + (-t124 * t272 + t125 * t271 + t163 * t287) * qJD(5)) * t312 + (t11 * t288 + t32 * (-t142 * t69 - t144 * t70 + t289) + t10 * (t163 * t97 + t204 * t70) + t26 * (t103 * t163 + t142 * t97 - t204 * t41) + t9 * (t204 * t69 - t396) + t25 * (-t103 * t159 - t144 * t97 + t204 * t42) - (t25 * t85 - t26 * t86) * t170 - (t32 * (t159 * t86 - t163 * t85) + (-t159 * t25 + t163 * t26) * t110) * qJD(5) - g(1) * t110 - g(2) * t85 - g(3) * t86) * m(6);];
tau = t1;
