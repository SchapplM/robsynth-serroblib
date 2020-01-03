% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPPR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:03
% EndTime: 2019-12-31 17:38:23
% DurationCPUTime: 18.29s
% Computational Cost: add. (7740->600), mult. (22222->887), div. (0->0), fcn. (23705->8), ass. (0->276)
t446 = Icges(4,4) + Icges(3,5);
t445 = Icges(3,6) - Icges(4,6);
t267 = sin(qJ(2));
t447 = (-Icges(3,4) + Icges(4,5)) * t267;
t269 = cos(qJ(2));
t438 = -t447 + (Icges(3,2) + Icges(4,3)) * t269;
t437 = -t446 * t267 - t445 * t269;
t444 = t447 + (Icges(3,1) + Icges(4,1)) * t269;
t264 = sin(pkin(7));
t386 = t264 * t269;
t255 = Icges(4,5) * t386;
t265 = cos(pkin(7));
t395 = Icges(3,4) * t269;
t321 = -Icges(3,2) * t267 + t395;
t387 = t264 * t267;
t443 = -Icges(4,3) * t387 + t264 * t321 - t445 * t265 - t255;
t384 = t265 * t269;
t256 = Icges(4,5) * t384;
t385 = t265 * t267;
t442 = -Icges(4,3) * t385 + t445 * t264 + t265 * t321 - t256;
t441 = t444 * t264 - t446 * t265;
t440 = t446 * t264 + t444 * t265;
t262 = t264 ^ 2;
t263 = t265 ^ 2;
t372 = t262 + t263;
t439 = t437 * qJD(2);
t324 = -Icges(3,1) * t267 - t395;
t368 = qJD(2) * t267;
t436 = (Icges(4,1) * t267 - Icges(4,5) * t269 - t324) * qJD(2) * t269 - t438 * t368;
t397 = sin(pkin(8));
t398 = cos(pkin(8));
t231 = t267 * t398 - t269 * t397;
t185 = t231 * t265;
t371 = qJD(2) * t264;
t154 = -qJD(5) * t185 + t371;
t410 = -t154 / 0.2e1;
t183 = t231 * t264;
t370 = qJD(2) * t265;
t155 = -qJD(5) * t183 - t370;
t408 = -t155 / 0.2e1;
t401 = -qJD(5) / 0.2e1;
t261 = qJD(3) * t267;
t250 = t264 * t261;
t236 = pkin(2) * t267 - qJ(3) * t269;
t298 = qJD(2) * t236;
t156 = -t264 * t298 + t250;
t252 = t265 * t261;
t157 = -t265 * t298 + t252;
t435 = t264 * t156 + t265 * t157 + t298 * t372 - t261;
t230 = t267 * t397 + t269 * t398;
t268 = cos(qJ(5));
t266 = sin(qJ(5));
t390 = Icges(6,4) * t266;
t322 = Icges(6,1) * t268 - t390;
t100 = Icges(6,5) * t230 + t231 * t322;
t182 = t230 * t264;
t148 = -t182 * t266 + t265 * t268;
t213 = t230 * qJD(2);
t158 = t264 * t213;
t160 = t265 * t213;
t312 = -t182 * t268 - t265 * t266;
t60 = -Icges(6,5) * t312 + Icges(6,6) * t148 - Icges(6,3) * t183;
t392 = Icges(6,4) * t312;
t62 = Icges(6,2) * t148 - Icges(6,6) * t183 - t392;
t140 = Icges(6,4) * t148;
t64 = -Icges(6,1) * t312 - Icges(6,5) * t183 + t140;
t19 = t148 * t62 - t183 * t60 - t312 * t64;
t184 = t230 * t265;
t150 = -t184 * t266 - t264 * t268;
t311 = -t184 * t268 + t264 * t266;
t61 = -Icges(6,5) * t311 + Icges(6,6) * t150 - Icges(6,3) * t185;
t391 = Icges(6,4) * t311;
t63 = Icges(6,2) * t150 - Icges(6,6) * t185 - t391;
t141 = Icges(6,4) * t150;
t65 = -Icges(6,1) * t311 - Icges(6,5) * t185 + t141;
t20 = t148 * t63 - t183 * t61 - t312 * t65;
t315 = Icges(6,5) * t268 - Icges(6,6) * t266;
t96 = Icges(6,3) * t230 + t231 * t315;
t389 = Icges(6,4) * t268;
t318 = -Icges(6,2) * t266 + t389;
t98 = Icges(6,6) * t230 + t231 * t318;
t33 = -t100 * t312 + t148 * t98 - t183 * t96;
t424 = qJD(2) * t231;
t364 = qJD(5) * t231;
t302 = t213 * t268 - t266 * t364;
t303 = -t213 * t266 - t268 * t364;
t55 = Icges(6,5) * t302 + Icges(6,6) * t303 - Icges(6,3) * t424;
t56 = Icges(6,4) * t302 + Icges(6,2) * t303 - Icges(6,6) * t424;
t57 = Icges(6,1) * t302 + Icges(6,4) * t303 - Icges(6,5) * t424;
t159 = t264 * t424;
t90 = qJD(5) * t312 + t159 * t266;
t91 = qJD(5) * t148 - t159 * t268;
t279 = (t100 * t91 + t148 * t56 - t158 * t96 - t183 * t55 - t312 * t57 + t90 * t98) * t230 - t158 * t19 - t160 * t20 - t424 * t33;
t38 = Icges(6,5) * t91 + Icges(6,6) * t90 - Icges(6,3) * t158;
t40 = Icges(6,4) * t91 + Icges(6,2) * t90 - Icges(6,6) * t158;
t42 = Icges(6,1) * t91 + Icges(6,4) * t90 - Icges(6,5) * t158;
t7 = t148 * t40 - t158 * t60 - t183 * t38 - t312 * t42 + t62 * t90 + t64 * t91;
t161 = t265 * t424;
t92 = qJD(5) * t311 + t161 * t266;
t93 = qJD(5) * t150 - t161 * t268;
t39 = Icges(6,5) * t93 + Icges(6,6) * t92 - Icges(6,3) * t160;
t41 = Icges(6,4) * t93 + Icges(6,2) * t92 - Icges(6,6) * t160;
t43 = Icges(6,1) * t93 + Icges(6,4) * t92 - Icges(6,5) * t160;
t8 = t148 * t41 - t158 * t61 - t183 * t39 - t312 * t43 + t63 * t90 + t65 * t91;
t434 = t279 * t401 + t7 * t408 + t8 * t410;
t433 = t437 * t264;
t432 = t437 * t265;
t431 = t436 * t265 + (t440 * t267 + t269 * t442) * qJD(2);
t430 = t436 * t264 + (t441 * t267 + t269 * t443) * qJD(2);
t429 = (-Icges(4,1) * t385 + t324 * t265 + t256 - t442) * t269 + (t265 * t438 - t440) * t267;
t260 = qJD(4) * t265;
t360 = pkin(3) * t368;
t228 = -t264 * t360 + t260;
t366 = qJD(4) * t264;
t229 = -t265 * t360 - t366;
t428 = t264 * t228 + t265 * t229 + t435;
t427 = (Icges(4,1) * t387 - t324 * t264 - t255 + t443) * t269 + (-t264 * t438 + t441) * t267;
t422 = -Icges(5,5) * t159 + Icges(5,6) * t158 - t264 * t439;
t421 = Icges(5,5) * t161 - Icges(5,6) * t160 + t265 * t439;
t420 = (-rSges(5,1) * t159 + rSges(5,2) * t158) * t264 + (-rSges(5,1) * t161 + rSges(5,2) * t160) * t265;
t237 = rSges(4,1) * t267 - rSges(4,3) * t269;
t419 = t372 * qJD(2) * t237;
t409 = t154 / 0.2e1;
t407 = t155 / 0.2e1;
t238 = rSges(3,1) * t267 + rSges(3,2) * t269;
t300 = qJD(2) * t238;
t365 = qJD(5) * t230;
t272 = t154 * (Icges(6,2) * t311 + t141 + t65) + t155 * (Icges(6,2) * t312 + t140 + t64) + t365 * (t100 + (-Icges(6,2) * t268 - t390) * t231);
t270 = qJD(2) ^ 2;
t10 = t150 * t41 - t160 * t61 - t185 * t39 - t311 * t43 + t63 * t92 + t65 * t93;
t21 = t150 * t62 - t185 * t60 - t311 * t64;
t22 = t150 * t63 - t185 * t61 - t311 * t65;
t34 = -t100 * t311 + t150 * t98 - t185 * t96;
t278 = (t100 * t93 + t150 * t56 - t160 * t96 - t185 * t55 - t311 * t57 + t92 * t98) * t230 - t158 * t21 - t160 * t22 - t424 * t34;
t9 = t150 * t40 - t160 * t60 - t185 * t38 - t311 * t42 + t62 * t92 + t64 * t93;
t411 = qJD(5) * t278 / 0.2e1 + t10 * t409 + t9 * t407;
t406 = -t372 * t368 / 0.2e1;
t405 = -t424 / 0.2e1;
t404 = -t269 / 0.2e1;
t403 = pkin(3) * t267;
t402 = pkin(3) * t269;
t394 = Icges(5,4) * t183;
t393 = Icges(5,4) * t185;
t239 = pkin(2) * t269 + qJ(3) * t267;
t226 = t239 * t264;
t227 = t239 * t265;
t377 = t264 * t226 + t265 * t227;
t367 = qJD(3) * t269;
t208 = qJD(2) * t239 - t367;
t240 = rSges(4,1) * t269 + rSges(4,3) * t267;
t376 = -t240 * qJD(2) - t208;
t375 = -t236 - t237;
t374 = -t239 - t240;
t373 = t250 + t260;
t362 = qJD(2) * qJD(3);
t361 = t270 * t402;
t358 = t156 * t371 + t157 * t370 + t267 * t362;
t356 = t269 * t362;
t355 = -t371 / 0.2e1;
t354 = t370 / 0.2e1;
t353 = t158 * t401;
t352 = t160 * t401;
t351 = qJD(5) * t405;
t350 = -t365 / 0.2e1;
t349 = t365 / 0.2e1;
t348 = -t236 - t403;
t347 = -t239 - t402;
t342 = qJD(2) * t376;
t341 = qJD(2) * t375;
t340 = t252 - t366;
t232 = pkin(3) * t386 + qJ(4) * t265;
t233 = pkin(3) * t384 - qJ(4) * t264;
t338 = t264 * t232 + t265 * t233 + t377;
t337 = t372 * t403;
t335 = -rSges(5,1) * t231 + rSges(5,2) * t230 + t348;
t333 = -pkin(4) * t231 - pkin(6) * t230 + t348;
t332 = -qJD(2) * t402 - t208;
t241 = rSges(3,1) * t269 - rSges(3,2) * t267;
t331 = rSges(6,1) * t268 - rSges(6,2) * t266;
t66 = -rSges(6,1) * t312 + rSges(6,2) * t148 - rSges(6,3) * t183;
t67 = -rSges(6,1) * t311 + rSges(6,2) * t150 - rSges(6,3) * t185;
t330 = t158 * t67 - t160 * t66;
t329 = -t266 * t62 + t268 * t64;
t328 = -t266 * t63 + t268 * t65;
t327 = t228 * t371 + t229 * t370 + t358;
t326 = t100 * t268 - t266 * t98;
t314 = (-Icges(5,5) * t183 + Icges(5,6) * t182) * t265 - (-Icges(5,5) * t185 + Icges(5,6) * t184) * t264;
t313 = t372 * t241;
t310 = t372 * t300;
t102 = rSges(6,3) * t230 + t231 * t331;
t309 = -t102 + t333;
t137 = rSges(5,1) * t213 + rSges(5,2) * t424;
t308 = -t137 + t332;
t307 = t226 * t371 + t227 * t370 + qJD(1) - t367;
t306 = qJD(2) * t335;
t305 = qJD(2) * t333;
t176 = -rSges(4,2) * t265 + t240 * t264;
t178 = rSges(4,2) * t264 + t240 * t265;
t58 = (t176 * t264 + t178 * t265) * qJD(2) + t307;
t304 = t58 * t237;
t139 = pkin(4) * t213 - pkin(6) * t424;
t59 = rSges(6,1) * t302 + rSges(6,2) * t303 - rSges(6,3) * t424;
t301 = -t139 + t332 - t59;
t44 = rSges(6,1) * t91 + rSges(6,2) * t90 - rSges(6,3) * t158;
t291 = -t102 * t158 - t230 * t44 + t424 * t66;
t45 = rSges(6,1) * t93 + rSges(6,2) * t92 - rSges(6,3) * t160;
t290 = t102 * t160 + t230 * t45 - t424 * t67;
t289 = -t361 + (-t137 - t208) * qJD(2);
t288 = -t361 + (-t139 - t208) * qJD(2);
t287 = t232 * t371 + t233 * t370 + t307;
t162 = Icges(5,4) * t182;
t113 = Icges(5,2) * t183 + Icges(5,6) * t265 + t162;
t163 = Icges(5,4) * t184;
t114 = Icges(5,2) * t185 - Icges(5,6) * t264 + t163;
t286 = (Icges(5,1) * t185 - t114 - t163) * t264 - (Icges(5,1) * t183 - t113 - t162) * t265;
t115 = Icges(5,1) * t182 + Icges(5,5) * t265 + t394;
t116 = Icges(5,1) * t184 - Icges(5,5) * t264 + t393;
t285 = (-Icges(5,2) * t184 + t116 + t393) * t264 - (-Icges(5,2) * t182 + t115 + t394) * t265;
t282 = (-Icges(6,5) * t266 - Icges(6,6) * t268) * t231 * t365 + (Icges(6,5) * t150 + Icges(6,6) * t311) * t154 + (Icges(6,5) * t148 + Icges(6,6) * t312) * t155;
t27 = t230 * t60 + t231 * t329;
t28 = t230 * t61 + t231 * t328;
t36 = t230 * t96 + t231 * t326;
t280 = (-t424 * t96 + t230 * t55 + t326 * t213 + (-t266 * t56 + t268 * t57 + (-t100 * t266 - t268 * t98) * qJD(5)) * t231) * t230 - t158 * t27 - t160 * t28 - t424 * t36;
t273 = (Icges(6,1) * t150 + t391 - t63) * t154 + (Icges(6,1) * t148 + t392 - t62) * t155 + ((-Icges(6,1) * t266 - t389) * t231 - t98) * t365;
t271 = (-Icges(6,3) * t184 - t185 * t315 + t328) * t154 + (-Icges(6,3) * t182 - t183 * t315 + t329) * t155 + (-Icges(6,3) * t231 + t230 * t315 + t326) * t365;
t253 = t265 * t367;
t251 = t264 * t367;
t249 = t265 * t356;
t248 = t264 * t356;
t138 = (-rSges(6,1) * t266 - rSges(6,2) * t268) * t231;
t133 = t265 * t341 + t252;
t132 = t264 * t341 + t250;
t131 = pkin(4) * t184 - pkin(6) * t185;
t130 = pkin(4) * t182 - pkin(6) * t183;
t121 = t310 * qJD(2);
t120 = rSges(5,1) * t184 + rSges(5,2) * t185 - rSges(5,3) * t264;
t119 = rSges(5,1) * t182 + rSges(5,2) * t183 + rSges(5,3) * t265;
t118 = -pkin(4) * t161 - pkin(6) * t160;
t117 = -pkin(4) * t159 - pkin(6) * t158;
t110 = -Icges(5,1) * t161 + Icges(5,4) * t160;
t109 = -Icges(5,1) * t159 + Icges(5,4) * t158;
t108 = -Icges(5,4) * t161 + Icges(5,2) * t160;
t107 = -Icges(5,4) * t159 + Icges(5,2) * t158;
t104 = t265 * t342 + t249;
t103 = t264 * t342 + t248;
t101 = -rSges(6,3) * t231 + t230 * t331;
t99 = -Icges(6,5) * t231 + t230 * t322;
t97 = -Icges(6,6) * t231 + t230 * t318;
t94 = qJD(2) * t313 + qJD(1);
t89 = rSges(6,1) * t150 + rSges(6,2) * t311;
t88 = rSges(6,1) * t148 + rSges(6,2) * t312;
t81 = -rSges(6,3) * t184 - t185 * t331;
t80 = -rSges(6,3) * t182 - t183 * t331;
t79 = -Icges(6,5) * t184 - t185 * t322;
t78 = -Icges(6,5) * t182 - t183 * t322;
t77 = -Icges(6,6) * t184 - t185 * t318;
t76 = -Icges(6,6) * t182 - t183 * t318;
t73 = t265 * t306 + t340;
t72 = t264 * t306 + t373;
t69 = t265 * t289 + t249;
t68 = t264 * t289 + t248;
t54 = -qJD(2) * t419 + t358;
t37 = (t119 * t264 + t120 * t265) * qJD(2) + t287;
t35 = qJD(2) * t420 + t327;
t32 = t102 * t155 + t265 * t305 - t365 * t66 + t340;
t31 = -t102 * t154 + t264 * t305 + t365 * t67 + t373;
t18 = t154 * t66 - t155 * t67 + (t130 * t264 + t131 * t265) * qJD(2) + t287;
t17 = qJD(5) * t291 + t155 * t59 + t265 * t288 + t249;
t16 = qJD(5) * t290 - t154 * t59 + t264 * t288 + t248;
t12 = t154 * t44 - t155 * t45 + t330 * qJD(5) + (t117 * t264 + t118 * t265) * qJD(2) + t327;
t11 = t154 * t28 + t155 * t27 + t36 * t365;
t6 = t154 * t22 + t155 * t21 + t34 * t365;
t5 = t154 * t20 + t155 * t19 + t33 * t365;
t4 = -t424 * t61 + t230 * t39 + t328 * t213 + (-t266 * t41 + t268 * t43 + (-t266 * t65 - t268 * t63) * qJD(5)) * t231;
t3 = -t424 * t60 + t230 * t38 + t329 * t213 + (-t266 * t40 + t268 * t42 + (-t266 * t64 - t268 * t62) * qJD(5)) * t231;
t1 = [-m(3) * t121 + m(4) * t54 + m(5) * t35 + m(6) * t12; ((t150 * t77 - t184 * t61 - t311 * t79) * t154 + (t150 * t76 - t184 * t60 - t311 * t78) * t155 + (-t22 * t184 - t21 * t182 + (t150 * t97 - t184 * t96 - t311 * t99) * t230 - t34 * t231) * qJD(5) - t271 * t185) * t410 + ((t148 * t77 - t182 * t61 - t312 * t79) * t154 + (t148 * t76 - t182 * t60 - t312 * t78) * t155 + (-t20 * t184 - t19 * t182 + (t148 * t97 - t182 * t96 - t312 * t99) * t230 - t33 * t231) * qJD(5) - t271 * t183) * t408 + ((-t27 * t182 - t28 * t184) * qJD(5) + ((-t266 * t77 + t268 * t79 - t61) * t154 + (-t266 * t76 + t268 * t78 - t60) * t155 - t36 * qJD(5)) * t231 + ((-t266 * t97 + t268 * t99 - t96) * t364 + t271) * t230) * t350 + t11 * t364 / 0.2e1 + t265 * t434 + (-t19 * t265 + t20 * t264) * t353 + (t264 * t8 - t265 * t7) * t407 + (t10 * t264 - t265 * t9) * t409 + t264 * t411 + (t264 * t4 - t265 * t3) * t349 + (t264 * t28 - t265 * t27) * t351 + (-t21 * t265 + t22 * t264) * t352 + (-t182 * t5 - t184 * t6) * t401 + (-t264 * (-t184 * t286 - t185 * t285 + t264 * t314) / 0.2e1 + t265 * (-t182 * t286 - t183 * t285 - t265 * t314) / 0.2e1) * t270 + (t12 * t338 + (t17 * t309 + t32 * t301 + t12 * (t131 + t67)) * t265 + (t16 * t309 + t31 * t301 + t12 * (t130 + t66)) * t264 - t32 * (t101 * t155 + t253) - t31 * (-t101 * t154 + t251) - (t32 * (-t102 * t182 - t230 * t80 + t231 * t66) + t31 * (t102 * t184 + t230 * t81 - t231 * t67)) * qJD(5) + (-t154 * t80 + t155 * t81 - (t182 * t67 - t184 * t66) * qJD(5) + (t118 + t45) * t265 + (t117 + t44) * t264 + t428) * t18) * m(6) + (t35 * t338 + (t35 * t120 + t308 * t73 + t335 * t69) * t265 + (t35 * t119 + t308 * t72 + t335 * t68) * t264 - t73 * t253 - t72 * t251 + (t420 + t428) * t37) * m(5) + (-t133 * t253 - t132 * t251 + t54 * t377 + (t104 * t375 + t133 * t376 + t54 * t178) * t265 + (t103 * t375 + t132 * t376 + t54 * t176) * t264 + (-t419 + t435) * t58) * m(4) + ((-t107 * t185 - t109 * t184 - t113 * t160 + t115 * t161 + t430 * t265 + (t427 - t432) * t354) * t265 + (t108 * t185 + t110 * t184 + t114 * t160 - t116 * t161 + t421 * t264 + t429 * t354 + (t422 - t431) * t265) * t264) * t371 + ((-t108 * t183 - t110 * t182 - t114 * t158 + t116 * t159 + t431 * t264 + (t429 - t433) * t355) * t264 + (t107 * t183 + t109 * t182 + t113 * t158 - t115 * t159 + t422 * t265 + t427 * t355 + (t421 - t430) * t264) * t265) * t370 + (-t121 * t313 - t94 * t310 + (t238 * t241 * t270 + t300 * t94) * t372) * m(3) + (t432 * t262 * t355 + t433 * t263 * t354 + (-(t31 * t264 + t32 * t265) * (-pkin(4) * t230 + pkin(6) * t231 + t347) + (t337 - t265 * (-pkin(4) * t185 - pkin(6) * t184) - t264 * (-pkin(4) * t183 - pkin(6) * t182)) * t18) * m(6) + (-(t72 * t264 + t73 * t265) * (-rSges(5,1) * t230 - rSges(5,2) * t231 + t347) + (t337 - t265 * (-rSges(5,1) * t185 + rSges(5,2) * t184) - t264 * (-rSges(5,1) * t183 + rSges(5,2) * t182)) * t37) * m(5) + (-(t133 * t374 - t265 * t304) * t265 - (t132 * t374 - t264 * t304) * t264) * m(4)) * qJD(2); 0.2e1 * (t12 * t404 + t18 * t406) * m(6) + 0.2e1 * (t35 * t404 + t37 * t406) * m(5) + 0.2e1 * (t404 * t54 + t406 * t58) * m(4) + 0.2e1 * (m(4) * (qJD(2) * t58 + t103 * t264 + t104 * t265) / 0.2e1 + m(5) * (qJD(2) * t37 + t264 * t68 + t265 * t69) / 0.2e1 + m(6) * (qJD(2) * t18 + t16 * t264 + t17 * t265) / 0.2e1) * t267; m(5) * (-t264 * t69 + t265 * t68) + m(6) * (t16 * t265 - t17 * t264); -t160 * t6 / 0.2e1 - t185 * t411 + (-t183 * t21 - t185 * t22 + t230 * t34) * t352 + (-t10 * t185 - t183 * t9 + t278) * t409 - t158 * t5 / 0.2e1 + t183 * t434 + (-t183 * t19 - t185 * t20 + t230 * t33) * t353 + (-t183 * t7 - t185 * t8 + t279) * t407 + t11 * t405 + t230 * (qJD(5) * t280 + t154 * t4 + t155 * t3) / 0.2e1 + (-t183 * t27 - t185 * t28 + t230 * t36) * t351 + (-t183 * t3 - t185 * t4 + t280) * t349 + (t150 * t272 - t185 * t282 - t273 * t311) * t410 + (t148 * t272 - t183 * t282 - t273 * t312) * t408 + (t282 * t230 + (-t266 * t272 + t273 * t268) * t231) * t350 + (t17 * (-t102 * t183 - t230 * t66) + t16 * (t102 * t185 + t230 * t67) + t12 * (t183 * t67 - t185 * t66) + (-t138 * t155 - t183 * t59 + t365 * t88 + t291) * t32 + (t138 * t154 + t185 * t59 - t365 * t89 + t290) * t31 + (-t154 * t88 + t155 * t89 + t183 * t45 - t185 * t44 + t330) * t18) * m(6);];
tauc = t1(:);
