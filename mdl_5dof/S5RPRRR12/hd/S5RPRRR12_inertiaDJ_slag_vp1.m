% Calculate time derivative of joint inertia matrix for
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR12_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR12_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR12_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR12_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:12:29
% EndTime: 2019-12-31 19:12:43
% DurationCPUTime: 7.17s
% Computational Cost: add. (16587->624), mult. (23498->880), div. (0->0), fcn. (22075->8), ass. (0->340)
t265 = sin(qJ(1));
t268 = cos(qJ(1));
t262 = qJ(3) + qJ(4);
t253 = cos(t262);
t252 = sin(t262);
t415 = Icges(5,4) * t252;
t302 = Icges(5,2) * t253 + t415;
t172 = Icges(5,6) * t268 + t265 * t302;
t414 = Icges(5,4) * t253;
t306 = Icges(5,1) * t252 + t414;
t174 = Icges(5,5) * t268 + t265 * t306;
t291 = t172 * t253 + t174 * t252;
t460 = t265 * t291;
t264 = sin(qJ(3));
t267 = cos(qJ(3));
t360 = qJD(3) * t267;
t363 = qJD(1) * t268;
t459 = t264 * t363 + t265 * t360;
t259 = qJD(3) + qJD(4);
t389 = t259 * t265;
t351 = t253 * t389;
t458 = t252 * t363 + t351;
t432 = rSges(6,3) + pkin(8);
t339 = t432 * t253;
t427 = pkin(4) * t252;
t457 = -t339 + t427;
t266 = cos(qJ(5));
t263 = sin(qJ(5));
t412 = Icges(6,4) * t266;
t301 = -Icges(6,2) * t263 + t412;
t395 = t252 * t259;
t413 = Icges(6,4) * t263;
t111 = -t301 * t395 + (Icges(6,6) * t259 + (-Icges(6,2) * t266 - t413) * qJD(5)) * t253;
t165 = Icges(6,6) * t252 + t253 * t301;
t305 = Icges(6,1) * t266 - t413;
t166 = Icges(6,5) * t252 + t253 * t305;
t456 = -t263 * t111 + (-t165 * t266 - t166 * t263) * qJD(5);
t380 = t266 * t268;
t385 = t263 * t265;
t204 = -t252 * t385 + t380;
t382 = t265 * t266;
t384 = t263 * t268;
t205 = t252 * t382 + t384;
t392 = t253 * t265;
t133 = Icges(6,5) * t205 + Icges(6,6) * t204 - Icges(6,3) * t392;
t135 = Icges(6,4) * t205 + Icges(6,2) * t204 - Icges(6,6) * t392;
t137 = Icges(6,1) * t205 + Icges(6,4) * t204 - Icges(6,5) * t392;
t50 = -t133 * t392 + t135 * t204 + t137 * t205;
t206 = t252 * t384 + t382;
t207 = -t252 * t380 + t385;
t391 = t253 * t268;
t134 = Icges(6,5) * t207 + Icges(6,6) * t206 + Icges(6,3) * t391;
t136 = Icges(6,4) * t207 + Icges(6,2) * t206 + Icges(6,6) * t391;
t138 = Icges(6,1) * t207 + Icges(6,4) * t206 + Icges(6,5) * t391;
t51 = -t134 * t392 + t136 * t204 + t138 * t205;
t38 = t265 * t51 + t268 * t50;
t298 = Icges(5,5) * t252 + Icges(5,6) * t253;
t170 = Icges(5,3) * t268 + t265 * t298;
t87 = t170 * t268 + t460;
t443 = -Icges(5,5) * t265 + t268 * t306;
t445 = -Icges(5,6) * t265 + t268 * t302;
t290 = -t252 * t443 - t253 * t445;
t282 = t290 * t265;
t447 = -Icges(5,3) * t265 + t268 * t298;
t88 = -t268 * t447 + t282;
t455 = -t265 * t88 - t268 * t87 - t38;
t212 = rSges(5,1) * t253 - rSges(5,2) * t252;
t387 = t259 * t268;
t454 = t212 * t387;
t417 = Icges(4,4) * t264;
t303 = Icges(4,2) * t267 + t417;
t188 = Icges(4,6) * t268 + t265 * t303;
t416 = Icges(4,4) * t267;
t307 = Icges(4,1) * t264 + t416;
t190 = Icges(4,5) * t268 + t265 * t307;
t289 = t188 * t267 + t190 * t264;
t453 = t268 * t289;
t452 = t268 * t291;
t317 = rSges(5,1) * t252 + rSges(5,2) * t253;
t451 = t268 * t317;
t318 = rSges(4,1) * t264 + rSges(4,2) * t267;
t283 = t268 * t318;
t314 = rSges(6,1) * t266 - rSges(6,2) * t263;
t113 = -t314 * t395 + (rSges(6,3) * t259 + (-rSges(6,1) * t263 - rSges(6,2) * t266) * qJD(5)) * t253;
t167 = rSges(6,3) * t252 + t253 * t314;
t450 = t265 * t113 + t167 * t363;
t210 = -Icges(5,2) * t252 + t414;
t323 = (t210 + t306) * t259;
t336 = t267 * t363;
t343 = rSges(4,1) * t459 + rSges(4,2) * t336;
t357 = -rSges(4,3) - pkin(1) - pkin(6);
t361 = qJD(3) * t264;
t370 = qJ(2) * t363 + qJD(2) * t265;
t105 = (-rSges(4,2) * t361 + qJD(1) * t357) * t265 + t343 + t370;
t423 = rSges(4,2) * t264;
t228 = rSges(4,1) * t267 - t423;
t251 = qJD(2) * t268;
t359 = qJD(3) * t268;
t106 = t251 + t228 * t359 + (t357 * t268 + (-qJ(2) - t318) * t265) * qJD(1);
t449 = -t105 * t268 + t106 * t265;
t321 = qJD(1) * t252 + qJD(5);
t350 = t253 * t387;
t448 = t265 * t321 - t350;
t299 = Icges(4,5) * t264 + Icges(4,6) * t267;
t446 = -Icges(4,3) * t265 + t268 * t299;
t444 = -Icges(4,6) * t265 + t268 * t303;
t442 = -Icges(4,5) * t265 + t268 * t307;
t182 = t302 * t259;
t209 = Icges(5,5) * t253 - Icges(5,6) * t252;
t211 = Icges(5,1) * t253 - t415;
t440 = t252 * t323 - t253 * (t211 * t259 - t182) + qJD(1) * t209;
t439 = 2 * m(4);
t438 = 2 * m(5);
t437 = 2 * m(6);
t260 = t265 ^ 2;
t261 = t268 ^ 2;
t436 = t265 / 0.2e1;
t435 = t268 / 0.2e1;
t434 = rSges(3,2) - pkin(1);
t433 = -rSges(5,3) - pkin(1);
t431 = m(4) * t228;
t430 = pkin(1) * t265;
t429 = pkin(3) * t264;
t428 = pkin(3) * t267;
t426 = pkin(4) * t253;
t258 = t268 * pkin(1);
t425 = t268 * pkin(6);
t297 = Icges(6,5) * t266 - Icges(6,6) * t263;
t110 = -t297 * t395 + (Icges(6,3) * t259 + (-Icges(6,5) * t263 - Icges(6,6) * t266) * qJD(5)) * t253;
t112 = -t305 * t395 + (Icges(6,5) * t259 + (-Icges(6,1) * t263 - t412) * qJD(5)) * t253;
t164 = Icges(6,3) * t252 + t253 * t297;
t388 = t259 * t266;
t393 = t253 * t259;
t403 = t165 * t263;
t275 = t253 * t266 * t112 + t164 * t393 + t395 * t403 + (-t166 * t388 + t110) * t252;
t83 = t164 * t252 + (t166 * t266 - t403) * t253;
t424 = (t253 * t456 + t275) * t252 + t83 * t393;
t296 = t135 * t263 - t137 * t266;
t322 = qJD(5) * t252 + qJD(1);
t119 = -t322 * t382 + (-t268 * t321 - t351) * t263;
t120 = t321 * t380 + (t253 * t388 - t263 * t322) * t265;
t335 = t253 * t363;
t349 = t252 * t389;
t278 = -t335 + t349;
t68 = Icges(6,5) * t120 + Icges(6,6) * t119 + Icges(6,3) * t278;
t70 = Icges(6,4) * t120 + Icges(6,2) * t119 + Icges(6,6) * t278;
t72 = Icges(6,1) * t120 + Icges(6,4) * t119 + Icges(6,5) * t278;
t20 = (t259 * t296 + t68) * t252 + (t133 * t259 - t263 * t70 + t266 * t72 + (-t135 * t266 - t137 * t263) * qJD(5)) * t253;
t422 = t20 * t268;
t295 = t136 * t263 - t138 * t266;
t286 = t322 * t268;
t117 = -t263 * t448 + t266 * t286;
t118 = t263 * t286 + t266 * t448;
t348 = t252 * t387;
t364 = qJD(1) * t265;
t277 = -t253 * t364 - t348;
t67 = Icges(6,5) * t118 + Icges(6,6) * t117 + Icges(6,3) * t277;
t69 = Icges(6,4) * t118 + Icges(6,2) * t117 + Icges(6,6) * t277;
t71 = Icges(6,1) * t118 + Icges(6,4) * t117 + Icges(6,5) * t277;
t21 = (t259 * t295 + t67) * t252 + (t134 * t259 - t263 * t69 + t266 * t71 + (-t136 * t266 - t138 * t263) * qJD(5)) * t253;
t421 = t21 * t265;
t420 = t265 * rSges(4,3);
t257 = t268 * rSges(4,3);
t256 = t268 * rSges(5,3);
t316 = rSges(6,1) * t118 + rSges(6,2) * t117;
t73 = rSges(6,3) * t277 + t316;
t419 = (t73 + t277 * pkin(8) + (t252 * t364 - t350) * pkin(4)) * t268;
t344 = pkin(4) * t458 + pkin(8) * t349;
t347 = t120 * rSges(6,1) + t119 * rSges(6,2) + rSges(6,3) * t349;
t74 = -rSges(6,3) * t335 + t347;
t418 = pkin(8) * t335 - t344 - t74;
t184 = t317 * t259;
t400 = t184 * t268;
t399 = t188 * t264;
t398 = t444 * t264;
t397 = t190 * t267;
t396 = t442 * t267;
t394 = t252 * t265;
t383 = t264 * t265;
t162 = t265 * t167;
t381 = t265 * t267;
t315 = -rSges(6,1) * t207 - rSges(6,2) * t206;
t140 = rSges(6,3) * t391 - t315;
t131 = t268 * t140;
t320 = pkin(8) * t253 - t427;
t199 = t320 * t268;
t379 = t268 * t199 + t131;
t374 = t205 * rSges(6,1) + t204 * rSges(6,2);
t139 = -rSges(6,3) * t392 + t374;
t231 = pkin(4) * t394;
t378 = pkin(8) * t392 - t139 - t231;
t377 = -t140 - t199;
t213 = pkin(8) * t252 + t426;
t376 = -t167 - t213;
t177 = rSges(5,1) * t394 + rSges(5,2) * t392 + t256;
t245 = pkin(3) * t383;
t269 = -pkin(7) - pkin(6);
t203 = t245 + (-pkin(6) - t269) * t268;
t375 = -t177 - t203;
t142 = t265 * t213 + t162;
t246 = pkin(3) * t381;
t354 = pkin(3) * t359;
t373 = qJD(1) * t246 + t264 * t354;
t362 = qJD(1) * t269;
t372 = -t267 * t354 - t268 * t362;
t371 = t265 * t269 + t268 * t429;
t369 = t265 * qJ(2) + t258;
t368 = t260 + t261;
t367 = qJD(1) * t170;
t186 = Icges(4,3) * t268 + t265 * t299;
t366 = qJD(1) * t186;
t356 = rSges(5,2) * t395;
t355 = pkin(3) * t361;
t61 = t133 * t252 - t253 * t296;
t79 = -t164 * t392 + t165 * t204 + t166 * t205;
t353 = t79 / 0.2e1 + t61 / 0.2e1;
t62 = t134 * t252 - t253 * t295;
t80 = t164 * t391 + t165 * t206 + t166 * t207;
t352 = -t80 / 0.2e1 - t62 / 0.2e1;
t346 = -t203 + t378;
t345 = rSges(5,1) * t458 + rSges(5,2) * t335;
t342 = pkin(3) * t459 + t265 * t362;
t341 = t251 - t372;
t192 = rSges(4,1) * t383 + rSges(4,2) * t381 + t257;
t255 = t268 * qJ(2);
t340 = t255 + t371;
t158 = t167 * t364;
t333 = -t364 / 0.2e1;
t332 = t363 / 0.2e1;
t331 = -qJ(2) - t429;
t330 = t184 * t368;
t219 = t318 * qJD(3);
t329 = t219 * t368;
t124 = qJD(1) * t172 - t210 * t387;
t328 = t259 * t443 - t124;
t125 = qJD(1) * t445 + t210 * t389;
t327 = t174 * t259 + t125;
t126 = qJD(1) * t174 - t211 * t387;
t326 = -t259 * t445 - t126;
t127 = qJD(1) * t443 + t211 * t389;
t325 = -t172 * t259 + t127;
t185 = t320 * t259;
t76 = t265 * t185 + t213 * t363 + t450;
t16 = t119 * t135 + t120 * t137 + t133 * t278 + t204 * t70 + t205 * t72 - t392 * t68;
t17 = t119 * t136 + t120 * t138 + t134 * t278 + t204 * t69 + t205 * t71 - t392 * t67;
t313 = t265 * t50 - t268 * t51;
t10 = -qJD(1) * t313 + t16 * t268 + t17 * t265;
t122 = -t209 * t387 + t367;
t123 = qJD(1) * t447 + t209 * t389;
t52 = t133 * t391 + t135 * t206 + t137 * t207;
t53 = t134 * t391 + t136 * t206 + t138 * t207;
t39 = t265 * t53 + t268 * t52;
t89 = t170 * t265 - t452;
t14 = t117 * t135 + t118 * t137 + t133 * t277 + t206 * t70 + t207 * t72 + t391 * t68;
t15 = t117 * t136 + t118 * t138 + t134 * t277 + t206 * t69 + t207 * t71 + t391 * t67;
t312 = t265 * t52 - t268 * t53;
t9 = -qJD(1) * t312 + t14 * t268 + t15 * t265;
t90 = -t265 * t447 - t268 * t290;
t319 = t39 * t363 + (t89 * t363 + t10 + (t268 * t123 + (t88 + t452) * qJD(1)) * t268) * t268 + (t9 + t90 * t363 + (t265 * t122 + (-t89 + t282) * qJD(1)) * t265 + ((-t393 * t443 + t395 * t445 + t123) * t265 + (t172 * t395 - t174 * t393 + t122 + t367) * t268 + ((t124 + t328) * t265 + (-t125 + t327) * t268) * t253 + ((t126 + t326) * t265 + (-t127 + t325) * t268) * t252 + (t90 - t87 + t460 + (-t170 + t290) * t268) * qJD(1)) * t268) * t265;
t311 = t62 * t265 + t61 * t268;
t310 = t61 * t265 - t62 * t268;
t309 = t342 + t370;
t308 = Icges(4,1) * t267 - t417;
t304 = -Icges(4,2) * t264 + t416;
t300 = Icges(4,5) * t267 - Icges(4,6) * t264;
t294 = t139 * t268 + t140 * t265;
t288 = -t264 * t442 - t267 * t444;
t287 = t210 * t253 + t211 * t252;
t285 = -t268 * t269 + t245 + t369;
t284 = rSges(3,3) * t268 + t265 * t434;
t75 = t158 + t213 * t364 + (-t113 - t185) * t268;
t281 = t288 * t265;
t280 = qJD(3) * t308;
t279 = qJD(3) * t304;
t276 = qJD(1) * t287 - t298 * t259;
t26 = t252 * t79 - t253 * t313;
t27 = t252 * t80 - t253 * t312;
t31 = t110 * t391 + t111 * t206 + t112 * t207 + t117 * t165 + t118 * t166 + t164 * t277;
t3 = (t259 * t312 + t31) * t252 + (-qJD(1) * t39 - t14 * t265 + t15 * t268 + t259 * t80) * t253;
t32 = -t110 * t392 + t111 * t204 + t112 * t205 + t119 * t165 + t120 * t166 + t164 * t278;
t4 = (t259 * t313 + t32) * t252 + (-qJD(1) * t38 - t16 * t265 + t17 * t268 + t259 * t79) * t253;
t274 = -t10 * t392 / 0.2e1 + t3 * t436 + t4 * t435 + t252 * (-qJD(1) * t310 + t421 + t422) / 0.2e1 + t27 * t332 - t39 * t348 / 0.2e1 + t311 * t393 / 0.2e1 + t9 * t391 / 0.2e1 + (t349 / 0.2e1 - t335 / 0.2e1) * t38 + (t253 * t39 + t26) * t333;
t273 = t455 * t364 + t319;
t144 = t265 * t433 + t340 + t451;
t145 = t285 + t177;
t85 = (qJD(1) * t433 - t356) * t265 + t309 + t345;
t86 = t454 + (t433 * t268 + (-t317 + t331) * t265) * qJD(1) + t341;
t272 = m(5) * (t265 * t86 - t268 * t85 + (t144 * t268 + t145 * t265) * qJD(1));
t108 = t212 * t364 + t373 + t400;
t240 = pkin(3) * t336;
t109 = t212 * t363 + t240 + (-t184 - t355) * t265;
t168 = t212 * t265 + t246;
t169 = (-t212 - t428) * t268;
t271 = m(5) * (-t108 * t268 + t109 * t265 + (t168 * t268 + t169 * t265) * qJD(1));
t270 = t422 / 0.2e1 + t421 / 0.2e1 + (t252 * t328 - t253 * t326 + t276 * t265 + t268 * t440 + t31) * t436 + (-t252 * t327 + t253 * t325 - t265 * t440 + t276 * t268 + t32) * t435 + (-t172 * t252 + t174 * t253 + t209 * t268 + t265 * t287 + t61 + t79) * t333 + (t209 * t265 + t252 * t445 - t253 * t443 - t268 * t287 + t62 + t80) * t332;
t202 = -pkin(6) * t265 - t371;
t197 = -rSges(3,2) * t268 + rSges(3,3) * t265 + t369;
t196 = t255 + t284;
t193 = t420 - t283;
t180 = t268 * t202;
t178 = rSges(5,3) * t265 - t451;
t163 = t268 * t178;
t161 = t251 + (t434 * t268 + (-rSges(3,3) - qJ(2)) * t265) * qJD(1);
t160 = qJD(1) * t284 + t370;
t157 = pkin(6) * t364 + t342;
t156 = t192 + t369 + t425;
t155 = t265 * t357 + t255 + t283;
t154 = t268 * ((t245 - t425) * qJD(1) + t372);
t147 = t300 * t265 * qJD(3) + qJD(1) * t446;
t146 = -t300 * t359 + t366;
t143 = t376 * t268;
t130 = (t376 - t428) * t268;
t129 = t246 + t142;
t128 = (-rSges(5,3) * qJD(1) - t356) * t265 + t345;
t121 = -t177 * t265 + t163;
t114 = t268 * (-t454 + (t265 * t317 + t256) * qJD(1));
t98 = -t265 * t446 - t268 * t288;
t97 = t186 * t265 - t453;
t96 = -t268 * t446 + t281;
t95 = t186 * t268 + t289 * t265;
t94 = -t140 * t252 + t167 * t391;
t93 = t139 * t252 + t162 * t253;
t92 = -t265 * t339 + t231 + t285 + t374;
t91 = t268 * t457 + t315 + t340 - t430;
t84 = t265 * t375 + t163 + t180;
t81 = t294 * t253;
t65 = t265 * t378 + t379;
t64 = -t265 * t355 + t240 + t76;
t63 = t75 + t373;
t60 = -t128 * t265 + t114 + (-t177 * t268 - t178 * t265) * qJD(1);
t55 = t265 * t346 + t180 + t379;
t46 = (t252 * t432 + t426) * t387 + (-t258 + (t331 - t457) * t265) * qJD(1) - t316 + t341;
t45 = (-t268 * t339 - t430) * qJD(1) + t309 + t344 + t347;
t44 = t114 + t154 + (-t128 - t157) * t265 + (t375 * t268 + (-t178 - t202) * t265) * qJD(1);
t43 = (-t162 * t259 + t74) * t252 + (t139 * t259 + t450) * t253;
t42 = (-t167 * t387 - t73) * t252 + (t113 * t268 - t140 * t259 - t158) * t253;
t28 = t294 * t395 + (-t265 * t73 - t268 * t74 + (t139 * t265 - t131) * qJD(1)) * t253;
t25 = t418 * t265 + (t265 * t377 + t268 * t378) * qJD(1) + t419;
t22 = t154 + (-t157 + t418) * t265 + (t346 * t268 + (-t202 + t377) * t265) * qJD(1) + t419;
t1 = [-t267 * t279 - t264 * t280 + t303 * t361 - t307 * t360 + t275 + (t45 * t92 + t46 * t91) * t437 + (t144 * t86 + t145 * t85) * t438 + (t105 * t156 + t106 * t155) * t439 + 0.2e1 * m(3) * (t160 * t197 + t161 * t196) + t252 * t182 - t211 * t395 + (-t323 + t456) * t253; m(6) * (t265 * t46 - t268 * t45 + (t265 * t92 + t268 * t91) * qJD(1)) + t272 + m(4) * ((t155 * t268 + t156 * t265) * qJD(1) + t449) + m(3) * (-t160 * t268 + t161 * t265 + (t196 * t268 + t197 * t265) * qJD(1)); 0; m(4) * (t449 * t228 - (t155 * t265 - t156 * t268) * t219) - (t260 / 0.2e1 + t261 / 0.2e1) * t299 * qJD(3) + t270 + m(6) * (t129 * t46 + t130 * t45 + t63 * t92 + t64 * t91) + m(5) * (t108 * t145 + t109 * t144 + t168 * t86 + t169 * t85) + ((t399 / 0.2e1 - t397 / 0.2e1 + t156 * t431) * t265 + (t155 * t431 + t398 / 0.2e1 - t396 / 0.2e1) * t268) * qJD(1) + (-qJD(3) * t289 - (qJD(1) * t444 + t265 * t279) * t264 + (qJD(1) * t442 + t265 * t280) * t267) * t435 + (-qJD(3) * t288 - (qJD(1) * t188 - t304 * t359) * t264 + (qJD(1) * t190 - t308 * t359) * t267) * t436; t271 + m(6) * (t265 * t64 - t268 * t63 + (t129 * t268 + t130 * t265) * qJD(1)) - m(4) * t329; (t129 * t64 + t130 * t63 + t55 * t22) * t437 + (t108 * t169 + t109 * t168 + t44 * t84) * t438 + ((-t192 * t265 + t193 * t268) * (-t265 * t343 + (-t228 * t261 + t260 * t423) * qJD(3) + ((-t192 + t257) * t268 + (-t193 + t283 + t420) * t265) * qJD(1)) - t228 * t329) * t439 + t268 * ((t268 * t147 + (t96 + t453) * qJD(1)) * t268 + (-t95 * qJD(1) + (-t360 * t442 + t361 * t444) * t265 + (t146 + (t397 - t399) * qJD(3) + (-t186 + t288) * qJD(1)) * t268) * t265) + (t265 * t98 + t268 * t97) * t363 + t265 * ((t265 * t146 + (-t97 + t281) * qJD(1)) * t265 + (t98 * qJD(1) + (t188 * t361 - t190 * t360 + t366) * t268 + (t147 + (t396 - t398) * qJD(3) + t289 * qJD(1)) * t265) * t268) + t319 + (-t265 * t96 - t268 * t95 + t455) * t364; m(6) * (t142 * t46 + t143 * t45 + t75 * t92 + t76 * t91) + t212 * t272 - m(5) * (t144 * t265 - t145 * t268) * t184 + t270; m(6) * (t76 * t265 - t75 * t268 + (t142 * t268 + t143 * t265) * qJD(1)) - m(5) * t330; m(6) * (t129 * t76 + t130 * t75 + t142 * t64 + t143 * t63 + t65 * t22 + t25 * t55) + m(5) * (-t168 * t184 * t265 + t121 * t44 + t169 * t400 + t60 * t84) + t212 * t271 + t273; (t121 * t60 - t212 * t330) * t438 + (t142 * t76 + t143 * t75 + t65 * t25) * t437 + t273; m(6) * (t42 * t91 + t43 * t92 + t45 * t93 + t46 * t94) + (t265 * t353 + t268 * t352) * t395 + ((t31 / 0.2e1 + t21 / 0.2e1) * t268 + (-t32 / 0.2e1 - t20 / 0.2e1) * t265 + (t265 * t352 - t268 * t353) * qJD(1)) * t253 + t424; m(6) * (t265 * t42 - t268 * t43 + (t265 * t93 + t268 * t94) * qJD(1)); t274 + m(6) * (t129 * t42 + t130 * t43 - t22 * t81 + t28 * t55 + t63 * t93 + t64 * t94); t274 + m(6) * (t142 * t42 + t143 * t43 - t25 * t81 + t28 * t65 + t75 * t93 + t76 * t94); (-t28 * t81 + t42 * t94 + t43 * t93) * t437 + ((t252 * t310 + t265 * t26 - t268 * t27) * t259 + t424) * t252 + (-t265 * t4 + t268 * t3 - t310 * t393 + (-t20 * t265 + t21 * t268 + t259 * t83) * t252 + (-t252 * t311 - t268 * t26 - t265 * t27) * qJD(1)) * t253;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
