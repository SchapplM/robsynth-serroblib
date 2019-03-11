% Calculate joint inertia matrix for
% S6RRRRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR15_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR15_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR15_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR15_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:34:19
% EndTime: 2019-03-10 00:34:40
% DurationCPUTime: 8.81s
% Computational Cost: add. (51322->723), mult. (142112->982), div. (0->0), fcn. (186583->14), ass. (0->331)
t308 = cos(pkin(6));
t313 = sin(qJ(1));
t315 = cos(qJ(2));
t380 = t313 * t315;
t312 = sin(qJ(2));
t316 = cos(qJ(1));
t382 = t312 * t316;
t294 = -t308 * t380 - t382;
t307 = sin(pkin(6));
t389 = cos(pkin(7));
t352 = t307 * t389;
t388 = sin(pkin(7));
t275 = -t294 * t388 + t313 * t352;
t379 = t315 * t316;
t381 = t313 * t312;
t292 = t308 * t379 - t381;
t274 = -t292 * t388 - t316 * t352;
t383 = t307 * t316;
t293 = t308 * t382 + t380;
t311 = sin(qJ(3));
t351 = t307 * t388;
t392 = cos(qJ(3));
t255 = t293 * t392 + (t292 * t389 - t316 * t351) * t311;
t310 = sin(qJ(4));
t391 = cos(qJ(4));
t229 = t255 * t391 + t274 * t310;
t295 = -t308 * t381 + t379;
t257 = t295 * t392 + (t294 * t389 + t313 * t351) * t311;
t231 = t257 * t391 + t275 * t310;
t273 = t308 * t388 * t311 + (t311 * t315 * t389 + t312 * t392) * t307;
t291 = t308 * t389 - t315 * t351;
t249 = t273 * t391 + t291 * t310;
t228 = t255 * t310 - t274 * t391;
t333 = t392 * t388;
t329 = t307 * t333;
t334 = t389 * t392;
t254 = -t292 * t334 + t293 * t311 + t316 * t329;
t309 = sin(qJ(6));
t314 = cos(qJ(6));
t179 = t228 * t314 - t254 * t309;
t180 = t228 * t309 + t254 * t314;
t117 = Icges(7,5) * t180 + Icges(7,6) * t179 + Icges(7,3) * t229;
t119 = Icges(7,4) * t180 + Icges(7,2) * t179 + Icges(7,6) * t229;
t121 = Icges(7,1) * t180 + Icges(7,4) * t179 + Icges(7,5) * t229;
t36 = t117 * t229 + t119 * t179 + t121 * t180;
t230 = t257 * t310 - t275 * t391;
t256 = -t294 * t334 + t295 * t311 - t313 * t329;
t181 = t230 * t314 - t256 * t309;
t182 = t230 * t309 + t256 * t314;
t118 = Icges(7,5) * t182 + Icges(7,6) * t181 + Icges(7,3) * t231;
t120 = Icges(7,4) * t182 + Icges(7,2) * t181 + Icges(7,6) * t231;
t122 = Icges(7,1) * t182 + Icges(7,4) * t181 + Icges(7,5) * t231;
t37 = t118 * t229 + t120 * t179 + t122 * t180;
t248 = t273 * t310 - t291 * t391;
t384 = t307 * t315;
t386 = t307 * t312;
t272 = -t308 * t333 + t311 * t386 - t334 * t384;
t219 = t248 * t314 - t272 * t309;
t220 = t248 * t309 + t272 * t314;
t138 = Icges(7,5) * t220 + Icges(7,6) * t219 + Icges(7,3) * t249;
t139 = Icges(7,4) * t220 + Icges(7,2) * t219 + Icges(7,6) * t249;
t140 = Icges(7,1) * t220 + Icges(7,4) * t219 + Icges(7,5) * t249;
t54 = t138 * t229 + t139 * t179 + t140 * t180;
t1 = t229 * t36 + t231 * t37 + t249 * t54;
t398 = t1 / 0.2e1;
t38 = t117 * t231 + t119 * t181 + t121 * t182;
t39 = t118 * t231 + t120 * t181 + t122 * t182;
t55 = t138 * t231 + t139 * t181 + t140 * t182;
t2 = t229 * t38 + t231 * t39 + t249 * t55;
t397 = t2 / 0.2e1;
t46 = t117 * t249 + t119 * t219 + t121 * t220;
t47 = t118 * t249 + t120 * t219 + t122 * t220;
t63 = t249 * t138 + t219 * t139 + t220 * t140;
t57 = t63 * t249;
t9 = t46 * t229 + t47 * t231 + t57;
t396 = t9 / 0.2e1;
t395 = t229 / 0.2e1;
t394 = t231 / 0.2e1;
t393 = t249 / 0.2e1;
t390 = pkin(5) * t254;
t261 = Icges(3,5) * t293 + Icges(3,6) * t292 - Icges(3,3) * t383;
t387 = t261 * t316;
t385 = t307 * t313;
t331 = -rSges(7,1) * t180 - rSges(7,2) * t179;
t123 = rSges(7,3) * t229 - t331;
t378 = pkin(12) * t229 + t123 + t390;
t124 = t182 * rSges(7,1) + t181 * rSges(7,2) + t231 * rSges(7,3);
t377 = t256 * pkin(5) + t231 * pkin(12) + t124;
t141 = rSges(7,1) * t220 + rSges(7,2) * t219 + rSges(7,3) * t249;
t376 = pkin(5) * t272 + pkin(12) * t249 + t141;
t332 = -rSges(6,1) * t254 - rSges(6,3) * t228;
t154 = -rSges(6,2) * t229 - t332;
t221 = t228 * qJ(5);
t172 = pkin(4) * t229 + t221;
t375 = -t154 - t172;
t155 = t256 * rSges(6,1) - t231 * rSges(6,2) + t230 * rSges(6,3);
t173 = t231 * pkin(4) + t230 * qJ(5);
t374 = -t155 - t173;
t156 = rSges(5,1) * t229 - rSges(5,2) * t228 + rSges(5,3) * t254;
t212 = pkin(3) * t255 + t254 * pkin(11);
t373 = -t156 - t212;
t204 = t275 * t212;
t372 = t275 * t172 + t204;
t213 = t257 * pkin(3) + t256 * pkin(11);
t205 = t291 * t213;
t371 = t291 * t173 + t205;
t193 = rSges(6,1) * t272 - rSges(6,2) * t249 + rSges(6,3) * t248;
t211 = pkin(4) * t249 + qJ(5) * t248;
t370 = -t193 - t211;
t194 = rSges(5,1) * t249 - rSges(5,2) * t248 + rSges(5,3) * t272;
t240 = pkin(3) * t273 + pkin(11) * t272;
t369 = -t194 - t240;
t218 = t274 * t240;
t368 = t274 * t211 + t218;
t260 = t295 * pkin(2) + t275 * pkin(10);
t258 = t308 * t260;
t367 = t308 * t213 + t258;
t259 = pkin(2) * t293 + t274 * pkin(10);
t366 = t259 * t385 + t260 * t383;
t365 = t316 * pkin(1) + pkin(9) * t385;
t364 = t46 / 0.2e1 + t54 / 0.2e1;
t363 = t47 / 0.2e1 + t55 / 0.2e1;
t362 = -t172 - t378;
t361 = -t173 - t377;
t360 = -t211 - t376;
t359 = -t212 + t375;
t186 = Icges(5,5) * t249 - Icges(5,6) * t248 + Icges(5,3) * t272;
t188 = Icges(5,4) * t249 - Icges(5,2) * t248 + Icges(5,6) * t272;
t190 = Icges(5,1) * t249 - Icges(5,4) * t248 + Icges(5,5) * t272;
t105 = t272 * t186 - t248 * t188 + t249 * t190;
t358 = t308 * t173 + t367;
t185 = Icges(6,5) * t272 - Icges(6,6) * t249 + Icges(6,3) * t248;
t187 = Icges(6,4) * t272 - Icges(6,2) * t249 + Icges(6,6) * t248;
t189 = Icges(6,1) * t272 - Icges(6,4) * t249 + Icges(6,5) * t248;
t104 = t248 * t185 - t249 * t187 + t272 * t189;
t357 = -t240 + t370;
t233 = Icges(4,5) * t273 - Icges(4,6) * t272 + Icges(4,3) * t291;
t234 = Icges(4,4) * t273 - Icges(4,2) * t272 + Icges(4,6) * t291;
t235 = Icges(4,1) * t273 - Icges(4,4) * t272 + Icges(4,5) * t291;
t132 = t291 * t233 - t272 * t234 + t273 * t235;
t157 = t231 * rSges(5,1) - t230 * rSges(5,2) + t256 * rSges(5,3);
t203 = t257 * rSges(4,1) - t256 * rSges(4,2) + t275 * rSges(4,3);
t282 = Icges(3,3) * t308 + (Icges(3,5) * t312 + Icges(3,6) * t315) * t307;
t283 = Icges(3,6) * t308 + (Icges(3,4) * t312 + Icges(3,2) * t315) * t307;
t284 = Icges(3,5) * t308 + (Icges(3,1) * t312 + Icges(3,4) * t315) * t307;
t356 = t308 * t282 + t283 * t384 + t284 * t386;
t268 = t295 * rSges(3,1) + t294 * rSges(3,2) + rSges(3,3) * t385;
t355 = -t313 * pkin(1) + pkin(9) * t383;
t236 = rSges(4,1) * t273 - rSges(4,2) * t272 + rSges(4,3) * t291;
t278 = pkin(2) * t386 + pkin(10) * t291;
t350 = t307 * (-t236 - t278);
t349 = -t212 + t362;
t348 = -t240 + t360;
t347 = t212 * t385 + t213 * t383 + t366;
t142 = Icges(6,5) * t254 - Icges(6,6) * t229 + Icges(6,3) * t228;
t146 = Icges(6,4) * t254 - Icges(6,2) * t229 + Icges(6,6) * t228;
t150 = Icges(6,1) * t254 - Icges(6,4) * t229 + Icges(6,5) * t228;
t67 = t142 * t228 - t146 * t229 + t150 * t254;
t143 = Icges(6,5) * t256 - Icges(6,6) * t231 + Icges(6,3) * t230;
t147 = Icges(6,4) * t256 - Icges(6,2) * t231 + Icges(6,6) * t230;
t151 = Icges(6,1) * t256 - Icges(6,4) * t231 + Icges(6,5) * t230;
t68 = t143 * t228 - t147 * t229 + t151 * t254;
t87 = t185 * t228 - t187 * t229 + t189 * t254;
t13 = t254 * t67 + t256 * t68 + t272 * t87;
t144 = Icges(5,5) * t229 - Icges(5,6) * t228 + Icges(5,3) * t254;
t148 = Icges(5,4) * t229 - Icges(5,2) * t228 + Icges(5,6) * t254;
t152 = Icges(5,1) * t229 - Icges(5,4) * t228 + Icges(5,5) * t254;
t71 = t144 * t254 - t148 * t228 + t152 * t229;
t145 = Icges(5,5) * t231 - Icges(5,6) * t230 + Icges(5,3) * t256;
t149 = Icges(5,4) * t231 - Icges(5,2) * t230 + Icges(5,6) * t256;
t153 = Icges(5,1) * t231 - Icges(5,4) * t230 + Icges(5,5) * t256;
t72 = t145 * t254 - t149 * t228 + t153 * t229;
t89 = t186 * t254 - t188 * t228 + t190 * t229;
t15 = t254 * t71 + t256 * t72 + t272 * t89;
t3 = t254 * t36 + t256 * t37 + t272 * t54;
t346 = t3 / 0.2e1 + t13 / 0.2e1 + t15 / 0.2e1;
t69 = t142 * t230 - t146 * t231 + t150 * t256;
t70 = t143 * t230 - t147 * t231 + t151 * t256;
t88 = t185 * t230 - t187 * t231 + t189 * t256;
t14 = t254 * t69 + t256 * t70 + t272 * t88;
t73 = t144 * t256 - t148 * t230 + t152 * t231;
t74 = t145 * t256 - t149 * t230 + t153 * t231;
t90 = t186 * t256 - t188 * t230 + t190 * t231;
t16 = t254 * t73 + t256 * t74 + t272 * t90;
t4 = t254 * t38 + t256 * t39 + t272 * t55;
t345 = t4 / 0.2e1 + t16 / 0.2e1 + t14 / 0.2e1;
t17 = t274 * t67 + t275 * t68 + t291 * t87;
t19 = t274 * t71 + t275 * t72 + t291 * t89;
t5 = t274 * t36 + t275 * t37 + t291 * t54;
t344 = t5 / 0.2e1 + t19 / 0.2e1 + t17 / 0.2e1;
t18 = t274 * t69 + t275 * t70 + t291 * t88;
t20 = t274 * t73 + t275 * t74 + t291 * t90;
t6 = t274 * t38 + t275 * t39 + t291 * t55;
t343 = t6 / 0.2e1 + t20 / 0.2e1 + t18 / 0.2e1;
t21 = t87 * t308 + (t313 * t68 - t316 * t67) * t307;
t23 = t89 * t308 + (t313 * t72 - t316 * t71) * t307;
t7 = t54 * t308 + (t313 * t37 - t316 * t36) * t307;
t342 = t7 / 0.2e1 + t21 / 0.2e1 + t23 / 0.2e1;
t22 = t88 * t308 + (t313 * t70 - t316 * t69) * t307;
t24 = t90 * t308 + (t313 * t74 - t316 * t73) * t307;
t8 = t55 * t308 + (t313 * t39 - t316 * t38) * t307;
t341 = t8 / 0.2e1 + t24 / 0.2e1 + t22 / 0.2e1;
t58 = t63 * t272;
t10 = t46 * t254 + t47 * t256 + t58;
t77 = t142 * t248 - t146 * t249 + t150 * t272;
t78 = t143 * t248 - t147 * t249 + t151 * t272;
t95 = t104 * t272;
t25 = t77 * t254 + t78 * t256 + t95;
t79 = t144 * t272 - t148 * t248 + t152 * t249;
t80 = t145 * t272 - t149 * t248 + t153 * t249;
t96 = t105 * t272;
t26 = t79 * t254 + t80 * t256 + t96;
t340 = t10 / 0.2e1 + t26 / 0.2e1 + t25 / 0.2e1;
t59 = t63 * t291;
t11 = t46 * t274 + t47 * t275 + t59;
t97 = t104 * t291;
t27 = t77 * t274 + t78 * t275 + t97;
t98 = t105 * t291;
t28 = t79 * t274 + t80 * t275 + t98;
t339 = t11 / 0.2e1 + t28 / 0.2e1 + t27 / 0.2e1;
t62 = t63 * t308;
t12 = t62 + (t47 * t313 - t46 * t316) * t307;
t102 = t104 * t308;
t29 = t102 + (t78 * t313 - t77 * t316) * t307;
t103 = t105 * t308;
t30 = t103 + (t80 * t313 - t79 * t316) * t307;
t338 = t12 / 0.2e1 + t30 / 0.2e1 + t29 / 0.2e1;
t335 = t307 * (-t278 + t369);
t330 = t307 * (-t278 + t357);
t328 = t172 * t385 + t173 * t383 + t347;
t327 = t307 * (-t278 + t348);
t202 = rSges(4,1) * t255 - rSges(4,2) * t254 + rSges(4,3) * t274;
t326 = -t259 + t355;
t267 = t293 * rSges(3,1) + t292 * rSges(3,2) - rSges(3,3) * t383;
t325 = t79 / 0.2e1 + t77 / 0.2e1 + t89 / 0.2e1 + t87 / 0.2e1 + t364;
t324 = t88 / 0.2e1 + t90 / 0.2e1 + t80 / 0.2e1 + t78 / 0.2e1 + t363;
t323 = t260 + t365;
t322 = -t212 + t326;
t196 = Icges(4,5) * t255 - Icges(4,6) * t254 + Icges(4,3) * t274;
t198 = Icges(4,4) * t255 - Icges(4,2) * t254 + Icges(4,6) * t274;
t200 = Icges(4,1) * t255 - Icges(4,4) * t254 + Icges(4,5) * t274;
t110 = t196 * t291 - t198 * t272 + t200 * t273;
t125 = t233 * t274 - t234 * t254 + t235 * t255;
t321 = t110 / 0.2e1 + t125 / 0.2e1 + t325;
t197 = Icges(4,5) * t257 - Icges(4,6) * t256 + Icges(4,3) * t275;
t199 = Icges(4,4) * t257 - Icges(4,2) * t256 + Icges(4,6) * t275;
t201 = Icges(4,1) * t257 - Icges(4,4) * t256 + Icges(4,5) * t275;
t111 = t197 * t291 - t199 * t272 + t201 * t273;
t126 = t233 * t275 - t234 * t256 + t235 * t257;
t320 = t126 / 0.2e1 + t111 / 0.2e1 + t324;
t319 = -t221 + t322;
t318 = t213 + t323;
t317 = t173 + t318;
t301 = rSges(2,1) * t316 - t313 * rSges(2,2);
t300 = -t313 * rSges(2,1) - rSges(2,2) * t316;
t285 = rSges(3,3) * t308 + (rSges(3,1) * t312 + rSges(3,2) * t315) * t307;
t266 = Icges(3,1) * t295 + Icges(3,4) * t294 + Icges(3,5) * t385;
t265 = Icges(3,1) * t293 + Icges(3,4) * t292 - Icges(3,5) * t383;
t264 = Icges(3,4) * t295 + Icges(3,2) * t294 + Icges(3,6) * t385;
t263 = Icges(3,4) * t293 + Icges(3,2) * t292 - Icges(3,6) * t383;
t262 = Icges(3,5) * t295 + Icges(3,6) * t294 + Icges(3,3) * t385;
t253 = t268 + t365;
t252 = -t267 + t355;
t239 = -t308 * t267 - t285 * t383;
t238 = t268 * t308 - t285 * t385;
t237 = t356 * t308;
t217 = (t267 * t313 + t268 * t316) * t307;
t216 = t282 * t385 + t283 * t294 + t284 * t295;
t215 = -t282 * t383 + t292 * t283 + t293 * t284;
t192 = t262 * t308 + (t264 * t315 + t266 * t312) * t307;
t191 = t261 * t308 + (t263 * t315 + t265 * t312) * t307;
t176 = t254 * t211;
t168 = t323 + t203;
t167 = -t202 + t326;
t160 = t272 * t173;
t158 = t256 * t172;
t137 = t203 * t291 - t236 * t275;
t136 = -t202 * t291 + t236 * t274;
t134 = (-t202 - t259) * t308 + t316 * t350;
t133 = t203 * t308 + t313 * t350 + t258;
t131 = t132 * t308;
t128 = t202 * t275 - t203 * t274;
t127 = t132 * t291;
t116 = (t202 * t313 + t203 * t316) * t307 + t366;
t115 = t318 + t157;
t114 = -t156 + t322;
t113 = t157 * t272 - t194 * t256;
t112 = -t156 * t272 + t194 * t254;
t109 = t197 * t275 - t199 * t256 + t201 * t257;
t108 = t196 * t275 - t198 * t256 + t200 * t257;
t107 = t197 * t274 - t199 * t254 + t201 * t255;
t106 = t196 * t274 - t198 * t254 + t200 * t255;
t101 = t156 * t256 - t157 * t254;
t100 = t317 + t155;
t99 = (rSges(6,2) - pkin(4)) * t229 + t319 + t332;
t94 = t157 * t291 + t275 * t369 + t205;
t93 = t194 * t274 + t291 * t373 + t218;
t92 = (-t259 + t373) * t308 + t316 * t335;
t91 = t157 * t308 + t313 * t335 + t367;
t86 = t124 * t249 - t141 * t231;
t85 = -t123 * t249 + t141 * t229;
t84 = t156 * t275 + t204 + (-t157 - t213) * t274;
t83 = (t156 * t313 + t157 * t316) * t307 + t347;
t82 = t155 * t272 + t256 * t370 + t160;
t81 = t193 * t254 + t272 * t375 + t176;
t76 = t317 + t377;
t75 = -t390 + (-rSges(7,3) - pkin(4) - pkin(12)) * t229 + t319 + t331;
t66 = (-t259 + t359) * t308 + t316 * t330;
t65 = t155 * t308 + t313 * t330 + t358;
t64 = t123 * t231 - t124 * t229;
t61 = t155 * t291 + t275 * t357 + t371;
t60 = t193 * t274 + t291 * t359 + t368;
t56 = t154 * t256 + t254 * t374 + t158;
t53 = (t154 * t313 + t155 * t316) * t307 + t328;
t52 = t154 * t275 + (-t213 + t374) * t274 + t372;
t51 = t131 + (-t110 * t316 + t111 * t313) * t307;
t50 = t110 * t274 + t111 * t275 + t127;
t49 = t256 * t360 + t272 * t377 + t160;
t48 = t254 * t376 + t272 * t362 + t176;
t45 = (-t259 + t349) * t308 + t316 * t327;
t44 = t308 * t377 + t313 * t327 + t358;
t43 = t126 * t308 + (-t108 * t316 + t109 * t313) * t307;
t42 = t125 * t308 + (-t106 * t316 + t107 * t313) * t307;
t41 = t275 * t348 + t291 * t377 + t371;
t40 = t274 * t376 + t291 * t349 + t368;
t35 = t108 * t274 + t109 * t275 + t126 * t291;
t34 = t106 * t274 + t107 * t275 + t125 * t291;
t33 = t254 * t361 + t256 * t378 + t158;
t32 = (t313 * t378 + t316 * t377) * t307 + t328;
t31 = t378 * t275 + (-t213 + t361) * t274 + t372;
t129 = [t356 + m(7) * (t75 ^ 2 + t76 ^ 2) + m(6) * (t100 ^ 2 + t99 ^ 2) + m(5) * (t114 ^ 2 + t115 ^ 2) + m(4) * (t167 ^ 2 + t168 ^ 2) + m(3) * (t252 ^ 2 + t253 ^ 2) + m(2) * (t300 ^ 2 + t301 ^ 2) + Icges(2,3) + t132 + t104 + t105 + t63; t103 + t237 + t102 + t131 + t62 + m(7) * (t44 * t76 + t45 * t75) + m(6) * (t100 * t65 + t66 * t99) + m(5) * (t114 * t92 + t115 * t91) + m(4) * (t133 * t168 + t134 * t167) + m(3) * (t238 * t253 + t239 * t252) + ((-t215 / 0.2e1 - t191 / 0.2e1 - t321) * t316 + (t192 / 0.2e1 + t216 / 0.2e1 + t320) * t313) * t307; (t12 + t30 + t29 + t51 + t237) * t308 + m(7) * (t32 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t53 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(5) * (t83 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(4) * (t116 ^ 2 + t133 ^ 2 + t134 ^ 2) + m(3) * (t217 ^ 2 + t238 ^ 2 + t239 ^ 2) + ((-t7 - t21 - t23 - t42 + (t292 * t263 + t293 * t265 - t307 * t387) * t383) * t316 + (t8 + t22 + t24 + t43 + ((t264 * t294 + t266 * t295 + (t262 * t313 - t387) * t307) * t313 + (t262 * t383 - t263 * t294 - t292 * t264 - t265 * t295 - t293 * t266) * t316) * t307) * t313 + ((-t191 - t215) * t316 + (t192 + t216) * t313) * t308) * t307; t59 + t97 + t98 + t127 + m(7) * (t40 * t75 + t41 * t76) + m(6) * (t100 * t61 + t60 * t99) + m(5) * (t114 * t93 + t115 * t94) + m(4) * (t136 * t167 + t137 * t168) + t320 * t275 + t321 * t274; (t50 / 0.2e1 + t339) * t308 + (t51 / 0.2e1 + t338) * t291 + (t43 / 0.2e1 + t341) * t275 + (t42 / 0.2e1 + t342) * t274 + m(7) * (t31 * t32 + t40 * t45 + t41 * t44) + m(6) * (t52 * t53 + t60 * t66 + t61 * t65) + m(5) * (t83 * t84 + t91 * t94 + t92 * t93) + m(4) * (t116 * t128 + t133 * t137 + t134 * t136) + ((-t34 / 0.2e1 - t344) * t316 + (t35 / 0.2e1 + t343) * t313) * t307; (t11 + t28 + t27 + t50) * t291 + (t6 + t20 + t18 + t35) * t275 + (t5 + t19 + t17 + t34) * t274 + m(7) * (t31 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(6) * (t52 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(5) * (t84 ^ 2 + t93 ^ 2 + t94 ^ 2) + m(4) * (t128 ^ 2 + t136 ^ 2 + t137 ^ 2); t58 + t95 + t96 + m(7) * (t48 * t75 + t49 * t76) + m(6) * (t100 * t82 + t81 * t99) + m(5) * (t112 * t114 + t113 * t115) + t324 * t256 + t325 * t254; t340 * t308 + t338 * t272 + t341 * t256 + t342 * t254 + m(7) * (t32 * t33 + t44 * t49 + t45 * t48) + m(6) * (t56 * t53 + t65 * t82 + t66 * t81) + m(5) * (t101 * t83 + t112 * t92 + t113 * t91) + (t313 * t345 - t316 * t346) * t307; t340 * t291 + t345 * t275 + t346 * t274 + t339 * t272 + t343 * t256 + t344 * t254 + m(7) * (t31 * t33 + t40 * t48 + t41 * t49) + m(6) * (t56 * t52 + t60 * t81 + t61 * t82) + m(5) * (t101 * t84 + t112 * t93 + t113 * t94); (t10 + t26 + t25) * t272 + (t4 + t14 + t16) * t256 + (t3 + t13 + t15) * t254 + m(7) * (t33 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t56 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(5) * (t101 ^ 2 + t112 ^ 2 + t113 ^ 2); m(7) * (t228 * t76 + t230 * t75) + m(6) * (t100 * t228 + t230 * t99); m(7) * (t228 * t44 + t230 * t45 + t248 * t32) + m(6) * (t228 * t65 + t230 * t66 + t248 * t53); m(7) * (t228 * t41 + t230 * t40 + t248 * t31) + m(6) * (t228 * t61 + t230 * t60 + t248 * t52); m(7) * (t228 * t49 + t230 * t48 + t248 * t33) + m(6) * (t228 * t82 + t230 * t81 + t248 * t56); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t228 ^ 2 + t230 ^ 2 + t248 ^ 2); m(7) * (t75 * t85 + t76 * t86) + t57 + t363 * t231 + t364 * t229; m(7) * (t64 * t32 + t44 * t86 + t45 * t85) + t8 * t394 + t7 * t395 + t308 * t396 + t12 * t393 + (t313 * t397 - t316 * t1 / 0.2e1) * t307; m(7) * (t64 * t31 + t40 * t85 + t41 * t86) + t291 * t396 + t6 * t394 + t5 * t395 + t275 * t397 + t11 * t393 + t274 * t398; t256 * t397 + t254 * t398 + t272 * t396 + t10 * t393 + m(7) * (t64 * t33 + t48 * t85 + t49 * t86) + t3 * t395 + t4 * t394; m(7) * (t228 * t86 + t230 * t85 + t248 * t64); t231 * t2 + t229 * t1 + t249 * t9 + m(7) * (t64 ^ 2 + t85 ^ 2 + t86 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t129(1) t129(2) t129(4) t129(7) t129(11) t129(16); t129(2) t129(3) t129(5) t129(8) t129(12) t129(17); t129(4) t129(5) t129(6) t129(9) t129(13) t129(18); t129(7) t129(8) t129(9) t129(10) t129(14) t129(19); t129(11) t129(12) t129(13) t129(14) t129(15) t129(20); t129(16) t129(17) t129(18) t129(19) t129(20) t129(21);];
Mq  = res;
