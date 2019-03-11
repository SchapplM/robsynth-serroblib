% Calculate time derivative of joint inertia matrix for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR6_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR6_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR6_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:50:43
% EndTime: 2019-03-09 01:50:59
% DurationCPUTime: 11.10s
% Computational Cost: add. (6112->633), mult. (16270->902), div. (0->0), fcn. (14873->6), ass. (0->305)
t204 = sin(qJ(4));
t207 = cos(qJ(4));
t334 = Icges(6,6) * t207;
t342 = Icges(5,4) * t207;
t388 = (-t334 - t342 + (Icges(5,2) + Icges(6,3)) * t204) * qJD(4);
t335 = Icges(6,6) * t204;
t343 = Icges(5,4) * t204;
t387 = (-t335 - t343 + (Icges(5,1) + Icges(6,2)) * t207) * qJD(4);
t205 = sin(qJ(1));
t208 = cos(qJ(1));
t248 = Icges(5,2) * t207 + t343;
t124 = Icges(5,6) * t208 + t205 * t248;
t242 = Icges(6,3) * t207 + t335;
t375 = -Icges(6,5) * t208 + t205 * t242;
t386 = t124 + t375;
t125 = -Icges(5,6) * t205 + t208 * t248;
t129 = -Icges(6,5) * t205 - t208 * t242;
t385 = -t125 + t129;
t250 = Icges(5,1) * t204 + t342;
t127 = Icges(5,5) * t208 + t205 * t250;
t243 = Icges(6,2) * t204 + t334;
t374 = -Icges(6,4) * t208 + t205 * t243;
t384 = -t127 - t374;
t128 = -Icges(5,5) * t205 + t208 * t250;
t131 = -Icges(6,4) * t205 - t208 * t243;
t383 = t128 - t131;
t306 = qJD(4) * t208;
t286 = t204 * t306;
t311 = qJD(1) * t207;
t382 = t205 * t311 + t286;
t206 = cos(qJ(6));
t203 = sin(qJ(6));
t340 = Icges(7,4) * t203;
t246 = Icges(7,2) * t206 + t340;
t123 = Icges(7,6) * t207 + t204 * t246;
t339 = Icges(7,4) * t206;
t249 = Icges(7,1) * t203 + t339;
t126 = Icges(7,5) * t207 + t204 * t249;
t381 = t123 * t206 + t126 * t203;
t380 = -qJD(4) / 0.2e1;
t326 = t204 * t208;
t324 = t207 * t208;
t147 = t205 * t203 - t206 * t324;
t148 = -t203 * t324 - t205 * t206;
t86 = t148 * rSges(7,1) + t147 * rSges(7,2) + rSges(7,3) * t326;
t379 = -pkin(8) * t326 - t86;
t201 = t205 ^ 2;
t202 = t208 ^ 2;
t315 = t201 + t202;
t199 = t208 * qJ(2);
t274 = rSges(5,1) * t204 + rSges(5,2) * t207;
t353 = -pkin(1) - qJ(3);
t227 = -t274 + t353;
t359 = -rSges(5,3) - pkin(7);
t209 = t205 * t227 + t208 * t359;
t106 = t199 + t209;
t316 = t208 * pkin(1) + t205 * qJ(2);
t292 = t208 * qJ(3) + t316;
t319 = rSges(5,1) * t326 + rSges(5,2) * t324;
t107 = t205 * t359 + t292 + t319;
t378 = t106 * t208 + t107 * t205;
t284 = t207 * t306;
t180 = rSges(5,1) * t284;
t310 = qJD(1) * t208;
t318 = qJ(2) * t310 + qJD(2) * t205;
t293 = qJD(3) * t208 + t318;
t59 = -rSges(5,2) * t286 + qJD(1) * t209 + t180 + t293;
t352 = rSges(5,2) * t204;
t174 = rSges(5,1) * t207 - t352;
t196 = qJD(2) * t208;
t312 = qJD(1) * t205;
t317 = pkin(7) * t312 + t196;
t60 = (-qJD(4) * t174 - qJD(3)) * t205 + ((rSges(5,3) - qJ(2)) * t205 + t227 * t208) * qJD(1) + t317;
t377 = t205 * t59 + t208 * t60;
t245 = Icges(5,5) * t204 + Icges(5,6) * t207;
t121 = Icges(5,3) * t208 + t205 * t245;
t376 = qJD(1) * t121;
t247 = Icges(6,4) * t204 + Icges(6,5) * t207;
t373 = -Icges(6,1) * t208 + t205 * t247;
t272 = rSges(7,1) * t203 + rSges(7,2) * t206;
t304 = qJD(6) * t204;
t103 = (rSges(7,1) * t206 - rSges(7,2) * t203) * t304 + (-rSges(7,3) * t204 + t207 * t272) * qJD(4);
t136 = rSges(7,3) * t207 + t204 * t272;
t290 = t204 * t312;
t277 = qJD(6) + t311;
t211 = t205 * t277 + t286;
t278 = qJD(6) * t207 + qJD(1);
t234 = t278 * t203;
t73 = t206 * t211 + t208 * t234;
t235 = t278 * t206;
t74 = t203 * t211 - t208 * t235;
t300 = t74 * rSges(7,1) + t73 * rSges(7,2) + rSges(7,3) * t284;
t43 = -rSges(7,3) * t290 + t300;
t21 = (-t136 * t306 + t43) * t207 + (-qJD(4) * t86 - t103 * t208 + t136 * t312) * t204;
t229 = t205 * t103 + t136 * t310;
t308 = qJD(4) * t205;
t307 = qJD(4) * t207;
t285 = t205 * t307;
t214 = t204 * t310 + t285;
t233 = t277 * t208;
t309 = qJD(4) * t204;
t71 = -t206 * t233 + (t206 * t309 + t234) * t205;
t287 = t204 * t308;
t72 = -t205 * t235 + (-t233 + t287) * t203;
t275 = t72 * rSges(7,1) + t71 * rSges(7,2);
t42 = rSges(7,3) * t214 + t275;
t325 = t205 * t207;
t149 = -t203 * t208 - t206 * t325;
t150 = -t203 * t325 + t206 * t208;
t273 = -t150 * rSges(7,1) - t149 * rSges(7,2);
t327 = t204 * t205;
t87 = rSges(7,3) * t327 - t273;
t22 = (t136 * t308 - t42) * t207 + (qJD(4) * t87 + t229) * t204;
t61 = t136 * t327 - t207 * t87;
t62 = -t136 * t326 + t207 * t86;
t372 = -qJD(1) * (t205 * t61 - t208 * t62) + t21 * t205 + t208 * t22;
t371 = 2 * m(5);
t370 = 2 * m(6);
t369 = 2 * m(7);
t368 = m(6) / 0.2e1;
t367 = m(7) / 0.2e1;
t366 = -pkin(5) - pkin(7);
t365 = -t205 / 0.2e1;
t364 = t207 / 0.2e1;
t363 = t208 / 0.2e1;
t362 = -rSges(6,1) - pkin(7);
t361 = rSges(3,2) - pkin(1);
t360 = rSges(6,2) - pkin(4);
t358 = m(5) * t174;
t357 = pkin(4) * t204;
t356 = pkin(5) * t208;
t355 = t205 * pkin(5);
t354 = -qJD(1) / 0.2e1;
t351 = rSges(6,2) * t204;
t350 = rSges(6,3) * t207;
t348 = t208 * rSges(6,1);
t346 = -rSges(6,3) - qJ(5);
t345 = t355 + t379;
t332 = qJ(5) * t204;
t143 = qJD(5) * t204 + (qJ(5) * t207 - t357) * qJD(4);
t172 = pkin(4) * t207 + t332;
t323 = t205 * t143 + t172 * t310;
t186 = qJ(5) * t325;
t151 = pkin(4) * t327 - t186;
t271 = t350 + t351;
t228 = t205 * t271;
t322 = t228 - t348 - t151;
t288 = t207 * t310;
t305 = qJD(5) * t207;
t321 = -qJ(5) * t288 - t205 * t305;
t320 = t186 + t199;
t122 = -Icges(5,3) * t205 + t208 * t245;
t314 = qJD(1) * t122;
t133 = -Icges(6,1) * t205 - t208 * t247;
t313 = qJD(1) * t133;
t303 = -rSges(7,3) - pkin(4) - pkin(8);
t302 = -rSges(4,3) + t353;
t244 = Icges(7,5) * t203 + Icges(7,6) * t206;
t120 = Icges(7,3) * t207 + t204 * t244;
t213 = t284 - t290;
t88 = (Icges(7,5) * t206 - Icges(7,6) * t203) * t304 + (-Icges(7,3) * t204 + t207 * t244) * qJD(4);
t91 = (-Icges(7,2) * t203 + t339) * t304 + (-Icges(7,6) * t204 + t207 * t246) * qJD(4);
t94 = (Icges(7,1) * t206 - t340) * t304 + (-Icges(7,5) * t204 + t207 * t249) * qJD(4);
t16 = t120 * t213 + t73 * t123 + t74 * t126 + t147 * t91 + t148 * t94 + t326 * t88;
t80 = Icges(7,4) * t148 + Icges(7,2) * t147 + Icges(7,6) * t326;
t82 = Icges(7,1) * t148 + Icges(7,4) * t147 + Icges(7,5) * t326;
t270 = t203 * t82 + t206 * t80;
t37 = Icges(7,5) * t74 + Icges(7,6) * t73 + Icges(7,3) * t213;
t39 = Icges(7,4) * t74 + Icges(7,2) * t73 + Icges(7,6) * t213;
t41 = Icges(7,1) * t74 + Icges(7,4) * t73 + Icges(7,5) * t213;
t78 = Icges(7,5) * t148 + Icges(7,6) * t147 + Icges(7,3) * t326;
t9 = (qJD(4) * t270 + t37) * t207 + (-qJD(4) * t78 + t203 * t41 + t206 * t39 + (-t203 * t80 + t206 * t82) * qJD(6)) * t204;
t301 = t16 / 0.2e1 + t9 / 0.2e1;
t81 = Icges(7,4) * t150 + Icges(7,2) * t149 + Icges(7,6) * t327;
t83 = Icges(7,1) * t150 + Icges(7,4) * t149 + Icges(7,5) * t327;
t269 = t203 * t83 + t206 * t81;
t36 = Icges(7,5) * t72 + Icges(7,6) * t71 + Icges(7,3) * t214;
t38 = Icges(7,4) * t72 + Icges(7,2) * t71 + Icges(7,6) * t214;
t40 = Icges(7,1) * t72 + Icges(7,4) * t71 + Icges(7,5) * t214;
t79 = Icges(7,5) * t150 + Icges(7,6) * t149 + Icges(7,3) * t327;
t10 = (qJD(4) * t269 + t36) * t207 + (-qJD(4) * t79 + t203 * t40 + t206 * t38 + (-t203 * t81 + t206 * t83) * qJD(6)) * t204;
t15 = t120 * t214 + t71 * t123 + t72 * t126 + t149 * t91 + t150 * t94 + t327 * t88;
t299 = t15 / 0.2e1 + t10 / 0.2e1;
t30 = t204 * t270 + t207 * t78;
t44 = t120 * t326 + t147 * t123 + t148 * t126;
t298 = t44 / 0.2e1 + t30 / 0.2e1;
t31 = t204 * t269 + t207 * t79;
t45 = t120 * t327 + t123 * t149 + t126 * t150;
t297 = t45 / 0.2e1 + t31 / 0.2e1;
t296 = -pkin(8) * t327 - t151 - t356 - t87;
t295 = pkin(4) * t284 + qJ(5) * t382;
t294 = -rSges(6,2) * t290 - rSges(6,3) * t382;
t291 = t362 * t208;
t283 = qJD(6) * t123 * t203;
t282 = t247 * qJD(4) / 0.2e1 + t245 * t380;
t281 = t307 / 0.2e1;
t280 = pkin(8) * t207 + t136;
t163 = t274 * qJD(4);
t279 = t163 * t315;
t276 = t317 - t321;
t191 = pkin(4) * t326;
t152 = -qJ(5) * t324 + t191;
t183 = pkin(8) * t284;
t230 = t204 * t303 + t353;
t210 = t230 * t205 + t366 * t208;
t222 = -t208 * t305 + t295;
t23 = qJD(1) * t210 + t183 + t222 + t293 + t300;
t24 = (-qJD(3) + (t207 * t303 - t332) * qJD(4)) * t205 + ((pkin(5) - qJ(2)) * t205 + t230 * t208) * qJD(1) - t275 + t276;
t264 = t205 * t23 + t208 * t24;
t25 = t147 * t80 + t148 * t82 + t326 * t78;
t26 = t147 * t81 + t148 * t83 + t326 * t79;
t263 = t205 * t26 + t208 * t25;
t17 = -t25 * t205 + t208 * t26;
t27 = t149 * t80 + t150 * t82 + t327 * t78;
t28 = t149 * t81 + t150 * t83 + t327 * t79;
t262 = t205 * t28 + t208 * t27;
t18 = -t27 * t205 + t208 * t28;
t261 = t205 * t31 + t208 * t30;
t260 = -t30 * t205 + t31 * t208;
t119 = t208 * t143;
t32 = t119 + (-pkin(8) * t309 + t103) * t208 + (-t172 - t280) * t312;
t33 = (-t287 + t288) * pkin(8) + t229 + t323;
t259 = t33 * t205 + t208 * t32;
t34 = (-rSges(6,2) * qJD(4) - qJD(5)) * t324 + (t291 + (t353 - t357) * t205) * qJD(1) + t293 - t294 + t295;
t212 = t204 * t360 + t350 + t353;
t35 = (-qJD(3) + (t204 * t346 + t207 * t360) * qJD(4)) * t205 + ((rSges(6,1) - qJ(2)) * t205 + t212 * t208) * qJD(1) + t276;
t258 = t205 * t34 + t208 * t35;
t47 = t210 + t273 + t320;
t48 = t366 * t205 + t152 + t292 - t379;
t257 = t205 * t48 + t208 * t47;
t256 = t205 * t62 + t208 * t61;
t162 = t271 * qJD(4);
t173 = -rSges(6,2) * t207 + rSges(6,3) * t204;
t63 = t162 * t208 + t119 + (-t172 - t173) * t312;
t64 = t205 * t162 + t173 * t310 + t323;
t254 = t64 * t205 + t208 * t63;
t65 = t205 * t212 + t291 + t320;
t66 = t191 + t362 * t205 + (t207 * t346 - t351) * t208 + t292;
t253 = t205 * t66 + t208 * t65;
t252 = -t205 * t86 + t208 * t87;
t241 = t124 * t207 + t127 * t204;
t240 = t125 * t207 + t128 * t204;
t239 = -t129 * t207 - t131 * t204;
t238 = t204 * t374 + t207 * t375;
t236 = (t368 + t367) * t309;
t232 = t206 * t126 * t304 + t207 * t88 + t381 * t307 + (t203 * t94 + t206 * t91) * t204;
t231 = rSges(3,3) * t208 + t205 * t361;
t226 = t241 * t208;
t225 = t240 * t205;
t224 = t239 * t205;
t223 = t238 * t208;
t219 = qJD(4) * (-Icges(6,4) * t207 + Icges(6,5) * t204);
t218 = qJD(4) * (Icges(5,5) * t207 - Icges(5,6) * t204);
t215 = rSges(4,2) * t208 + t205 * t302;
t155 = t208 * t172;
t154 = t205 * t172;
t142 = -rSges(3,2) * t208 + t205 * rSges(3,3) + t316;
t141 = t199 + t231;
t139 = -t205 * rSges(6,1) - t208 * t271;
t138 = -rSges(5,3) * t205 + t319;
t137 = rSges(5,3) * t208 + t205 * t274;
t135 = t152 * t312;
t117 = t205 * rSges(4,2) + rSges(4,3) * t208 + t292;
t116 = t199 + t215;
t115 = t173 * t208 + t155;
t114 = t173 * t205 + t154;
t109 = t196 + (t361 * t208 + (-rSges(3,3) - qJ(2)) * t205) * qJD(1);
t108 = qJD(1) * t231 + t318;
t105 = -qJD(3) * t205 + t196 + ((-rSges(4,2) - qJ(2)) * t205 + t302 * t208) * qJD(1);
t104 = qJD(1) * t215 + t293;
t102 = t373 * qJD(1) + t208 * t219;
t101 = t205 * t219 + t313;
t90 = t205 * t218 + t314;
t89 = t208 * t218 - t376;
t85 = pkin(4) * t214 + qJ(5) * t287 + t321;
t84 = -pkin(4) * t290 + t222;
t77 = t208 * t280 + t155;
t76 = t205 * t280 + t154;
t58 = t205 * t238 - t208 * t373;
t57 = t133 * t208 + t224;
t56 = t205 * t373 + t223;
t55 = -t205 * t133 + t208 * t239;
t54 = -t205 * t122 + t208 * t240;
t53 = -t205 * t121 + t226;
t52 = t122 * t208 + t225;
t51 = t121 * t208 + t205 * t241;
t50 = t120 * t207 + t204 * t381;
t49 = (-t139 - t152) * t208 + t322 * t205;
t46 = t252 * t204;
t29 = (-t152 + t345) * t208 + t296 * t205;
t20 = t135 + (qJD(1) * t139 - t85 + (rSges(6,1) * qJD(1) + rSges(6,2) * t307 - rSges(6,3) * t309) * t205) * t205 + (rSges(6,2) * t284 - t84 + (t228 + t322 + t348) * qJD(1) + t294) * t208;
t19 = ((-qJD(4) * t120 - t283) * t204 + t232) * t207;
t14 = t252 * t307 + (-t205 * t43 + t208 * t42 + (-t205 * t87 - t208 * t86) * qJD(1)) * t204;
t13 = t204 * t262 + t45 * t207;
t12 = t204 * t263 + t44 * t207;
t11 = t135 + (-t183 - t43 - t84 + (t296 + t356) * qJD(1)) * t208 + (-pkin(8) * t285 - t42 - t85 + (-t345 + t355) * qJD(1)) * t205;
t8 = t79 * t284 + t147 * t38 + t148 * t40 + t73 * t81 + t74 * t83 + (t208 * t36 - t312 * t79) * t204;
t7 = t78 * t284 + t147 * t39 + t148 * t41 + t73 * t80 + t74 * t82 + (t208 * t37 - t312 * t78) * t204;
t6 = t79 * t285 + t149 * t38 + t150 * t40 + t71 * t81 + t72 * t83 + (t205 * t36 + t310 * t79) * t204;
t5 = t78 * t285 + t149 * t39 + t150 * t41 + t71 * t80 + t72 * t82 + (t205 * t37 + t310 * t78) * t204;
t4 = -qJD(1) * t263 - t7 * t205 + t208 * t8;
t3 = -qJD(1) * t262 - t5 * t205 + t208 * t6;
t2 = (qJD(4) * t263 + t16) * t207 + (qJD(1) * t17 - qJD(4) * t44 + t205 * t8 + t208 * t7) * t204;
t1 = (qJD(4) * t262 + t15) * t207 + (qJD(1) * t18 - qJD(4) * t45 + t205 * t6 + t208 * t5) * t204;
t67 = [t232 + (t23 * t48 + t24 * t47) * t369 + (t34 * t66 + t35 * t65) * t370 + (t106 * t60 + t107 * t59) * t371 + 0.2e1 * m(4) * (t104 * t117 + t105 * t116) + 0.2e1 * m(3) * (t108 * t142 + t109 * t141) + (-t250 - t243) * t307 + t388 * t207 + (t248 + t242 - t120) * t309 + (-t283 - t387) * t204; m(7) * (qJD(1) * t257 + t205 * t24 - t208 * t23) + m(6) * (qJD(1) * t253 + t205 * t35 - t208 * t34) + m(5) * (t378 * qJD(1) + t205 * t60 - t208 * t59) + m(4) * (-t104 * t208 + t205 * t105 + (t116 * t208 + t117 * t205) * qJD(1)) + m(3) * (-t108 * t208 + t205 * t109 + (t141 * t208 + t142 * t205) * qJD(1)); 0; m(7) * ((-t205 * t47 + t208 * t48) * qJD(1) + t264) + m(6) * ((-t205 * t65 + t208 * t66) * qJD(1) + t258) + m(5) * ((-t106 * t205 + t107 * t208) * qJD(1) + t377) + m(4) * (t205 * t104 + t105 * t208 + (-t116 * t205 + t117 * t208) * qJD(1)); 0; 0; (t208 * t282 + t299) * t208 + (t205 * t282 - t301) * t205 + m(5) * (-t378 * t163 + t377 * t174) + m(6) * (t114 * t34 + t115 * t35 + t63 * t65 + t64 * t66) + m(7) * (t23 * t76 + t24 * t77 + t32 * t47 + t33 * t48) + ((t385 * qJD(4) + t387 * t208) * t365 + (-t386 * qJD(4) + t387 * t205) * t363 + (t383 * t363 + t384 * t365) * qJD(1)) * t207 + ((-t383 * qJD(4) + t388 * t208) * t365 + (t384 * qJD(4) + t388 * t205) * t363 + (t385 * t363 + t386 * t365) * qJD(1)) * t204 + ((t107 * t358 + (t131 / 0.2e1 - t128 / 0.2e1) * t207 + (-t129 / 0.2e1 + t125 / 0.2e1) * t204 - t298) * t208 + (-t106 * t358 + (-t374 / 0.2e1 - t127 / 0.2e1) * t207 + (t375 / 0.2e1 + t124 / 0.2e1) * t204 - t297) * t205) * qJD(1); m(6) * (t63 * t205 - t208 * t64 + (t114 * t205 + t115 * t208) * qJD(1)) + m(7) * (t32 * t205 - t208 * t33 + (t205 * t76 + t208 * t77) * qJD(1)); m(6) * ((t114 * t208 - t115 * t205) * qJD(1) + t254) + m(7) * ((-t205 * t77 + t208 * t76) * qJD(1) + t259) - m(5) * t279; t208 * t3 - t205 * t4 + (t11 * t29 + t32 * t77 + t33 * t76) * t369 + (t114 * t64 + t115 * t63 + t20 * t49) * t370 - t205 * ((t205 * t89 + (-t53 + t225) * qJD(1)) * t205 + (-t54 * qJD(1) + (-t124 * t309 + t127 * t307 - t376) * t208 + (-t90 + (t125 * t204 - t128 * t207) * qJD(4) + (t122 - t241) * qJD(1)) * t205) * t208) + (-t174 * t279 + (-t208 * t180 + (-t174 * t201 + t202 * t352) * qJD(4) + (t315 * rSges(5,3) - t208 * t137 + t205 * t138) * qJD(1)) * (-t205 * t137 - t138 * t208)) * t371 - t205 * ((t205 * t102 + (-t56 + t224) * qJD(1)) * t205 + (-t55 * qJD(1) + (t307 * t374 - t309 * t375) * t208 + (-t101 + (-t129 * t204 + t131 * t207) * qJD(4) + (t133 - t238) * qJD(1)) * t205) * t208) + t208 * ((t208 * t90 + (-t52 + t226) * qJD(1)) * t208 + (-t51 * qJD(1) + (t125 * t309 - t128 * t307 + t314) * t205 + (-t89 + (-t124 * t204 + t127 * t207) * qJD(4) + (-t121 - t240) * qJD(1)) * t208) * t205) + t208 * ((t208 * t101 + (-t57 + t223) * qJD(1)) * t208 + (-t58 * qJD(1) + (-t129 * t309 + t131 * t307 + t313) * t205 + (-t102 + (-t204 * t375 + t207 * t374) * qJD(4) - t239 * qJD(1)) * t208) * t205) + (-t18 + (-t51 - t58) * t208 + (t52 + t57) * t205) * t312 + (-t17 + (-t53 - t56) * t208 + (t54 + t55) * t205) * t310; 0.2e1 * (t253 * t368 + t257 * t367) * t309 + 0.2e1 * ((-t310 * t48 + t312 * t47 - t264) * t367 + (-t310 * t66 + t312 * t65 - t258) * t368) * t207; 0; 0.2e1 * t315 * t236; 0.2e1 * ((t306 * t77 + t308 * t76 + t11) * t367 + (t114 * t308 + t115 * t306 + t20) * t368) * t204 + 0.2e1 * ((qJD(4) * t29 - t310 * t76 + t312 * t77 - t259) * t367 + (qJD(4) * t49 - t114 * t310 + t115 * t312 - t254) * t368) * t207; 0.4e1 * (0.1e1 - t315) * t207 * t236; m(7) * (t21 * t48 + t22 * t47 + t23 * t62 + t24 * t61) + t19 + (t205 * t297 + t208 * t298) * t307 + (-qJD(4) * t50 + t301 * t208 + t299 * t205 + (-t205 * t298 + t208 * t297) * qJD(1)) * t204; m(7) * (qJD(1) * t256 + t22 * t205 - t208 * t21); m(7) * t372; m(7) * (t11 * t46 + t14 * t29 + t21 * t76 + t22 * t77 + t32 * t61 + t33 * t62) + (t1 / 0.2e1 + t17 * t281 + (-qJD(1) * t30 + t10) * t364 + t12 * t354) * t208 + (t13 * t354 + t18 * t281 + (-qJD(1) * t31 - t9) * t364 - t2 / 0.2e1) * t205 + (t205 * t3 / 0.2e1 + t4 * t363 + t260 * t380 + (t17 * t365 + t18 * t363) * qJD(1)) * t204; m(7) * ((qJD(4) * t256 + t14) * t204 + (qJD(4) * t46 - t372) * t207); (t14 * t46 + t21 * t62 + t22 * t61) * t369 + (t19 + (t208 * t12 + t205 * t13 + t207 * t261) * qJD(4)) * t207 + (t208 * t2 + t205 * t1 + t207 * (t10 * t205 + t208 * t9) + (-t204 * t261 - 0.2e1 * t50 * t207) * qJD(4) + (-t205 * t12 + t208 * t13 + t207 * t260) * qJD(1)) * t204;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t67(1) t67(2) t67(4) t67(7) t67(11) t67(16); t67(2) t67(3) t67(5) t67(8) t67(12) t67(17); t67(4) t67(5) t67(6) t67(9) t67(13) t67(18); t67(7) t67(8) t67(9) t67(10) t67(14) t67(19); t67(11) t67(12) t67(13) t67(14) t67(15) t67(20); t67(16) t67(17) t67(18) t67(19) t67(20) t67(21);];
Mq  = res;
