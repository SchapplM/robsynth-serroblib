% Calculate time derivative of joint inertia matrix for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR16_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR16_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR16_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:45
% EndTime: 2019-12-31 18:39:02
% DurationCPUTime: 11.00s
% Computational Cost: add. (5663->584), mult. (15280->838), div. (0->0), fcn. (14127->6), ass. (0->289)
t206 = cos(qJ(3));
t203 = sin(qJ(3));
t328 = Icges(5,6) * t203;
t336 = Icges(4,4) * t203;
t377 = -t328 - t336 + (Icges(4,1) + Icges(5,2)) * t206;
t327 = Icges(5,6) * t206;
t335 = Icges(4,4) * t206;
t376 = t327 + t335 + (-Icges(4,2) - Icges(5,3)) * t203;
t375 = t376 * qJD(3);
t374 = t377 * qJD(3);
t204 = sin(qJ(1));
t207 = cos(qJ(1));
t239 = Icges(4,2) * t206 + t336;
t119 = Icges(4,6) * t207 + t204 * t239;
t229 = Icges(5,3) * t206 + t328;
t125 = Icges(5,5) * t207 - t204 * t229;
t373 = -t119 + t125;
t242 = Icges(4,1) * t203 + t335;
t122 = Icges(4,5) * t207 + t204 * t242;
t231 = Icges(5,2) * t203 + t327;
t127 = Icges(5,4) * t207 - t204 * t231;
t372 = t122 - t127;
t126 = Icges(5,4) * t204 + t207 * t231;
t361 = -Icges(4,5) * t204 + t207 * t242;
t371 = t361 + t126;
t124 = Icges(5,5) * t204 + t207 * t229;
t362 = -Icges(4,6) * t204 + t207 * t239;
t370 = t362 + t124;
t302 = qJD(3) * t206;
t278 = t204 * t302;
t305 = qJD(1) * t207;
t211 = t203 * t305 + t278;
t205 = cos(qJ(5));
t202 = sin(qJ(5));
t333 = Icges(6,4) * t202;
t236 = Icges(6,2) * t205 + t333;
t118 = Icges(6,6) * t206 + t203 * t236;
t332 = Icges(6,4) * t205;
t241 = Icges(6,1) * t202 + t332;
t121 = Icges(6,5) * t206 + t203 * t241;
t369 = t118 * t205 + t121 * t202;
t368 = -qJD(3) / 0.2e1;
t270 = -qJD(1) * t206 - qJD(5);
t223 = t270 * t207;
t271 = qJD(5) * t206 + qJD(1);
t304 = qJD(3) * t203;
t71 = t205 * t223 + (t202 * t271 + t205 * t304) * t204;
t303 = qJD(3) * t204;
t280 = t203 * t303;
t319 = t204 * t205;
t72 = -t271 * t319 + (t223 + t280) * t202;
t40 = t72 * rSges(6,1) + t71 * rSges(6,2) + rSges(6,3) * t211;
t367 = -pkin(7) * t211 - t40;
t225 = t125 * t206 + t127 * t203;
t366 = t207 * t225;
t228 = t119 * t206 + t122 * t203;
t365 = t207 * t228;
t266 = rSges(4,1) * t203 + rSges(4,2) * t206;
t219 = t207 * t266;
t281 = t206 * t305;
t285 = rSges(4,1) * t211 + rSges(4,2) * t281;
t354 = -pkin(1) - pkin(6);
t296 = -rSges(4,3) + t354;
t311 = qJ(2) * t305 + qJD(2) * t204;
t63 = (-rSges(4,2) * t304 + qJD(1) * t296) * t204 + t285 + t311;
t346 = rSges(4,2) * t203;
t169 = rSges(4,1) * t206 - t346;
t192 = qJD(2) * t207;
t301 = qJD(3) * t207;
t64 = t192 + t169 * t301 + (t296 * t207 + (-qJ(2) - t266) * t204) * qJD(1);
t364 = t204 * t64 - t207 * t63;
t321 = t203 * t204;
t188 = pkin(3) * t321;
t318 = t204 * t206;
t289 = qJ(4) * t318;
t146 = t188 - t289;
t198 = t207 * pkin(4);
t322 = t202 * t207;
t144 = -t205 * t318 - t322;
t317 = t205 * t207;
t145 = -t202 * t318 + t317;
t87 = rSges(6,1) * t145 + rSges(6,2) * t144 + rSges(6,3) * t321;
t290 = -pkin(7) * t321 - t146 - t198 - t87;
t234 = Icges(4,5) * t203 + Icges(4,6) * t206;
t363 = -Icges(4,3) * t204 + t207 * t234;
t237 = Icges(5,4) * t203 + Icges(5,5) * t206;
t128 = Icges(5,1) * t204 + t207 * t237;
t264 = rSges(6,1) * t202 + rSges(6,2) * t205;
t299 = qJD(5) * t203;
t103 = (rSges(6,1) * t205 - rSges(6,2) * t202) * t299 + (-rSges(6,3) * t203 + t206 * t264) * qJD(3);
t130 = rSges(6,3) * t206 + t203 * t264;
t306 = qJD(1) * t204;
t277 = t206 * t301;
t283 = t203 * t306;
t210 = -t277 + t283;
t279 = t203 * t301;
t209 = t204 * t270 - t279;
t73 = t205 * t209 - t271 * t322;
t74 = t202 * t209 + t271 * t317;
t268 = t74 * rSges(6,1) + t73 * rSges(6,2);
t41 = rSges(6,3) * t210 + t268;
t316 = t206 * t207;
t142 = -t202 * t204 + t205 * t316;
t143 = t202 * t316 + t319;
t265 = -t143 * rSges(6,1) - t142 * rSges(6,2);
t320 = t203 * t207;
t86 = -rSges(6,3) * t320 - t265;
t21 = (-t130 * t301 - t41) * t206 + (qJD(3) * t86 - t103 * t207 + t130 * t306) * t203;
t220 = t103 * t204 + t130 * t305;
t22 = (-t130 * t303 + t40) * t206 + (-qJD(3) * t87 - t220) * t203;
t59 = -t130 * t321 + t206 * t87;
t60 = -t130 * t320 - t206 * t86;
t360 = qJD(1) * (t204 * t59 + t207 * t60) + t21 * t204 - t207 * t22;
t359 = 2 * m(4);
t358 = 2 * m(5);
t357 = 2 * m(6);
t200 = t204 ^ 2;
t201 = t207 ^ 2;
t356 = m(5) / 0.2e1;
t355 = m(6) / 0.2e1;
t353 = t204 / 0.2e1;
t352 = t206 / 0.2e1;
t351 = t207 / 0.2e1;
t350 = rSges(3,2) - pkin(1);
t349 = rSges(6,3) + pkin(7);
t348 = m(4) * t169;
t347 = t204 * pkin(4);
t345 = rSges(5,2) * t203;
t344 = rSges(5,2) * t206;
t343 = rSges(5,3) * t206;
t342 = t204 * rSges(5,1);
t341 = t204 * rSges(4,3);
t196 = t207 * rSges(5,1);
t195 = t207 * rSges(4,3);
t295 = pkin(7) * t320;
t338 = -t295 + t347 + t86;
t325 = qJ(4) * t206;
t138 = qJD(4) * t203 + (-pkin(3) * t203 + t325) * qJD(3);
t167 = pkin(3) * t206 + qJ(4) * t203;
t315 = t138 * t204 + t167 * t305;
t262 = t343 + t345;
t314 = t204 * t262 - t146 - t196;
t189 = pkin(3) * t320;
t194 = t207 * qJ(2);
t312 = t189 + t194;
t310 = pkin(1) * t207 + qJ(2) * t204;
t309 = t200 + t201;
t116 = Icges(4,3) * t207 + t204 * t234;
t308 = qJD(1) * t116;
t129 = Icges(5,1) * t207 - t204 * t237;
t307 = qJD(1) * t129;
t300 = qJD(4) * t206;
t298 = -pkin(4) + t354;
t297 = -rSges(5,1) + t354;
t233 = Icges(6,5) * t202 + Icges(6,6) * t205;
t115 = Icges(6,3) * t206 + t203 * t233;
t88 = (Icges(6,5) * t205 - Icges(6,6) * t202) * t299 + (-Icges(6,3) * t203 + t206 * t233) * qJD(3);
t91 = (-Icges(6,2) * t202 + t332) * t299 + (-Icges(6,6) * t203 + t206 * t236) * qJD(3);
t94 = (Icges(6,1) * t205 - t333) * t299 + (-Icges(6,5) * t203 + t206 * t241) * qJD(3);
t16 = t115 * t210 + t118 * t73 + t121 * t74 + t142 * t91 + t143 * t94 - t320 * t88;
t81 = Icges(6,4) * t143 + Icges(6,2) * t142 - Icges(6,6) * t320;
t83 = Icges(6,1) * t143 + Icges(6,4) * t142 - Icges(6,5) * t320;
t261 = t202 * t83 + t205 * t81;
t35 = Icges(6,5) * t74 + Icges(6,6) * t73 + Icges(6,3) * t210;
t37 = Icges(6,4) * t74 + Icges(6,2) * t73 + Icges(6,6) * t210;
t39 = Icges(6,1) * t74 + Icges(6,4) * t73 + Icges(6,5) * t210;
t79 = Icges(6,5) * t143 + Icges(6,6) * t142 - Icges(6,3) * t320;
t9 = (qJD(3) * t261 + t35) * t206 + (-qJD(3) * t79 + t202 * t39 + t205 * t37 + (-t202 * t81 + t205 * t83) * qJD(5)) * t203;
t294 = -t16 / 0.2e1 - t9 / 0.2e1;
t82 = Icges(6,4) * t145 + Icges(6,2) * t144 + Icges(6,6) * t321;
t84 = Icges(6,1) * t145 + Icges(6,4) * t144 + Icges(6,5) * t321;
t260 = t202 * t84 + t205 * t82;
t34 = Icges(6,5) * t72 + Icges(6,6) * t71 + Icges(6,3) * t211;
t36 = Icges(6,4) * t72 + Icges(6,2) * t71 + Icges(6,6) * t211;
t38 = Icges(6,1) * t72 + Icges(6,4) * t71 + Icges(6,5) * t211;
t80 = Icges(6,5) * t145 + Icges(6,6) * t144 + Icges(6,3) * t321;
t10 = (qJD(3) * t260 + t34) * t206 + (-qJD(3) * t80 + t202 * t38 + t205 * t36 + (-t202 * t82 + t205 * t84) * qJD(5)) * t203;
t15 = t115 * t211 + t71 * t118 + t72 * t121 + t144 * t91 + t145 * t94 + t321 * t88;
t293 = t15 / 0.2e1 + t10 / 0.2e1;
t30 = t203 * t261 + t206 * t79;
t44 = -t115 * t320 + t118 * t142 + t121 * t143;
t292 = t44 / 0.2e1 + t30 / 0.2e1;
t31 = t203 * t260 + t206 * t80;
t45 = t115 * t321 + t118 * t144 + t121 * t145;
t291 = t45 / 0.2e1 + t31 / 0.2e1;
t288 = qJ(4) * t316;
t287 = pkin(3) * t211 + qJ(4) * t280;
t286 = pkin(3) * t277 + qJ(4) * t279 + qJD(1) * t289;
t131 = rSges(4,1) * t321 + rSges(4,2) * t318 + t195;
t284 = pkin(6) * t207 + t310;
t276 = qJD(5) * t202 * t118;
t275 = t237 * qJD(3) / 0.2e1 + t234 * t368;
t274 = pkin(7) * t206 + t130;
t158 = t266 * qJD(3);
t273 = t158 * t309;
t272 = t298 * t204;
t269 = t192 + t286;
t267 = t274 * t204;
t263 = -rSges(5,3) * t203 + t344;
t245 = t287 + t311;
t23 = -t204 * t300 + (t272 - t288) * qJD(1) + t245 - t367;
t24 = (qJD(3) * t349 - qJD(4)) * t316 + (t298 * t207 + (-qJ(2) + (-pkin(3) - t349) * t203) * t204) * qJD(1) - t268 + t269;
t255 = t204 * t24 - t207 * t23;
t25 = t142 * t81 + t143 * t83 - t320 * t79;
t26 = t142 * t82 + t143 * t84 - t320 * t80;
t254 = t204 * t26 - t207 * t25;
t17 = t204 * t25 + t207 * t26;
t27 = t144 * t81 + t145 * t83 + t321 * t79;
t28 = t144 * t82 + t145 * t84 + t321 * t80;
t253 = t204 * t28 - t207 * t27;
t18 = t204 * t27 + t207 * t28;
t252 = t204 * t31 - t207 * t30;
t251 = t204 * t30 + t207 * t31;
t148 = t167 * t306;
t32 = t148 + qJD(1) * t267 + (pkin(7) * t304 - t103 - t138) * t207;
t33 = (-t280 + t281) * pkin(7) + t220 + t315;
t250 = t204 * t33 - t207 * t32;
t173 = rSges(5,3) * t280;
t218 = -t345 + (-rSges(5,3) - qJ(4)) * t206;
t208 = t204 * t297 + t207 * t218;
t42 = t173 + (-rSges(5,2) * qJD(3) - qJD(4)) * t318 + t208 * qJD(1) + t245;
t43 = (-qJD(3) * t263 - t300) * t207 + (t297 * t207 + (t343 - qJ(2) + (rSges(5,2) - pkin(3)) * t203) * t204) * qJD(1) + t269;
t249 = t204 * t43 - t207 * t42;
t157 = t262 * qJD(3);
t61 = -t263 * t306 + t148 + (-t138 - t157) * t207;
t62 = t157 * t204 - t263 * t305 + t315;
t247 = t204 * t62 - t207 * t61;
t246 = t204 * t86 + t207 * t87;
t238 = Icges(5,4) * t206 - Icges(5,5) * t203;
t235 = Icges(4,5) * t206 - Icges(4,6) * t203;
t227 = -t203 * t361 - t206 * t362;
t226 = t124 * t206 + t126 * t203;
t224 = (t356 + t355) * t304;
t222 = t205 * t121 * t299 + t206 * t88 + t369 * t302 + (t202 * t94 + t205 * t91) * t203;
t221 = rSges(3,3) * t207 + t204 * t350;
t217 = t227 * t204;
t216 = t226 * t204;
t150 = t204 * t167;
t147 = -t189 + t288;
t137 = t207 * t147;
t136 = -rSges(3,2) * t207 + rSges(3,3) * t204 + t310;
t135 = t194 + t221;
t133 = t207 * t262 + t342;
t132 = t341 - t219;
t113 = (-t167 + t263) * t207;
t112 = -t204 * t263 + t150;
t107 = t192 + (t350 * t207 + (-rSges(3,3) - qJ(2)) * t204) * qJD(1);
t106 = qJD(1) * t221 + t311;
t105 = t284 + t131;
t104 = t204 * t296 + t194 + t219;
t102 = t238 * t301 + t307;
t101 = -qJD(1) * t128 - t238 * t303;
t90 = qJD(1) * t363 + t235 * t303;
t89 = -t235 * t301 + t308;
t85 = (-qJ(4) * t305 - qJD(4) * t204) * t206 + t287;
t78 = (-t167 - t274) * t207;
t77 = t150 + t267;
t75 = t207 * (pkin(3) * t283 + t207 * t300 - t286);
t70 = t204 * t218 + t188 + t196 + t284;
t69 = t208 + t312;
t58 = t129 * t207 - t204 * t225;
t57 = t128 * t207 - t216;
t56 = t129 * t204 + t366;
t55 = t204 * t128 + t207 * t226;
t54 = -t204 * t363 - t207 * t227;
t53 = t116 * t204 - t365;
t52 = -t207 * t363 + t217;
t51 = t116 * t207 + t204 * t228;
t50 = t115 * t206 + t203 * t369;
t49 = t284 - t290;
t48 = (t203 * t349 - t325) * t207 + t272 + t265 + t312;
t47 = t133 * t207 + t204 * t314 + t137;
t46 = t246 * t203;
t29 = t204 * t290 + t207 * t338 + t137;
t20 = t75 + (-t85 - t173) * t204 + (t200 * t344 + t201 * t263) * qJD(3) + ((t314 + t196) * t207 + (-t133 - t147 + t342) * t204) * qJD(1);
t19 = ((-qJD(3) * t115 - t276) * t203 + t222) * t206;
t14 = t246 * t302 + (t204 * t41 + t207 * t40 + (-t204 * t87 + t207 * t86) * qJD(1)) * t203;
t13 = t203 * t253 + t206 * t45;
t12 = t203 * t254 + t206 * t44;
t11 = t75 + (-pkin(7) * t277 + t41) * t207 + (-t85 + t367) * t204 + ((t290 + t198) * t207 + (-t147 + t295 - t338 + t347) * t204) * qJD(1);
t8 = -t80 * t277 + t142 * t36 + t143 * t38 + t73 * t82 + t74 * t84 + (-t207 * t34 + t306 * t80) * t203;
t7 = -t79 * t277 + t142 * t37 + t143 * t39 + t73 * t81 + t74 * t83 + (-t207 * t35 + t306 * t79) * t203;
t6 = t80 * t278 + t144 * t36 + t145 * t38 + t71 * t82 + t72 * t84 + (t204 * t34 + t305 * t80) * t203;
t5 = t79 * t278 + t144 * t37 + t145 * t39 + t71 * t81 + t72 * t83 + (t204 * t35 + t305 * t79) * t203;
t4 = -qJD(1) * t254 + t204 * t7 + t207 * t8;
t3 = -qJD(1) * t253 + t204 * t5 + t207 * t6;
t2 = (qJD(3) * t254 + t16) * t206 + (qJD(1) * t17 - qJD(3) * t44 + t204 * t8 - t207 * t7) * t203;
t1 = (qJD(3) * t253 + t15) * t206 + (qJD(1) * t18 - qJD(3) * t45 + t204 * t6 - t207 * t5) * t203;
t65 = [(t23 * t49 + t24 * t48) * t357 + (t42 * t70 + t43 * t69) * t358 + (t104 * t64 + t105 * t63) * t359 + 0.2e1 * m(3) * (t106 * t136 + t107 * t135) + t222 + (-t242 - t231) * t302 - t375 * t206 + (t239 + t229 - t115) * t304 + (-t276 - t374) * t203; m(6) * ((t204 * t49 + t207 * t48) * qJD(1) + t255) + m(5) * ((t204 * t70 + t207 * t69) * qJD(1) + t249) + m(4) * ((t104 * t207 + t105 * t204) * qJD(1) + t364) + m(3) * (-t106 * t207 + t204 * t107 + (t135 * t207 + t136 * t204) * qJD(1)); 0; (t207 * t275 + t293) * t207 + (t204 * t275 - t294) * t204 + m(4) * (t364 * t169 - (t104 * t204 - t105 * t207) * t158) + m(6) * (t23 * t78 + t24 * t77 + t32 * t49 + t33 * t48) + m(5) * (t112 * t43 + t113 * t42 + t61 * t70 + t62 * t69) + ((t370 * qJD(3) - t301 * t377) * t353 + (t373 * qJD(3) + t204 * t374) * t351 + (t351 * t371 + t353 * t372) * qJD(1)) * t206 + ((t371 * qJD(3) + t301 * t376) * t353 + (-t372 * qJD(3) - t204 * t375) * t351 + (-t351 * t370 + t353 * t373) * qJD(1)) * t203 + ((t105 * t348 + (-t122 / 0.2e1 + t127 / 0.2e1) * t206 + (t119 / 0.2e1 - t125 / 0.2e1) * t203 - t291) * t204 + (t104 * t348 + (-t126 / 0.2e1 - t361 / 0.2e1) * t206 + (t124 / 0.2e1 + t362 / 0.2e1) * t203 + t292) * t207) * qJD(1); m(5) * ((t112 * t207 + t113 * t204) * qJD(1) + t247) + m(6) * ((t204 * t78 + t207 * t77) * qJD(1) + t250) - m(4) * t273; (t11 * t29 + t32 * t78 + t33 * t77) * t357 + t204 * t4 + t207 * t3 + (t112 * t62 + t113 * t61 + t20 * t47) * t358 + t207 * ((t101 * t207 + (t57 - t366) * qJD(1)) * t207 + (-t58 * qJD(1) + (t124 * t304 - t126 * t302) * t204 + (t102 + (t125 * t203 - t127 * t206) * qJD(3) + (-t129 - t226) * qJD(1)) * t207) * t204) + t204 * ((t204 * t102 + (-t56 - t216) * qJD(1)) * t204 + (t55 * qJD(1) + (-t125 * t304 + t127 * t302 + t307) * t207 + (t101 + (-t124 * t203 + t126 * t206) * qJD(3) - t225 * qJD(1)) * t204) * t207) + t204 * ((t204 * t89 + (-t53 + t217) * qJD(1)) * t204 + (t54 * qJD(1) + (t119 * t304 - t122 * t302 + t308) * t207 + (t90 + (-t203 * t362 + t206 * t361) * qJD(3) + t228 * qJD(1)) * t204) * t207) + t207 * ((t207 * t90 + (t52 + t365) * qJD(1)) * t207 + (-t51 * qJD(1) + (-t302 * t361 + t304 * t362) * t204 + (t89 + (-t119 * t203 + t122 * t206) * qJD(3) + (-t116 + t227) * qJD(1)) * t207) * t204) + ((-t131 * t204 + t132 * t207) * (-t204 * t285 + (-t169 * t201 + t200 * t346) * qJD(3) + ((-t131 + t195) * t207 + (-t132 + t219 + t341) * t204) * qJD(1)) - t169 * t273) * t359 + (-t18 + (-t51 - t58) * t207 + (-t52 - t57) * t204) * t306 + (t17 + (t53 + t56) * t207 + (t54 + t55) * t204) * t305; 0.2e1 * ((t204 * t48 - t207 * t49) * t355 + (t204 * t69 - t207 * t70) * t356) * t304 + 0.2e1 * ((-t305 * t48 - t306 * t49 - t255) * t355 + (-t305 * t69 - t306 * t70 - t249) * t356) * t206; 0.2e1 * t309 * t224; 0.2e1 * ((-t301 * t78 + t303 * t77 + t11) * t355 + (t112 * t303 - t113 * t301 + t20) * t356) * t203 + 0.2e1 * ((qJD(3) * t29 - t305 * t77 - t306 * t78 - t250) * t355 + (qJD(3) * t47 - t112 * t305 - t113 * t306 - t247) * t356) * t206; 0.4e1 * (0.1e1 - t309) * t206 * t224; m(6) * (t21 * t48 + t22 * t49 + t23 * t59 + t24 * t60) + t19 + (t204 * t291 - t207 * t292) * t302 + (-qJD(3) * t50 + t294 * t207 + t293 * t204 + (t204 * t292 + t207 * t291) * qJD(1)) * t203; m(6) * t360; m(6) * (t11 * t46 + t14 * t29 + t21 * t77 + t22 * t78 + t32 * t59 + t33 * t60) + (-t17 * t302 / 0.2e1 + t1 / 0.2e1 + qJD(1) * t12 / 0.2e1 + (qJD(1) * t30 + t10) * t352) * t207 + (t18 * t302 / 0.2e1 - qJD(1) * t13 / 0.2e1 + t2 / 0.2e1 + (-qJD(1) * t31 + t9) * t352) * t204 + (-t207 * t4 / 0.2e1 + t3 * t353 + t251 * t368 + (t17 * t353 + t18 * t351) * qJD(1)) * t203; m(6) * ((t14 + (t204 * t60 - t207 * t59) * qJD(3)) * t203 + (qJD(3) * t46 - t360) * t206); (t14 * t46 + t21 * t60 + t22 * t59) * t357 + (t19 + (-t12 * t207 + t13 * t204 + t206 * t252) * qJD(3)) * t206 + (t204 * t1 - t207 * t2 + t206 * (t10 * t204 - t207 * t9) + (-t203 * t252 - 0.2e1 * t206 * t50) * qJD(3) + (t12 * t204 + t13 * t207 + t206 * t251) * qJD(1)) * t203;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t65(1), t65(2), t65(4), t65(7), t65(11); t65(2), t65(3), t65(5), t65(8), t65(12); t65(4), t65(5), t65(6), t65(9), t65(13); t65(7), t65(8), t65(9), t65(10), t65(14); t65(11), t65(12), t65(13), t65(14), t65(15);];
Mq = res;
