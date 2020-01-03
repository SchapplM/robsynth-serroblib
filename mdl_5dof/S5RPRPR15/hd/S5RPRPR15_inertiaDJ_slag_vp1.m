% Calculate time derivative of joint inertia matrix for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR15_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR15_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR15_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR15_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:39
% EndTime: 2019-12-31 18:36:57
% DurationCPUTime: 9.42s
% Computational Cost: add. (10027->649), mult. (16650->939), div. (0->0), fcn. (15627->8), ass. (0->306)
t392 = Icges(5,5) / 0.2e1;
t391 = Icges(5,6) / 0.2e1;
t390 = Icges(5,3) / 0.2e1;
t222 = sin(qJ(1));
t364 = -t222 / 0.2e1;
t389 = -qJD(1) / 0.2e1;
t388 = rSges(5,3) + qJ(4);
t221 = sin(qJ(3));
t223 = cos(qJ(3));
t310 = qJD(3) * t223;
t285 = t222 * t310;
t224 = cos(qJ(1));
t313 = qJD(1) * t224;
t387 = t221 * t313 + t285;
t309 = qJD(3) * t224;
t286 = t221 * t309;
t314 = qJD(1) * t222;
t386 = t223 * t314 + t286;
t311 = qJD(3) * t222;
t287 = t221 * t311;
t288 = t223 * t313;
t228 = t287 - t288;
t215 = pkin(8) + qJ(5);
t206 = sin(t215);
t207 = cos(t215);
t344 = Icges(6,4) * t207;
t249 = -Icges(6,2) * t206 + t344;
t135 = Icges(6,6) * t221 + t223 * t249;
t345 = Icges(6,4) * t206;
t253 = Icges(6,1) * t207 - t345;
t136 = Icges(6,5) * t221 + t223 * t253;
t385 = -t135 * t206 + t136 * t207;
t218 = sin(pkin(8));
t219 = cos(pkin(8));
t254 = Icges(5,1) * t219 - Icges(5,4) * t218;
t143 = Icges(5,5) * t221 + t223 * t254;
t384 = -t143 / 0.2e1;
t333 = t218 * t224;
t199 = pkin(4) * t333;
t204 = pkin(4) * t219 + pkin(3);
t220 = -pkin(7) - qJ(4);
t325 = t222 * t223;
t331 = t221 * t222;
t329 = t222 * t206;
t146 = t207 * t224 - t221 * t329;
t328 = t222 * t207;
t147 = t206 * t224 + t221 * t328;
t90 = t147 * rSges(6,1) + t146 * rSges(6,2) - rSges(6,3) * t325;
t298 = -t204 * t331 - t220 * t325 - t199 - t90;
t278 = qJD(5) * t221 + qJD(1);
t238 = t278 * t222;
t277 = qJD(1) * t221 + qJD(5);
t373 = t224 * t277 + t285;
t81 = -t206 * t373 - t207 * t238;
t82 = -t206 * t238 + t207 * t373;
t47 = t82 * rSges(6,1) + t81 * rSges(6,2) + rSges(6,3) * t228;
t383 = t204 * t387 + t220 * t288 + t47;
t275 = rSges(4,1) * t221 + rSges(4,2) * t223;
t233 = t224 * t275;
t347 = Icges(4,4) * t221;
t251 = Icges(4,2) * t223 + t347;
t152 = Icges(4,6) * t224 + t222 * t251;
t346 = Icges(4,4) * t223;
t255 = Icges(4,1) * t221 + t346;
t154 = Icges(4,5) * t224 + t222 * t255;
t242 = t152 * t223 + t154 * t221;
t232 = t242 * t224;
t358 = pkin(3) - t204;
t382 = t358 * t221;
t327 = t222 * t218;
t164 = t219 * t224 - t221 * t327;
t326 = t222 * t219;
t165 = t221 * t326 + t333;
t202 = pkin(3) * t331;
t381 = -t165 * rSges(5,1) - t164 * rSges(5,2) - t202;
t380 = -t135 * t207 - t136 * t206;
t292 = rSges(4,1) * t387 + rSges(4,2) * t288;
t366 = -pkin(1) - pkin(6);
t306 = -rSges(4,3) + t366;
t312 = qJD(3) * t221;
t318 = qJ(2) * t313 + qJD(2) * t222;
t74 = (-rSges(4,2) * t312 + qJD(1) * t306) * t222 + t292 + t318;
t356 = rSges(4,2) * t221;
t186 = rSges(4,1) * t223 - t356;
t209 = qJD(2) * t224;
t75 = t209 + t186 * t309 + (t306 * t224 + (-qJ(2) - t275) * t222) * qJD(1);
t379 = t222 * t75 - t224 * t74;
t284 = t223 * t309;
t378 = t222 * t277 - t284;
t336 = t204 * t221;
t350 = rSges(6,3) - t220;
t377 = -t223 * t350 + t336;
t247 = Icges(4,5) * t221 + Icges(4,6) * t223;
t376 = -Icges(4,3) * t222 + t224 * t247;
t375 = -Icges(4,6) * t222 + t224 * t251;
t374 = -Icges(4,5) * t222 + t224 * t255;
t270 = rSges(6,1) * t207 - rSges(6,2) * t206;
t138 = rSges(6,3) * t221 + t223 * t270;
t239 = t224 * t278;
t79 = -t206 * t378 + t207 * t239;
t80 = t206 * t239 + t207 * t378;
t276 = t80 * rSges(6,1) + t79 * rSges(6,2);
t46 = -rSges(6,3) * t386 + t276;
t330 = t221 * t224;
t148 = t206 * t330 + t328;
t149 = -t207 * t330 + t329;
t271 = -t149 * rSges(6,1) - t148 * rSges(6,2);
t324 = t223 * t224;
t91 = rSges(6,3) * t324 - t271;
t307 = qJD(5) * t223;
t95 = (-rSges(6,1) * t206 - rSges(6,2) * t207) * t307 + (rSges(6,3) * t223 - t221 * t270) * qJD(3);
t21 = (-t138 * t309 - t46) * t221 + (-qJD(3) * t91 - t138 * t314 + t224 * t95) * t223;
t22 = (-t138 * t311 + t47) * t221 + (qJD(3) * t90 + t138 * t313 + t222 * t95) * t223;
t58 = t138 * t325 + t221 * t90;
t59 = t138 * t324 - t221 * t91;
t372 = qJD(1) * (t222 * t58 + t224 * t59) + t21 * t222 - t22 * t224;
t371 = 2 * m(4);
t370 = 2 * m(5);
t369 = 2 * m(6);
t216 = t222 ^ 2;
t217 = t224 ^ 2;
t368 = m(5) / 0.2e1;
t367 = m(6) / 0.2e1;
t365 = t221 / 0.2e1;
t362 = t224 / 0.2e1;
t361 = rSges(3,2) - pkin(1);
t360 = m(4) * t186;
t359 = pkin(3) * t221;
t245 = Icges(6,5) * t207 - Icges(6,6) * t206;
t134 = Icges(6,3) * t221 + t223 * t245;
t92 = (-Icges(6,5) * t206 - Icges(6,6) * t207) * t307 + (Icges(6,3) * t223 - t221 * t245) * qJD(3);
t94 = (-Icges(6,1) * t206 - t344) * t307 + (Icges(6,5) * t223 - t221 * t253) * qJD(3);
t226 = t223 * t207 * t94 + t134 * t310 + t221 * t92 - t312 * t385;
t93 = (-Icges(6,2) * t207 - t345) * t307 + (Icges(6,6) * t223 - t221 * t249) * qJD(3);
t354 = t206 * t93;
t53 = t134 * t221 + t223 * t385;
t357 = ((t380 * qJD(5) - t354) * t223 + t226) * t221 + t53 * t310;
t355 = rSges(5,3) * t223;
t353 = t222 * rSges(4,3);
t212 = t224 * rSges(4,3);
t203 = pkin(3) * t330;
t323 = -qJ(4) - t220;
t281 = t323 * t223;
t305 = pkin(4) * t327;
t349 = t305 + t203 + (t281 - t336) * t224 + t91;
t348 = (t281 + t382) * qJD(3) + t95;
t335 = t204 * t223;
t334 = t218 * t221;
t332 = t219 * t221;
t322 = t325 * t388 + t381;
t321 = t221 * t323 - t223 * t358 + t138;
t163 = qJD(4) * t221 + (qJ(4) * t223 - t359) * qJD(3);
t185 = pkin(3) * t223 + qJ(4) * t221;
t320 = t222 * t163 + t185 * t313;
t317 = t224 * pkin(1) + t222 * qJ(2);
t316 = t216 + t217;
t150 = Icges(4,3) * t224 + t222 * t247;
t315 = qJD(1) * t150;
t308 = qJD(4) * t223;
t97 = Icges(5,5) * t165 + Icges(5,6) * t164 - Icges(5,3) * t325;
t304 = t97 * t325;
t166 = t218 * t330 + t326;
t234 = t219 * t330 - t327;
t98 = -Icges(5,5) * t234 + Icges(5,6) * t166 + Icges(5,3) * t324;
t303 = t98 * t325;
t302 = t97 * t324;
t301 = t98 * t324;
t86 = Icges(6,4) * t147 + Icges(6,2) * t146 - Icges(6,6) * t325;
t88 = Icges(6,1) * t147 + Icges(6,4) * t146 - Icges(6,5) * t325;
t269 = t206 * t86 - t207 * t88;
t84 = Icges(6,5) * t147 + Icges(6,6) * t146 - Icges(6,3) * t325;
t32 = t221 * t84 - t223 * t269;
t49 = -t134 * t325 + t135 * t146 + t136 * t147;
t300 = t49 / 0.2e1 + t32 / 0.2e1;
t87 = Icges(6,4) * t149 + Icges(6,2) * t148 + Icges(6,6) * t324;
t89 = Icges(6,1) * t149 + Icges(6,4) * t148 + Icges(6,5) * t324;
t268 = t206 * t87 - t207 * t89;
t85 = Icges(6,5) * t149 + Icges(6,6) * t148 + Icges(6,3) * t324;
t33 = t221 * t85 - t223 * t268;
t50 = t134 * t324 + t148 * t135 + t149 * t136;
t299 = -t50 / 0.2e1 - t33 / 0.2e1;
t124 = -qJD(1) * t166 - t218 * t285;
t125 = qJD(1) * t234 + t219 * t285;
t297 = -t125 * rSges(5,1) - t124 * rSges(5,2) - rSges(5,3) * t287;
t294 = pkin(3) * t387 + qJ(4) * t287;
t293 = pkin(3) * t284 + qJ(4) * t386;
t156 = rSges(4,1) * t331 + rSges(4,2) * t325 + t212;
t291 = t224 * pkin(6) + t317;
t283 = t388 * t223;
t179 = t275 * qJD(3);
t282 = t179 * t316;
t280 = qJD(1) * t321;
t279 = -pkin(4) * t218 + t366;
t122 = qJD(1) * t164 + t218 * t284;
t123 = qJD(1) * t165 - t219 * t284;
t274 = -t123 * rSges(5,1) - t122 * rSges(5,2);
t273 = rSges(5,1) * t234 - t166 * rSges(5,2);
t272 = rSges(5,1) * t219 - rSges(5,2) * t218;
t24 = (qJD(1) * t279 - t220 * t312 - t308) * t222 + t318 + t383;
t25 = t209 + (-t308 + (t221 * t350 + t335) * qJD(3)) * t224 + (t279 * t224 + (-qJ(2) - t377) * t222) * qJD(1) - t276;
t266 = t222 * t25 - t224 * t24;
t26 = t146 * t86 + t147 * t88 - t325 * t84;
t27 = t146 * t87 + t147 * t89 - t325 * t85;
t17 = t27 * t222 + t224 * t26;
t265 = t222 * t26 - t224 * t27;
t28 = t148 * t86 + t149 * t88 + t324 * t84;
t29 = t148 * t87 + t149 * t89 + t324 * t85;
t18 = t29 * t222 + t224 * t28;
t264 = t222 * t28 - t224 * t29;
t171 = t185 * t314;
t30 = t171 + t222 * t280 + (-t163 - t348) * t224;
t31 = t222 * t348 + t224 * t280 + t320;
t263 = t31 * t222 - t224 * t30;
t262 = t33 * t222 + t32 * t224;
t261 = t222 * t32 - t224 * t33;
t225 = t222 * t366 - t224 * t283;
t38 = qJD(1) * t225 - t222 * t308 + t294 - t297 + t318;
t39 = t209 + (rSges(5,3) * t312 - t308) * t224 + (t366 * t224 + (-qJ(2) + t355 - t359) * t222) * qJD(1) + t274 + t293;
t260 = t222 * t39 - t224 * t38;
t137 = (-t221 * t272 + t355) * qJD(3);
t145 = rSges(5,3) * t221 + t223 * t272;
t56 = t145 * t314 + t171 + (-t137 - t163) * t224;
t57 = t222 * t137 + t145 * t313 + t320;
t259 = t57 * t222 - t224 * t56;
t257 = t222 * t91 + t224 * t90;
t256 = Icges(4,1) * t223 - t347;
t252 = -Icges(4,2) * t221 + t346;
t250 = Icges(5,4) * t219 - Icges(5,2) * t218;
t248 = Icges(4,5) * t223 - Icges(4,6) * t221;
t246 = Icges(5,5) * t219 - Icges(5,6) * t218;
t241 = -t221 * t374 - t223 * t375;
t240 = (t368 + t367) * t312;
t235 = rSges(3,3) * t224 + t222 * t361;
t231 = t241 * t222;
t230 = qJD(3) * t256;
t229 = qJD(3) * t252;
t211 = t224 * qJ(2);
t175 = t222 * t185;
t170 = qJ(4) * t324 - t203;
t160 = t224 * t170;
t159 = -rSges(3,2) * t224 + t222 * rSges(3,3) + t317;
t158 = t211 + t235;
t157 = t353 - t233;
t142 = Icges(5,6) * t221 + t223 * t250;
t133 = (Icges(5,5) * t223 - t221 * t254) * qJD(3);
t132 = (Icges(5,6) * t223 - t221 * t250) * qJD(3);
t128 = t209 + (t361 * t224 + (-rSges(3,3) - qJ(2)) * t222) * qJD(1);
t127 = qJD(1) * t235 + t318;
t118 = t291 + t156;
t117 = t222 * t306 + t211 + t233;
t115 = (-t145 - t185) * t224;
t114 = t145 * t222 + t175;
t109 = qJD(1) * t376 + t248 * t311;
t108 = -t248 * t309 + t315;
t107 = (-qJ(4) * t313 - qJD(4) * t222) * t223 + t294;
t104 = rSges(5,3) * t324 - t273;
t102 = -Icges(5,1) * t234 + Icges(5,4) * t166 + Icges(5,5) * t324;
t101 = Icges(5,1) * t165 + Icges(5,4) * t164 - Icges(5,5) * t325;
t100 = -Icges(5,4) * t234 + Icges(5,2) * t166 + Icges(5,6) * t324;
t99 = Icges(5,4) * t165 + Icges(5,2) * t164 - Icges(5,6) * t325;
t96 = t224 * (qJD(1) * t202 + t224 * t308 - t293);
t73 = (-t185 - t321) * t224;
t72 = t222 * t321 + t175;
t71 = -t222 * t283 + t291 - t381;
t70 = t203 + t211 + t225 + t273;
t69 = -t222 * t376 - t224 * t241;
t68 = t222 * t150 - t232;
t67 = -t224 * t376 + t231;
t66 = t150 * t224 + t242 * t222;
t65 = Icges(5,1) * t125 + Icges(5,4) * t124 + Icges(5,5) * t228;
t64 = Icges(5,1) * t123 + Icges(5,4) * t122 - Icges(5,5) * t386;
t63 = Icges(5,4) * t125 + Icges(5,2) * t124 + Icges(5,6) * t228;
t62 = Icges(5,4) * t123 + Icges(5,2) * t122 - Icges(5,6) * t386;
t55 = t291 - t298;
t54 = t279 * t222 + t224 * t377 + t211 + t271;
t51 = t257 * t223;
t48 = t104 * t224 + t222 * t322 + t160;
t45 = Icges(6,1) * t82 + Icges(6,4) * t81 + Icges(6,5) * t228;
t44 = Icges(6,1) * t80 + Icges(6,4) * t79 - Icges(6,5) * t386;
t43 = Icges(6,4) * t82 + Icges(6,2) * t81 + Icges(6,6) * t228;
t42 = Icges(6,4) * t80 + Icges(6,2) * t79 - Icges(6,6) * t386;
t41 = Icges(6,5) * t82 + Icges(6,6) * t81 + Icges(6,3) * t228;
t40 = Icges(6,5) * t80 + Icges(6,6) * t79 - Icges(6,3) * t386;
t37 = t166 * t100 - t102 * t234 + t301;
t36 = -t101 * t234 + t166 * t99 + t302;
t35 = t100 * t164 + t102 * t165 - t303;
t34 = t101 * t165 + t164 * t99 - t304;
t23 = t222 * t298 + t224 * t349 + t160;
t19 = t96 + t224 * (-rSges(5,3) * t286 - t274) + (-t107 + t297) * t222 + (t322 * t224 + (-t104 - t170) * t222) * qJD(1);
t16 = t134 * t228 + t81 * t135 + t82 * t136 + t146 * t93 + t147 * t94 - t325 * t92;
t15 = -t134 * t386 + t79 * t135 + t80 * t136 + t148 * t93 + t149 * t94 + t324 * t92;
t14 = t257 * t312 + (-t222 * t46 - t224 * t47 + (t222 * t90 - t224 * t91) * qJD(1)) * t223;
t13 = t50 * t221 - t223 * t264;
t12 = t49 * t221 - t223 * t265;
t11 = (qJD(3) * t268 + t40) * t221 + (qJD(3) * t85 - t206 * t42 + t207 * t44 + (-t206 * t89 - t207 * t87) * qJD(5)) * t223;
t10 = (qJD(3) * t269 + t41) * t221 + (qJD(3) * t84 - t206 * t43 + t207 * t45 + (-t206 * t88 - t207 * t86) * qJD(5)) * t223;
t9 = -t85 * t288 + t146 * t42 + t147 * t44 + t81 * t87 + t82 * t89 + (-t223 * t40 + t312 * t85) * t222;
t8 = -t84 * t288 + t146 * t43 + t147 * t45 + t81 * t86 + t82 * t88 + (-t223 * t41 + t312 * t84) * t222;
t7 = -t85 * t286 + t148 * t42 + t149 * t44 + t79 * t87 + t80 * t89 + (t224 * t40 - t314 * t85) * t223;
t6 = -t84 * t286 + t148 * t43 + t149 * t45 + t79 * t86 + t80 * t88 + (t224 * t41 - t314 * t84) * t223;
t5 = t96 + (t46 + (t220 * t221 - t335) * t309 + t293) * t224 + (t220 * t287 - t107 + t294 - t383) * t222 + ((t298 + t199) * t224 + (t305 - t170 + ((-qJ(4) + t220) * t223 - t382) * t224 - t349) * t222) * qJD(1);
t4 = -qJD(1) * t265 + t9 * t222 + t224 * t8;
t3 = -qJD(1) * t264 + t7 * t222 + t224 * t6;
t2 = (qJD(3) * t265 + t16) * t221 + (-qJD(1) * t17 + qJD(3) * t49 - t222 * t8 + t224 * t9) * t223;
t1 = (qJD(3) * t264 + t15) * t221 + (-qJD(1) * t18 + qJD(3) * t50 - t222 * t6 + t224 * t7) * t223;
t20 = [-t255 * t310 + t226 + (t24 * t55 + t25 * t54) * t369 + (t38 * t71 + t39 * t70) * t370 + (t117 * t75 + t118 * t74) * t371 + 0.2e1 * m(3) * (t127 * t159 + t128 * t158) + t380 * t307 + (Icges(5,3) * t310 - t230) * t221 + (t218 * t142 - t219 * t143 - t221 * t246 + t251) * t312 + (Icges(5,3) * t312 - t218 * t132 + t219 * t133 + t246 * t310 - t229 - t354) * t223; m(6) * ((t222 * t55 + t224 * t54) * qJD(1) + t266) + m(5) * ((t222 * t71 + t224 * t70) * qJD(1) + t260) + m(4) * ((t117 * t224 + t118 * t222) * qJD(1) + t379) + m(3) * (-t127 * t224 + t222 * t128 + (t158 * t224 + t159 * t222) * qJD(1)); 0; m(4) * (t379 * t186 - (t117 * t222 - t118 * t224) * t179) + m(5) * (t114 * t39 + t115 * t38 + t56 * t71 + t57 * t70) + m(6) * (t24 * t73 + t25 * t72 + t30 * t55 + t31 * t54) + ((t124 * t391 + t125 * t392 + t228 * t390 + t229 * t364 + t375 * t389) * t224 + (t152 * t389 + t252 * t309 / 0.2e1 + t123 * t392 + t122 * t391 - t386 * t390) * t222) * t221 + ((t118 * t360 - t142 * t164 / 0.2e1 + t165 * t384 + (-t97 / 0.2e1 + t152 / 0.2e1) * t221 + (-t101 * t219 / 0.2e1 + t218 * t99 / 0.2e1) * t223 - t300) * t222 + (t166 * t142 / 0.2e1 + t234 * t384 + t117 * t360 + (t98 / 0.2e1 + t375 / 0.2e1) * t221 + (-t100 * t218 / 0.2e1 + t102 * t219 / 0.2e1 - t374 / 0.2e1) * t223 - t299) * t224) * qJD(1) + (-(t217 / 0.2e1 + t216 / 0.2e1) * t247 + t241 * t364 - t232 / 0.2e1) * qJD(3) + (t11 + t15 + t122 * t142 + t123 * t143 + t166 * t132 - t234 * t133 + (-t218 * t62 + t219 * t64 - t256 * t309) * t223 + (t100 * t334 - t102 * t332 + t223 * t98) * qJD(3)) * t222 / 0.2e1 + (t10 + t16 + t124 * t142 + t125 * t143 + t164 * t132 + t165 * t133 + (qJD(1) * t374 - t218 * t63 + t219 * t65 + t222 * t230) * t223 + (-t101 * t332 + t223 * t97 + t334 * t99) * qJD(3)) * t362; m(5) * ((t114 * t224 + t115 * t222) * qJD(1) + t259) + m(6) * ((t222 * t73 + t224 * t72) * qJD(1) + t263) - m(4) * t282; t224 * t4 + (t23 * t5 + t30 * t73 + t31 * t72) * t369 + t222 * t3 + (t114 * t57 + t115 * t56 + t19 * t48) * t370 + t224 * ((t224 * t109 + (t67 + t232) * qJD(1)) * t224 + (-t66 * qJD(1) + (-t310 * t374 + t312 * t375) * t222 + (t108 + (-t152 * t221 + t154 * t223) * qJD(3) + (-t150 + t241) * qJD(1)) * t224) * t222) + t224 * ((t125 * t101 + t124 * t99 + t164 * t63 + t165 * t65 + (t35 - t302) * qJD(1)) * t224 + (t124 * t100 + t125 * t102 + t164 * t62 + t165 * t64 + (-t34 - t301) * qJD(1)) * t222) + ((-t222 * t156 + t157 * t224) * (-t222 * t292 + (-t186 * t217 + t216 * t356) * qJD(3) + ((-t156 + t212) * t224 + (-t157 + t233 + t353) * t222) * qJD(1)) - t186 * t282) * t371 + t222 * ((t222 * t108 + (-t68 + t231) * qJD(1)) * t222 + (t69 * qJD(1) + (t152 * t312 - t154 * t310 + t315) * t224 + (t109 + (-t221 * t375 + t223 * t374) * qJD(3) + t242 * qJD(1)) * t222) * t224) + t222 * ((t122 * t100 + t123 * t102 + t166 * t62 - t234 * t64 + (-t36 - t303) * qJD(1)) * t222 + (t123 * t101 + t122 * t99 + t166 * t63 - t234 * t65 + (t37 - t304) * qJD(1)) * t224) + (-t17 + (-t34 - t66) * t224 + (-t35 - t67) * t222) * t314 + (t18 + (t36 + t68) * t224 + (t37 + t69) * t222) * t313; 0.2e1 * ((t222 * t54 - t224 * t55) * t367 + (t222 * t70 - t224 * t71) * t368) * t312 + 0.2e1 * ((-t313 * t54 - t314 * t55 - t266) * t367 + (-t313 * t70 - t314 * t71 - t260) * t368) * t223; 0.2e1 * t316 * t240; 0.2e1 * ((-t309 * t73 + t311 * t72 + t5) * t367 + (t114 * t311 - t115 * t309 + t19) * t368) * t221 + 0.2e1 * ((qJD(3) * t23 - t313 * t72 - t314 * t73 - t263) * t367 + (qJD(3) * t48 - t114 * t313 - t115 * t314 - t259) * t368) * t223; 0.4e1 * (0.1e1 - t316) * t223 * t240; m(6) * (t21 * t54 + t22 * t55 + t24 * t58 + t25 * t59) + (t222 * t300 + t224 * t299) * t312 + ((t15 / 0.2e1 + t11 / 0.2e1) * t224 + (-t16 / 0.2e1 - t10 / 0.2e1) * t222 + (t222 * t299 - t224 * t300) * qJD(1)) * t223 + t357; m(6) * t372; m(6) * (t14 * t23 + t21 * t72 + t22 * t73 + t30 * t58 + t31 * t59 - t5 * t51) + (t2 / 0.2e1 + (t33 * qJD(1) + t10) * t365 - t18 * t312 / 0.2e1 + qJD(1) * t13 / 0.2e1) * t224 + (t12 * t389 + (-qJD(1) * t32 + t11) * t365 + t17 * t312 / 0.2e1 + t1 / 0.2e1) * t222 + (qJD(3) * t262 / 0.2e1 + t3 * t362 + t4 * t364 + (t18 * t364 - t224 * t17 / 0.2e1) * qJD(1)) * t223; m(6) * ((t14 + (t222 * t59 - t224 * t58) * qJD(3)) * t221 + (-qJD(3) * t51 - t372) * t223); (-t14 * t51 + t21 * t59 + t22 * t58) * t369 + ((t222 * t12 - t224 * t13 + t221 * t261) * qJD(3) + t357) * t221 + (-t222 * t2 + t224 * t1 + t221 * (-t10 * t222 + t11 * t224) + (t53 * t221 - t223 * t261) * qJD(3) + (-t224 * t12 - t222 * t13 - t221 * t262) * qJD(1)) * t223;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t20(1), t20(2), t20(4), t20(7), t20(11); t20(2), t20(3), t20(5), t20(8), t20(12); t20(4), t20(5), t20(6), t20(9), t20(13); t20(7), t20(8), t20(9), t20(10), t20(14); t20(11), t20(12), t20(13), t20(14), t20(15);];
Mq = res;
