% Calculate time derivative of joint inertia matrix for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRPR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_inertiaDJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR2_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRPR2_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRPR2_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:49:10
% EndTime: 2019-03-08 18:49:32
% DurationCPUTime: 16.20s
% Computational Cost: add. (89542->904), mult. (263835->1199), div. (0->0), fcn. (326970->14), ass. (0->370)
t327 = sin(pkin(11));
t330 = cos(pkin(6));
t326 = sin(pkin(12));
t329 = cos(pkin(11));
t398 = t329 * t326;
t402 = cos(pkin(12));
t320 = t327 * t402 + t330 * t398;
t332 = sin(qJ(3));
t359 = t330 * t402;
t399 = t327 * t326;
t340 = -t329 * t359 + t399;
t403 = cos(pkin(7));
t337 = t340 * t403;
t328 = sin(pkin(6));
t401 = sin(pkin(7));
t360 = t328 * t401;
t406 = cos(qJ(3));
t302 = t320 * t406 + (-t329 * t360 - t337) * t332;
t361 = t328 * t403;
t314 = -t329 * t361 + t340 * t401;
t404 = sin(qJ(4));
t405 = cos(qJ(4));
t279 = t302 * t404 - t314 * t405;
t280 = t302 * t405 + t314 * t404;
t345 = t406 * t401;
t343 = t328 * t345;
t301 = t320 * t332 + t329 * t343 + t337 * t406;
t186 = Icges(6,5) * t301 - Icges(6,6) * t280 + Icges(6,3) * t279;
t192 = Icges(5,4) * t280 - Icges(5,2) * t279 + Icges(5,6) * t301;
t448 = t186 - t192;
t321 = t329 * t402 - t330 * t399;
t339 = t327 * t359 + t398;
t336 = t339 * t403;
t304 = t321 * t406 + (t327 * t360 - t336) * t332;
t315 = t327 * t361 + t339 * t401;
t281 = t304 * t404 - t315 * t405;
t282 = t304 * t405 + t315 * t404;
t303 = t321 * t332 - t327 * t343 + t336 * t406;
t187 = Icges(6,5) * t303 - Icges(6,6) * t282 + Icges(6,3) * t281;
t193 = Icges(5,4) * t282 - Icges(5,2) * t281 + Icges(5,6) * t303;
t447 = t187 - t193;
t188 = Icges(5,5) * t280 - Icges(5,6) * t279 + Icges(5,3) * t301;
t194 = Icges(6,1) * t301 - Icges(6,4) * t280 + Icges(6,5) * t279;
t446 = t188 + t194;
t189 = Icges(5,5) * t282 - Icges(5,6) * t281 + Icges(5,3) * t303;
t195 = Icges(6,1) * t303 - Icges(6,4) * t282 + Icges(6,5) * t281;
t445 = t189 + t195;
t190 = Icges(6,4) * t301 - Icges(6,2) * t280 + Icges(6,6) * t279;
t196 = Icges(5,1) * t280 - Icges(5,4) * t279 + Icges(5,5) * t301;
t444 = -t190 + t196;
t191 = Icges(6,4) * t303 - Icges(6,2) * t282 + Icges(6,6) * t281;
t197 = Icges(5,1) * t282 - Icges(5,4) * t281 + Icges(5,5) * t303;
t443 = -t191 + t197;
t344 = t403 * t402;
t313 = t330 * t401 * t332 + (t326 * t406 + t332 * t344) * t328;
t319 = t330 * t403 - t360 * t402;
t305 = t313 * t404 - t319 * t405;
t306 = t313 * t405 + t319 * t404;
t437 = (-t326 * t332 + t406 * t344) * t328 + t330 * t345;
t248 = -Icges(6,5) * t437 - Icges(6,6) * t306 + Icges(6,3) * t305;
t251 = Icges(5,4) * t306 - Icges(5,2) * t305 - Icges(5,6) * t437;
t442 = t248 - t251;
t249 = Icges(5,5) * t306 - Icges(5,6) * t305 - Icges(5,3) * t437;
t252 = -Icges(6,1) * t437 - Icges(6,4) * t306 + Icges(6,5) * t305;
t441 = t249 + t252;
t250 = -Icges(6,4) * t437 - Icges(6,2) * t306 + Icges(6,6) * t305;
t253 = Icges(5,1) * t306 - Icges(5,4) * t305 - Icges(5,5) * t437;
t440 = -t250 + t253;
t435 = t279 * t448 + t444 * t280 + t446 * t301;
t434 = t279 * t447 + t280 * t443 + t301 * t445;
t433 = t281 * t448 + t444 * t282 + t446 * t303;
t432 = t281 * t447 + t282 * t443 + t303 * t445;
t309 = t313 * qJD(3);
t427 = t305 * t442 + t306 * t440 - t437 * t441;
t439 = t427 * t309;
t431 = t305 * t448 + t444 * t306 - t446 * t437;
t430 = t305 * t447 + t306 * t443 - t437 * t445;
t429 = t279 * t442 + t280 * t440 + t301 * t441;
t428 = t281 * t442 + t282 * t440 + t303 * t441;
t419 = m(7) / 0.2e1;
t372 = t419 + m(6) / 0.2e1;
t438 = 0.2e1 * t372;
t436 = 2 * m(4);
t297 = t301 * qJD(3);
t229 = qJD(4) * t280 - t297 * t404;
t362 = qJD(4) * t404;
t363 = qJD(4) * t405;
t230 = -t297 * t405 - t302 * t362 + t314 * t363;
t149 = pkin(4) * t230 + qJ(5) * t229 + qJD(5) * t279;
t225 = pkin(4) * t280 + qJ(5) * t279;
t300 = t304 * qJD(3);
t394 = t303 * t149 + t300 * t225;
t308 = t437 * qJD(3);
t277 = t306 * qJD(4) + t308 * t404;
t298 = t302 * qJD(3);
t331 = sin(qJ(6));
t333 = cos(qJ(6));
t234 = t279 * t331 + t301 * t333;
t159 = -qJD(6) * t234 + t229 * t333 - t298 * t331;
t233 = t279 * t333 - t301 * t331;
t160 = qJD(6) * t233 + t229 * t331 + t298 * t333;
t105 = Icges(7,5) * t160 + Icges(7,6) * t159 + Icges(7,3) * t230;
t107 = Icges(7,4) * t160 + Icges(7,2) * t159 + Icges(7,6) * t230;
t109 = Icges(7,1) * t160 + Icges(7,4) * t159 + Icges(7,5) * t230;
t151 = Icges(7,5) * t234 + Icges(7,6) * t233 + Icges(7,3) * t280;
t153 = Icges(7,4) * t234 + Icges(7,2) * t233 + Icges(7,6) * t280;
t155 = Icges(7,1) * t234 + Icges(7,4) * t233 + Icges(7,5) * t280;
t22 = t105 * t280 + t107 * t233 + t109 * t234 + t151 * t230 + t153 * t159 + t155 * t160;
t299 = t303 * qJD(3);
t231 = qJD(4) * t282 - t299 * t404;
t236 = t281 * t331 + t303 * t333;
t161 = -qJD(6) * t236 + t231 * t333 - t300 * t331;
t235 = t281 * t333 - t303 * t331;
t162 = qJD(6) * t235 + t231 * t331 + t300 * t333;
t232 = -t299 * t405 - t304 * t362 + t315 * t363;
t106 = Icges(7,5) * t162 + Icges(7,6) * t161 + Icges(7,3) * t232;
t108 = Icges(7,4) * t162 + Icges(7,2) * t161 + Icges(7,6) * t232;
t110 = Icges(7,1) * t162 + Icges(7,4) * t161 + Icges(7,5) * t232;
t152 = Icges(7,5) * t236 + Icges(7,6) * t235 + Icges(7,3) * t282;
t154 = Icges(7,4) * t236 + Icges(7,2) * t235 + Icges(7,6) * t282;
t156 = Icges(7,1) * t236 + Icges(7,4) * t235 + Icges(7,5) * t282;
t23 = t106 * t280 + t108 * t233 + t110 * t234 + t152 * t230 + t154 * t159 + t156 * t160;
t284 = t305 * t331 - t333 * t437;
t209 = -qJD(6) * t284 + t277 * t333 - t309 * t331;
t283 = t305 * t333 + t331 * t437;
t210 = qJD(6) * t283 + t277 * t331 + t309 * t333;
t278 = t308 * t405 - t313 * t362 + t319 * t363;
t134 = Icges(7,5) * t210 + Icges(7,6) * t209 + Icges(7,3) * t278;
t135 = Icges(7,4) * t210 + Icges(7,2) * t209 + Icges(7,6) * t278;
t136 = Icges(7,1) * t210 + Icges(7,4) * t209 + Icges(7,5) * t278;
t202 = Icges(7,5) * t284 + Icges(7,6) * t283 + Icges(7,3) * t306;
t203 = Icges(7,4) * t284 + Icges(7,2) * t283 + Icges(7,6) * t306;
t204 = Icges(7,1) * t284 + Icges(7,4) * t283 + Icges(7,5) * t306;
t40 = t134 * t280 + t135 * t233 + t136 * t234 + t159 * t203 + t160 * t204 + t202 * t230;
t71 = t151 * t280 + t153 * t233 + t155 * t234;
t72 = t152 * t280 + t154 * t233 + t156 * t234;
t89 = t202 * t280 + t203 * t233 + t204 * t234;
t3 = t22 * t301 + t23 * t303 + t298 * t71 + t300 * t72 + t309 * t89 - t40 * t437;
t163 = Icges(6,5) * t298 - Icges(6,6) * t230 + Icges(6,3) * t229;
t167 = Icges(6,4) * t298 - Icges(6,2) * t230 + Icges(6,6) * t229;
t171 = Icges(6,1) * t298 - Icges(6,4) * t230 + Icges(6,5) * t229;
t50 = t163 * t279 - t167 * t280 + t171 * t301 + t186 * t229 - t190 * t230 + t194 * t298;
t164 = Icges(6,5) * t300 - Icges(6,6) * t232 + Icges(6,3) * t231;
t168 = Icges(6,4) * t300 - Icges(6,2) * t232 + Icges(6,6) * t231;
t172 = Icges(6,1) * t300 - Icges(6,4) * t232 + Icges(6,5) * t231;
t51 = t164 * t279 - t168 * t280 + t172 * t301 + t187 * t229 - t191 * t230 + t195 * t298;
t165 = Icges(5,5) * t230 - Icges(5,6) * t229 + Icges(5,3) * t298;
t169 = Icges(5,4) * t230 - Icges(5,2) * t229 + Icges(5,6) * t298;
t173 = Icges(5,1) * t230 - Icges(5,4) * t229 + Icges(5,5) * t298;
t54 = t165 * t301 - t169 * t279 + t173 * t280 + t188 * t298 - t192 * t229 + t196 * t230;
t166 = Icges(5,5) * t232 - Icges(5,6) * t231 + Icges(5,3) * t300;
t170 = Icges(5,4) * t232 - Icges(5,2) * t231 + Icges(5,6) * t300;
t174 = Icges(5,1) * t232 - Icges(5,4) * t231 + Icges(5,5) * t300;
t55 = t166 * t301 - t170 * t279 + t174 * t280 + t189 * t298 - t193 * t229 + t197 * t230;
t211 = Icges(6,5) * t309 - Icges(6,6) * t278 + Icges(6,3) * t277;
t213 = Icges(6,4) * t309 - Icges(6,2) * t278 + Icges(6,6) * t277;
t215 = Icges(6,1) * t309 - Icges(6,4) * t278 + Icges(6,5) * t277;
t67 = t211 * t279 - t213 * t280 + t215 * t301 + t229 * t248 - t230 * t250 + t252 * t298;
t212 = Icges(5,5) * t278 - Icges(5,6) * t277 + Icges(5,3) * t309;
t214 = Icges(5,4) * t278 - Icges(5,2) * t277 + Icges(5,6) * t309;
t216 = Icges(5,1) * t278 - Icges(5,4) * t277 + Icges(5,5) * t309;
t69 = t212 * t301 - t214 * t279 + t216 * t280 - t229 * t251 + t230 * t253 + t249 * t298;
t425 = t3 - (t67 + t69) * t437 + t429 * t309 + (t51 + t55) * t303 + (t50 + t54) * t301 + t434 * t300 + t435 * t298;
t24 = t105 * t282 + t107 * t235 + t109 * t236 + t151 * t232 + t153 * t161 + t155 * t162;
t25 = t106 * t282 + t108 * t235 + t110 * t236 + t152 * t232 + t154 * t161 + t156 * t162;
t41 = t134 * t282 + t135 * t235 + t136 * t236 + t161 * t203 + t162 * t204 + t202 * t232;
t73 = t151 * t282 + t153 * t235 + t155 * t236;
t74 = t152 * t282 + t154 * t235 + t156 * t236;
t90 = t202 * t282 + t203 * t235 + t204 * t236;
t4 = t24 * t301 + t25 * t303 + t298 * t73 + t300 * t74 + t309 * t90 - t41 * t437;
t52 = t163 * t281 - t167 * t282 + t171 * t303 + t186 * t231 - t190 * t232 + t194 * t300;
t53 = t164 * t281 - t168 * t282 + t172 * t303 + t187 * t231 - t191 * t232 + t195 * t300;
t56 = t165 * t303 - t169 * t281 + t173 * t282 + t188 * t300 - t192 * t231 + t196 * t232;
t57 = t166 * t303 - t170 * t281 + t174 * t282 + t189 * t300 - t193 * t231 + t197 * t232;
t68 = t211 * t281 - t213 * t282 + t215 * t303 + t231 * t248 - t232 * t250 + t252 * t300;
t70 = t212 * t303 - t214 * t281 + t216 * t282 - t231 * t251 + t232 * t253 + t249 * t300;
t424 = t4 - (t68 + t70) * t437 + t428 * t309 + (t53 + t57) * t303 + (t52 + t56) * t301 + t432 * t300 + t433 * t298;
t104 = t202 * t306 + t203 * t283 + t204 * t284;
t27 = t105 * t306 + t107 * t283 + t109 * t284 + t151 * t278 + t153 * t209 + t155 * t210;
t28 = t106 * t306 + t108 * t283 + t110 * t284 + t152 * t278 + t154 * t209 + t156 * t210;
t43 = t134 * t306 + t135 * t283 + t136 * t284 + t202 * t278 + t203 * t209 + t204 * t210;
t82 = t151 * t306 + t153 * t283 + t155 * t284;
t83 = t152 * t306 + t154 * t283 + t156 * t284;
t6 = t104 * t309 + t27 * t301 + t28 * t303 + t298 * t82 + t300 * t83 - t43 * t437;
t60 = t163 * t305 - t167 * t306 - t171 * t437 + t186 * t277 - t190 * t278 + t194 * t309;
t61 = t164 * t305 - t168 * t306 - t172 * t437 + t187 * t277 - t191 * t278 + t195 * t309;
t62 = -t165 * t437 - t169 * t305 + t173 * t306 + t188 * t309 - t192 * t277 + t196 * t278;
t63 = -t166 * t437 - t170 * t305 + t174 * t306 + t189 * t309 - t193 * t277 + t197 * t278;
t79 = t211 * t305 - t213 * t306 - t215 * t437 + t248 * t277 - t250 * t278 + t252 * t309;
t80 = -t212 * t437 - t214 * t305 + t216 * t306 + t249 * t309 - t251 * t277 + t253 * t278;
t423 = t6 - (t79 + t80) * t437 + t439 + (t61 + t63) * t303 + (t60 + t62) * t301 + t430 * t300 + t431 * t298;
t422 = -0.2e1 * t314;
t421 = 0.2e1 * t330;
t418 = t230 / 0.2e1;
t417 = t232 / 0.2e1;
t416 = t278 / 0.2e1;
t415 = t280 / 0.2e1;
t414 = t282 / 0.2e1;
t413 = t298 / 0.2e1;
t412 = t300 / 0.2e1;
t411 = t306 / 0.2e1;
t410 = t309 / 0.2e1;
t409 = t314 / 0.2e1;
t408 = t315 / 0.2e1;
t407 = t319 / 0.2e1;
t111 = rSges(7,1) * t160 + rSges(7,2) * t159 + rSges(7,3) * t230;
t397 = pkin(5) * t298 + pkin(10) * t230 + t111;
t112 = rSges(7,1) * t162 + rSges(7,2) * t161 + rSges(7,3) * t232;
t396 = pkin(5) * t300 + pkin(10) * t232 + t112;
t137 = rSges(7,1) * t210 + rSges(7,2) * t209 + rSges(7,3) * t278;
t395 = pkin(5) * t309 + pkin(10) * t278 + t137;
t150 = pkin(4) * t232 + qJ(5) * t231 + qJD(5) * t281;
t226 = pkin(4) * t282 + qJ(5) * t281;
t393 = -t150 * t437 + t309 * t226;
t146 = t315 * t149;
t269 = -pkin(3) * t297 + pkin(9) * t298;
t239 = t315 * t269;
t392 = t146 + t239;
t270 = -pkin(3) * t299 + pkin(9) * t300;
t258 = t319 * t270;
t391 = t319 * t150 + t258;
t175 = rSges(6,1) * t298 - rSges(6,2) * t230 + rSges(6,3) * t229;
t390 = -t149 - t175;
t176 = rSges(6,1) * t300 - rSges(6,2) * t232 + rSges(6,3) * t231;
t389 = -t150 - t176;
t157 = rSges(7,1) * t234 + rSges(7,2) * t233 + rSges(7,3) * t280;
t388 = pkin(5) * t301 + pkin(10) * t280 + t157;
t158 = rSges(7,1) * t236 + rSges(7,2) * t235 + rSges(7,3) * t282;
t387 = pkin(5) * t303 + pkin(10) * t282 + t158;
t182 = pkin(4) * t278 + qJ(5) * t277 + qJD(5) * t305;
t273 = pkin(4) * t306 + qJ(5) * t305;
t386 = t301 * t182 + t298 * t273;
t293 = pkin(3) * t308 + pkin(9) * t309;
t274 = t314 * t293;
t385 = t314 * t182 + t274;
t217 = rSges(6,1) * t309 - rSges(6,2) * t278 + rSges(6,3) * t277;
t384 = -t182 - t217;
t198 = rSges(6,1) * t301 - rSges(6,2) * t280 + rSges(6,3) * t279;
t383 = -t198 - t225;
t199 = rSges(6,1) * t303 - rSges(6,2) * t282 + rSges(6,3) * t281;
t382 = -t199 - t226;
t207 = rSges(7,1) * t284 + rSges(7,2) * t283 + rSges(7,3) * t306;
t381 = -pkin(5) * t437 + pkin(10) * t306 + t207;
t271 = pkin(3) * t302 + pkin(9) * t301;
t256 = t315 * t271;
t380 = t315 * t225 + t256;
t272 = pkin(3) * t304 + pkin(9) * t303;
t260 = t319 * t272;
t379 = t319 * t226 + t260;
t378 = t270 * t422 + 0.2e1 * t239;
t254 = -rSges(6,1) * t437 - rSges(6,2) * t306 + rSges(6,3) * t305;
t377 = -t254 - t273;
t294 = pkin(3) * t313 - pkin(9) * t437;
t275 = t314 * t294;
t376 = t314 * t273 + t275;
t373 = 0.2e1 * t328;
t371 = -t149 - t397;
t370 = -t150 - t396;
t369 = -t182 - t395;
t368 = -t225 - t388;
t367 = -t226 - t387;
t366 = -t273 - t381;
t358 = 2 * m(5);
t357 = 0.2e1 * m(6);
t355 = 0.2e1 * m(7);
t26 = t298 * t367 + t300 * t388 + t301 * t370 + t303 * t397 + t394;
t45 = t175 * t303 + t198 * t300 + t298 * t382 + t301 * t389 + t394;
t352 = m(6) * t45 + m(7) * t26;
t35 = t298 * t381 + t301 * t395 + t309 * t368 - t371 * t437 + t386;
t58 = t217 * t301 + t254 * t298 + t309 * t383 - t390 * t437 + t386;
t351 = m(6) * t58 + m(7) * t35;
t36 = t300 * t366 + t303 * t369 + t309 * t387 - t396 * t437 + t393;
t59 = -t176 * t437 + t199 * t309 + t300 * t377 + t303 * t384 + t393;
t350 = m(6) * t59 + m(7) * t36;
t42 = t397 * t315 + (-t270 + t370) * t314 + t392;
t64 = t175 * t315 + (-t270 + t389) * t314 + t392;
t349 = m(6) * t64 + m(7) * t42;
t46 = t395 * t314 + (-t269 + t371) * t319 + t385;
t75 = t217 * t314 + (-t269 + t390) * t319 + t385;
t348 = m(6) * t75 + m(7) * t46;
t47 = t396 * t319 + (-t293 + t369) * t315 + t391;
t76 = t176 * t319 + (-t293 + t384) * t315 + t391;
t347 = m(6) * t76 + m(7) * t47;
t44 = t111 * t282 - t112 * t280 + t157 * t232 - t158 * t230;
t177 = rSges(5,1) * t230 - rSges(5,2) * t229 + rSges(5,3) * t298;
t335 = m(5) * t177 + m(6) * t175 + m(7) * t397;
t178 = rSges(5,1) * t232 - rSges(5,2) * t231 + rSges(5,3) * t300;
t334 = -m(5) * t178 - m(6) * t176 - m(7) * t396;
t292 = rSges(4,1) * t308 - rSges(4,2) * t309;
t291 = Icges(4,1) * t308 - Icges(4,4) * t309;
t290 = Icges(4,4) * t308 - Icges(4,2) * t309;
t289 = Icges(4,5) * t308 - Icges(4,6) * t309;
t288 = rSges(4,1) * t313 + rSges(4,2) * t437 + rSges(4,3) * t319;
t287 = Icges(4,1) * t313 + Icges(4,4) * t437 + Icges(4,5) * t319;
t286 = Icges(4,4) * t313 + Icges(4,2) * t437 + Icges(4,6) * t319;
t268 = -rSges(4,1) * t299 - rSges(4,2) * t300;
t267 = -rSges(4,1) * t297 - rSges(4,2) * t298;
t266 = -Icges(4,1) * t299 - Icges(4,4) * t300;
t265 = -Icges(4,1) * t297 - Icges(4,4) * t298;
t264 = -Icges(4,4) * t299 - Icges(4,2) * t300;
t263 = -Icges(4,4) * t297 - Icges(4,2) * t298;
t262 = -Icges(4,5) * t299 - Icges(4,6) * t300;
t261 = -Icges(4,5) * t297 - Icges(4,6) * t298;
t255 = rSges(5,1) * t306 - rSges(5,2) * t305 - rSges(5,3) * t437;
t247 = rSges(4,1) * t304 - rSges(4,2) * t303 + rSges(4,3) * t315;
t246 = rSges(4,1) * t302 - rSges(4,2) * t301 + rSges(4,3) * t314;
t245 = Icges(4,1) * t304 - Icges(4,4) * t303 + Icges(4,5) * t315;
t244 = Icges(4,1) * t302 - Icges(4,4) * t301 + Icges(4,5) * t314;
t243 = Icges(4,4) * t304 - Icges(4,2) * t303 + Icges(4,6) * t315;
t242 = Icges(4,4) * t302 - Icges(4,2) * t301 + Icges(4,6) * t314;
t228 = t301 * t273;
t222 = t437 * t226;
t218 = rSges(5,1) * t278 - rSges(5,2) * t277 + rSges(5,3) * t309;
t208 = t303 * t225;
t206 = t268 * t319 - t292 * t315;
t205 = -t267 * t319 + t292 * t314;
t201 = rSges(5,1) * t282 - rSges(5,2) * t281 + rSges(5,3) * t303;
t200 = rSges(5,1) * t280 - rSges(5,2) * t279 + rSges(5,3) * t301;
t179 = t267 * t315 - t268 * t314;
t139 = -t201 * t437 - t255 * t303;
t138 = t200 * t437 + t255 * t301;
t131 = t200 * t303 - t201 * t301;
t129 = t201 * t319 + t260 + (-t255 - t294) * t315;
t128 = t255 * t314 + t275 + (-t200 - t271) * t319;
t123 = t158 * t306 - t207 * t282;
t122 = -t157 * t306 + t207 * t280;
t121 = t200 * t315 + t256 + (-t201 - t272) * t314;
t120 = -t199 * t437 + t303 * t377 - t222;
t119 = t254 * t301 - t383 * t437 + t228;
t114 = t178 * t319 + t258 + (-t218 - t293) * t315;
t113 = t218 * t314 + t274 + (-t177 - t269) * t319;
t95 = t157 * t282 - t158 * t280;
t94 = t199 * t319 + (-t294 + t377) * t315 + t379;
t93 = t254 * t314 + (-t271 + t383) * t319 + t376;
t92 = t198 * t303 + t301 * t382 + t208;
t91 = t177 * t315 + t239 + (-t178 - t270) * t314;
t88 = -t178 * t437 + t201 * t309 - t218 * t303 - t255 * t300;
t87 = t177 * t437 - t200 * t309 + t218 * t301 + t255 * t298;
t86 = t198 * t315 + (-t272 + t382) * t314 + t380;
t85 = t303 * t366 - t387 * t437 - t222;
t84 = t301 * t381 - t368 * t437 + t228;
t81 = t177 * t303 - t178 * t301 + t200 * t300 - t201 * t298;
t78 = t387 * t319 + (-t294 + t366) * t315 + t379;
t77 = t381 * t314 + (-t271 + t368) * t319 + t376;
t66 = t301 * t367 + t303 * t388 + t208;
t65 = t388 * t315 + (-t272 + t367) * t314 + t380;
t49 = t112 * t306 - t137 * t282 + t158 * t278 - t207 * t232;
t48 = -t111 * t306 + t137 * t280 - t157 * t278 + t207 * t230;
t39 = t104 * t319 + t314 * t82 + t315 * t83;
t38 = -t104 * t437 + t301 * t82 + t303 * t83;
t37 = t104 * t306 + t280 * t82 + t282 * t83;
t34 = t314 * t73 + t315 * t74 + t319 * t90;
t33 = t314 * t71 + t315 * t72 + t319 * t89;
t32 = t301 * t73 + t303 * t74 - t437 * t90;
t31 = t301 * t71 + t303 * t72 - t437 * t89;
t30 = t280 * t73 + t282 * t74 + t306 * t90;
t29 = t280 * t71 + t282 * t72 + t306 * t89;
t21 = t314 * t62 + t315 * t63 + t319 * t80;
t20 = t314 * t60 + t315 * t61 + t319 * t79;
t19 = t314 * t56 + t315 * t57 + t319 * t70;
t18 = t314 * t54 + t315 * t55 + t319 * t69;
t17 = t314 * t52 + t315 * t53 + t319 * t68;
t16 = t314 * t50 + t315 * t51 + t319 * t67;
t9 = t27 * t314 + t28 * t315 + t319 * t43;
t8 = t24 * t314 + t25 * t315 + t319 * t41;
t7 = t22 * t314 + t23 * t315 + t319 * t40;
t5 = t104 * t278 + t230 * t82 + t232 * t83 + t27 * t280 + t28 * t282 + t306 * t43;
t2 = t230 * t73 + t232 * t74 + t24 * t280 + t25 * t282 + t278 * t90 + t306 * t41;
t1 = t22 * t280 + t23 * t282 + t230 * t71 + t232 * t72 + t278 * t89 + t306 * t40;
t10 = [0; 0; 0; m(5) * t378 / 0.2e1 + (m(4) * t267 + t335) * t315 + (-m(4) * t268 + t334) * t314 + t372 * (t150 * t422 + 0.2e1 * t146 + t378); (m(4) * t179 + m(5) * t91 + t349) * t330 + ((-m(4) * t206 - m(5) * t114 - t347) * t329 + (m(4) * t205 + m(5) * t113 + t348) * t327) * t328; (t42 * t65 + t46 * t77 + t47 * t78) * t355 + (t64 * t86 + t75 * t93 + t76 * t94) * t357 + (t113 * t128 + t114 * t129 + t121 * t91) * t358 + (t9 + t21 + t20 + (-t205 * t246 + t206 * t247) * t436 + (-t286 * t309 + t287 * t308 + t289 * t319 + t290 * t437 + t291 * t313) * t319) * t319 + (t8 + t19 + t17 + (t179 * t246 - t206 * t288) * t436 + (-t243 * t300 - t245 * t299 + t262 * t315 - t264 * t303 + t266 * t304) * t315 + (-t243 * t309 + t245 * t308 + t262 * t319 + t264 * t437 + t266 * t313 - t286 * t300 - t287 * t299 + t289 * t315 - t290 * t303 + t291 * t304) * t319) * t315 + (t7 + t16 + t18 + (-t179 * t247 + t205 * t288) * t436 + (-t242 * t298 - t244 * t297 + t261 * t314 - t263 * t301 + t265 * t302) * t314 + (-t242 * t309 + t244 * t308 + t261 * t319 + t263 * t437 + t265 * t313 - t286 * t298 - t287 * t297 + t289 * t314 - t290 * t301 + t291 * t302) * t319 + (-t242 * t300 - t243 * t298 - t244 * t299 - t245 * t297 + t261 * t315 + t262 * t314 - t263 * t303 - t264 * t301 + t265 * t304 + t266 * t302) * t315) * t314; t335 * t303 + t334 * t301 + (m(5) * t200 + m(6) * t198 + m(7) * t388) * t300 + (-m(5) * t201 - m(6) * t199 - m(7) * t387) * t298 + (-t301 * t150 - t298 * t226 + t394) * t438; (m(5) * t81 + t352) * t330 + ((-m(5) * t88 - t350) * t329 + (m(5) * t87 + t351) * t327) * t328; (t26 * t65 + t35 * t77 + t36 * t78 + t42 * t66 + t46 * t84 + t47 * t85) * m(7) - (t21 / 0.2e1 + t20 / 0.2e1 + t9 / 0.2e1) * t437 + (t19 / 0.2e1 + t17 / 0.2e1 + t8 / 0.2e1) * t303 + (t18 / 0.2e1 + t16 / 0.2e1 + t7 / 0.2e1) * t301 + m(6) * (t119 * t75 + t120 * t76 + t45 * t86 + t58 * t93 + t59 * t94 + t64 * t92) + m(5) * (t113 * t138 + t114 * t139 + t121 * t81 + t128 * t87 + t129 * t88 + t131 * t91) + (t435 * t314 + t434 * t315 + t429 * t319 + t33) * t413 + (t433 * t314 + t432 * t315 + t428 * t319 + t34) * t412 + (t431 * t314 + t430 * t315 + t427 * t319 + t39) * t410 + t425 * t409 + t424 * t408 + t423 * t407; (t26 * t66 + t35 * t84 + t36 * t85) * t355 + t309 * t38 + (t119 * t58 + t120 * t59 + t45 * t92) * t357 + (t131 * t81 + t138 * t87 + t139 * t88) * t358 - (t423 + t439) * t437 + (t430 * t309 + t424) * t303 + (t431 * t309 + t425) * t301 + (t433 * t301 + t432 * t303 - t428 * t437 + t32) * t300 + (t435 * t301 + t434 * t303 - t429 * t437 + t31) * t298; t277 * t438; t372 * (t277 * t421 + (-t229 * t329 + t231 * t327) * t373); t349 * t305 + t348 * t281 + t347 * t279 + (m(6) * t86 + m(7) * t65) * t277 + (m(6) * t93 + m(7) * t77) * t231 + (m(6) * t94 + m(7) * t78) * t229; t352 * t305 + t351 * t281 + t350 * t279 + (m(6) * t92 + m(7) * t66) * t277 + (m(6) * t119 + m(7) * t84) * t231 + (m(6) * t120 + m(7) * t85) * t229; 0.4e1 * t372 * (t229 * t279 + t231 * t281 + t277 * t305); t44 * m(7); (t44 * t421 + (t327 * t48 - t329 * t49) * t373) * t419; (t122 * t46 + t123 * t47 + t42 * t95 + t44 * t65 + t48 * t77 + t49 * t78) * m(7) + t34 * t417 + t8 * t414 + t39 * t416 + t9 * t411 + t2 * t408 + t1 * t409 + t5 * t407 + t33 * t418 + t7 * t415; (t122 * t35 + t123 * t36 + t26 * t95 + t44 * t66 + t48 * t84 + t49 * t85) * m(7) + t30 * t412 + t303 * t2 / 0.2e1 + t31 * t418 + t3 * t415 + t32 * t417 + t4 * t414 + t29 * t413 + t301 * t1 / 0.2e1 + t38 * t416 + t6 * t411 + t37 * t410 - t437 * t5 / 0.2e1; (t122 * t231 + t123 * t229 + t277 * t95 + t279 * t49 + t281 * t48 + t305 * t44) * m(7); t232 * t30 + t282 * t2 + t230 * t29 + t280 * t1 + t278 * t37 + t306 * t5 + (t122 * t48 + t123 * t49 + t44 * t95) * t355;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
