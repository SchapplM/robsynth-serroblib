% Calculate time derivative of joint inertia matrix for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRPR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_inertiaDJ_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRPR1_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRPR1_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:45:27
% EndTime: 2019-03-08 18:45:56
% DurationCPUTime: 18.32s
% Computational Cost: add. (107690->932), mult. (296359->1243), div. (0->0), fcn. (367466->16), ass. (0->384)
t341 = sin(pkin(11));
t345 = cos(pkin(6));
t340 = sin(pkin(12));
t344 = cos(pkin(11));
t414 = t344 * t340;
t425 = cos(pkin(12));
t329 = t341 * t425 + t345 * t414;
t348 = sin(qJ(3));
t377 = t345 * t425;
t415 = t341 * t340;
t355 = -t344 * t377 + t415;
t426 = cos(pkin(7));
t352 = t355 * t426;
t342 = sin(pkin(6));
t424 = sin(pkin(7));
t378 = t342 * t424;
t429 = cos(qJ(3));
t312 = t329 * t429 + (-t344 * t378 - t352) * t348;
t379 = t342 * t426;
t323 = -t344 * t379 + t355 * t424;
t347 = sin(qJ(4));
t428 = cos(qJ(4));
t295 = t312 * t428 + t323 * t347;
t363 = t429 * t424;
t361 = t342 * t363;
t311 = t329 * t348 + t344 * t361 + t352 * t429;
t339 = sin(pkin(13));
t343 = cos(pkin(13));
t249 = -t295 * t339 + t311 * t343;
t421 = t311 * t339;
t250 = t295 * t343 + t421;
t359 = -t312 * t347 + t323 * t428;
t172 = Icges(6,5) * t250 + Icges(6,6) * t249 - Icges(6,3) * t359;
t213 = Icges(5,4) * t295 + Icges(5,2) * t359 + Icges(5,6) * t311;
t465 = -t172 + t213;
t330 = t344 * t425 - t345 * t415;
t354 = t341 * t377 + t414;
t351 = t354 * t426;
t314 = t330 * t429 + (t341 * t378 - t351) * t348;
t324 = t341 * t379 + t354 * t424;
t297 = t314 * t428 + t324 * t347;
t313 = t330 * t348 - t341 * t361 + t351 * t429;
t251 = -t297 * t339 + t313 * t343;
t420 = t313 * t339;
t252 = t297 * t343 + t420;
t358 = -t314 * t347 + t324 * t428;
t173 = Icges(6,5) * t252 + Icges(6,6) * t251 - Icges(6,3) * t358;
t214 = Icges(5,4) * t297 + Icges(5,2) * t358 + Icges(5,6) * t313;
t464 = -t173 + t214;
t362 = t426 * t425;
t322 = t345 * t424 * t348 + (t340 * t429 + t348 * t362) * t342;
t328 = t345 * t426 - t378 * t425;
t316 = t322 * t428 + t328 * t347;
t460 = t342 * (-t340 * t348 + t429 * t362) + t345 * t363;
t292 = -t316 * t339 - t343 * t460;
t417 = t460 * t339;
t293 = t316 * t343 - t417;
t357 = -t322 * t347 + t328 * t428;
t208 = Icges(6,5) * t293 + Icges(6,6) * t292 - Icges(6,3) * t357;
t263 = Icges(5,4) * t316 + Icges(5,2) * t357 - Icges(5,6) * t460;
t463 = -t208 + t263;
t174 = Icges(6,4) * t250 + Icges(6,2) * t249 - Icges(6,6) * t359;
t176 = Icges(6,1) * t250 + Icges(6,4) * t249 - Icges(6,5) * t359;
t211 = Icges(5,5) * t295 + Icges(5,6) * t359 + Icges(5,3) * t311;
t215 = Icges(5,1) * t295 + Icges(5,4) * t359 + Icges(5,5) * t311;
t458 = t174 * t249 + t176 * t250 + t211 * t311 + t215 * t295 + t359 * t465;
t175 = Icges(6,4) * t252 + Icges(6,2) * t251 - Icges(6,6) * t358;
t177 = Icges(6,1) * t252 + Icges(6,4) * t251 - Icges(6,5) * t358;
t212 = Icges(5,5) * t297 + Icges(5,6) * t358 + Icges(5,3) * t313;
t216 = Icges(5,1) * t297 + Icges(5,4) * t358 + Icges(5,5) * t313;
t457 = t175 * t249 + t177 * t250 + t212 * t311 + t216 * t295 + t359 * t464;
t456 = t174 * t251 + t176 * t252 + t211 * t313 + t215 * t297 + t358 * t465;
t455 = t175 * t251 + t177 * t252 + t212 * t313 + t216 * t297 + t358 * t464;
t454 = t174 * t292 + t176 * t293 - t211 * t460 + t215 * t316 + t357 * t465;
t453 = t175 * t292 + t177 * t293 - t212 * t460 + t216 * t316 + t357 * t464;
t318 = t322 * qJD(3);
t209 = Icges(6,4) * t293 + Icges(6,2) * t292 - Icges(6,6) * t357;
t210 = Icges(6,1) * t293 + Icges(6,4) * t292 - Icges(6,5) * t357;
t262 = Icges(5,5) * t316 + Icges(5,6) * t357 - Icges(5,3) * t460;
t264 = Icges(5,1) * t316 + Icges(5,4) * t357 - Icges(5,5) * t460;
t450 = t209 * t292 + t210 * t293 - t262 * t460 + t264 * t316 + t357 * t463;
t462 = t318 * t450;
t452 = t209 * t249 + t210 * t250 + t262 * t311 + t264 * t295 + t359 * t463;
t451 = t209 * t251 + t210 * t252 + t262 * t313 + t264 * t297 + t358 * t463;
t442 = m(7) / 0.2e1;
t387 = m(6) / 0.2e1 + t442;
t461 = 0.2e1 * t387;
t459 = 2 * m(4);
t307 = t311 * qJD(3);
t245 = qJD(4) * t295 - t307 * t347;
t246 = qJD(4) * t359 - t307 * t428;
t170 = t246 * pkin(4) + t245 * qJ(5) - qJD(5) * t359;
t237 = t295 * pkin(4) - qJ(5) * t359;
t310 = t314 * qJD(3);
t407 = t313 * t170 + t310 * t237;
t317 = t460 * qJD(3);
t290 = t316 * qJD(4) + t317 * t347;
t338 = pkin(13) + qJ(6);
t336 = sin(t338);
t337 = cos(t338);
t242 = t295 * t337 + t311 * t336;
t308 = t312 * qJD(3);
t180 = -qJD(6) * t242 - t246 * t336 + t308 * t337;
t241 = -t295 * t336 + t311 * t337;
t181 = qJD(6) * t241 + t246 * t337 + t308 * t336;
t112 = Icges(7,5) * t181 + Icges(7,6) * t180 + Icges(7,3) * t245;
t114 = Icges(7,4) * t181 + Icges(7,2) * t180 + Icges(7,6) * t245;
t116 = Icges(7,1) * t181 + Icges(7,4) * t180 + Icges(7,5) * t245;
t162 = Icges(7,5) * t242 + Icges(7,6) * t241 - Icges(7,3) * t359;
t164 = Icges(7,4) * t242 + Icges(7,2) * t241 - Icges(7,6) * t359;
t166 = Icges(7,1) * t242 + Icges(7,4) * t241 - Icges(7,5) * t359;
t25 = -t112 * t359 + t114 * t241 + t116 * t242 + t162 * t245 + t164 * t180 + t166 * t181;
t244 = t297 * t337 + t313 * t336;
t309 = t313 * qJD(3);
t248 = qJD(4) * t358 - t309 * t428;
t182 = -qJD(6) * t244 - t248 * t336 + t310 * t337;
t243 = -t297 * t336 + t313 * t337;
t183 = qJD(6) * t243 + t248 * t337 + t310 * t336;
t247 = qJD(4) * t297 - t309 * t347;
t113 = Icges(7,5) * t183 + Icges(7,6) * t182 + Icges(7,3) * t247;
t115 = Icges(7,4) * t183 + Icges(7,2) * t182 + Icges(7,6) * t247;
t117 = Icges(7,1) * t183 + Icges(7,4) * t182 + Icges(7,5) * t247;
t163 = Icges(7,5) * t244 + Icges(7,6) * t243 - Icges(7,3) * t358;
t165 = Icges(7,4) * t244 + Icges(7,2) * t243 - Icges(7,6) * t358;
t167 = Icges(7,1) * t244 + Icges(7,4) * t243 - Icges(7,5) * t358;
t26 = -t113 * t359 + t115 * t241 + t117 * t242 + t163 * t245 + t165 * t180 + t167 * t181;
t287 = t316 * t337 - t336 * t460;
t291 = qJD(4) * t357 + t317 * t428;
t223 = -qJD(6) * t287 - t291 * t336 + t318 * t337;
t286 = -t316 * t336 - t337 * t460;
t224 = qJD(6) * t286 + t291 * t337 + t318 * t336;
t143 = Icges(7,5) * t224 + Icges(7,6) * t223 + Icges(7,3) * t290;
t144 = Icges(7,4) * t224 + Icges(7,2) * t223 + Icges(7,6) * t290;
t145 = Icges(7,1) * t224 + Icges(7,4) * t223 + Icges(7,5) * t290;
t200 = Icges(7,5) * t287 + Icges(7,6) * t286 - Icges(7,3) * t357;
t201 = Icges(7,4) * t287 + Icges(7,2) * t286 - Icges(7,6) * t357;
t202 = Icges(7,1) * t287 + Icges(7,4) * t286 - Icges(7,5) * t357;
t45 = -t143 * t359 + t144 * t241 + t145 * t242 + t180 * t201 + t181 * t202 + t200 * t245;
t78 = -t162 * t359 + t164 * t241 + t166 * t242;
t79 = -t163 * t359 + t165 * t241 + t167 * t242;
t98 = -t200 * t359 + t201 * t241 + t202 * t242;
t3 = t25 * t311 + t26 * t313 + t308 * t78 + t310 * t79 + t318 * t98 - t45 * t460;
t229 = -t246 * t339 + t308 * t343;
t423 = t308 * t339;
t230 = t246 * t343 + t423;
t134 = Icges(6,5) * t230 + Icges(6,6) * t229 + Icges(6,3) * t245;
t136 = Icges(6,4) * t230 + Icges(6,2) * t229 + Icges(6,6) * t245;
t138 = Icges(6,1) * t230 + Icges(6,4) * t229 + Icges(6,5) * t245;
t41 = -t134 * t359 + t136 * t249 + t138 * t250 + t172 * t245 + t174 * t229 + t176 * t230;
t231 = -t248 * t339 + t310 * t343;
t422 = t310 * t339;
t232 = t248 * t343 + t422;
t135 = Icges(6,5) * t232 + Icges(6,6) * t231 + Icges(6,3) * t247;
t137 = Icges(6,4) * t232 + Icges(6,2) * t231 + Icges(6,6) * t247;
t139 = Icges(6,1) * t232 + Icges(6,4) * t231 + Icges(6,5) * t247;
t42 = -t135 * t359 + t137 * t249 + t139 * t250 + t173 * t245 + t175 * t229 + t177 * t230;
t269 = -t291 * t339 + t318 * t343;
t418 = t318 * t339;
t270 = t291 * t343 + t418;
t192 = Icges(6,5) * t270 + Icges(6,6) * t269 + Icges(6,3) * t290;
t193 = Icges(6,4) * t270 + Icges(6,2) * t269 + Icges(6,6) * t290;
t194 = Icges(6,1) * t270 + Icges(6,4) * t269 + Icges(6,5) * t290;
t56 = -t192 * t359 + t193 * t249 + t194 * t250 + t208 * t245 + t209 * t229 + t210 * t230;
t184 = Icges(5,5) * t246 - Icges(5,6) * t245 + Icges(5,3) * t308;
t186 = Icges(5,4) * t246 - Icges(5,2) * t245 + Icges(5,6) * t308;
t188 = Icges(5,1) * t246 - Icges(5,4) * t245 + Icges(5,5) * t308;
t61 = t184 * t311 + t186 * t359 + t188 * t295 + t211 * t308 - t213 * t245 + t215 * t246;
t185 = Icges(5,5) * t248 - Icges(5,6) * t247 + Icges(5,3) * t310;
t187 = Icges(5,4) * t248 - Icges(5,2) * t247 + Icges(5,6) * t310;
t189 = Icges(5,1) * t248 - Icges(5,4) * t247 + Icges(5,5) * t310;
t62 = t185 * t311 + t187 * t359 + t189 * t295 + t212 * t308 - t214 * t245 + t216 * t246;
t225 = Icges(5,5) * t291 - Icges(5,6) * t290 + Icges(5,3) * t318;
t226 = Icges(5,4) * t291 - Icges(5,2) * t290 + Icges(5,6) * t318;
t227 = Icges(5,1) * t291 - Icges(5,4) * t290 + Icges(5,5) * t318;
t76 = t225 * t311 + t226 * t359 + t227 * t295 - t245 * t263 + t246 * t264 + t262 * t308;
t448 = t3 - (t56 + t76) * t460 + t452 * t318 + (t42 + t62) * t313 + (t41 + t61) * t311 + t457 * t310 + t458 * t308;
t27 = -t112 * t358 + t114 * t243 + t116 * t244 + t162 * t247 + t164 * t182 + t166 * t183;
t28 = -t113 * t358 + t115 * t243 + t117 * t244 + t163 * t247 + t165 * t182 + t167 * t183;
t46 = -t143 * t358 + t144 * t243 + t145 * t244 + t182 * t201 + t183 * t202 + t200 * t247;
t80 = -t162 * t358 + t164 * t243 + t166 * t244;
t81 = -t163 * t358 + t165 * t243 + t167 * t244;
t99 = -t200 * t358 + t201 * t243 + t202 * t244;
t4 = t27 * t311 + t28 * t313 + t308 * t80 + t310 * t81 + t318 * t99 - t46 * t460;
t43 = -t134 * t358 + t136 * t251 + t138 * t252 + t172 * t247 + t174 * t231 + t176 * t232;
t44 = -t135 * t358 + t137 * t251 + t139 * t252 + t173 * t247 + t175 * t231 + t177 * t232;
t57 = -t192 * t358 + t193 * t251 + t194 * t252 + t208 * t247 + t209 * t231 + t210 * t232;
t63 = t184 * t313 + t186 * t358 + t188 * t297 + t211 * t310 - t213 * t247 + t215 * t248;
t64 = t185 * t313 + t187 * t358 + t189 * t297 + t212 * t310 - t214 * t247 + t216 * t248;
t77 = t225 * t313 + t226 * t358 + t227 * t297 - t247 * t263 + t248 * t264 + t262 * t310;
t447 = t4 - (t57 + t77) * t460 + t451 * t318 + (t44 + t64) * t313 + (t43 + t63) * t311 + t455 * t310 + t456 * t308;
t50 = -t134 * t357 + t136 * t292 + t138 * t293 + t172 * t290 + t174 * t269 + t176 * t270;
t51 = -t135 * t357 + t137 * t292 + t139 * t293 + t173 * t290 + t175 * t269 + t177 * t270;
t106 = -t200 * t357 + t201 * t286 + t202 * t287;
t29 = -t112 * t357 + t114 * t286 + t116 * t287 + t162 * t290 + t164 * t223 + t166 * t224;
t30 = -t113 * t357 + t115 * t286 + t117 * t287 + t163 * t290 + t165 * t223 + t167 * t224;
t52 = -t143 * t357 + t144 * t286 + t145 * t287 + t200 * t290 + t201 * t223 + t202 * t224;
t87 = -t162 * t357 + t164 * t286 + t166 * t287;
t88 = -t163 * t357 + t165 * t286 + t167 * t287;
t6 = t106 * t318 + t29 * t311 + t30 * t313 + t308 * t87 + t310 * t88 - t460 * t52;
t65 = -t192 * t357 + t193 * t292 + t194 * t293 + t208 * t290 + t209 * t269 + t210 * t270;
t68 = -t184 * t460 + t186 * t357 + t188 * t316 + t211 * t318 - t213 * t290 + t215 * t291;
t69 = -t185 * t460 + t187 * t357 + t189 * t316 + t212 * t318 - t214 * t290 + t216 * t291;
t89 = -t225 * t460 + t226 * t357 + t227 * t316 + t262 * t318 - t263 * t290 + t264 * t291;
t446 = t6 - (t65 + t89) * t460 + t462 + (t51 + t69) * t313 + (t50 + t68) * t311 + t453 * t310 + t454 * t308;
t445 = -0.2e1 * t323;
t444 = 0.2e1 * t345;
t441 = t245 / 0.2e1;
t440 = t247 / 0.2e1;
t439 = t290 / 0.2e1;
t438 = -t359 / 0.2e1;
t437 = -t358 / 0.2e1;
t436 = t308 / 0.2e1;
t435 = t310 / 0.2e1;
t434 = -t357 / 0.2e1;
t433 = t318 / 0.2e1;
t432 = t323 / 0.2e1;
t431 = t324 / 0.2e1;
t430 = t328 / 0.2e1;
t427 = pkin(5) * t343;
t118 = rSges(7,1) * t181 + rSges(7,2) * t180 + rSges(7,3) * t245;
t412 = pkin(5) * t423 + pkin(10) * t245 + t246 * t427 + t118;
t119 = rSges(7,1) * t183 + rSges(7,2) * t182 + rSges(7,3) * t247;
t411 = pkin(5) * t422 + pkin(10) * t247 + t248 * t427 + t119;
t140 = rSges(6,1) * t230 + rSges(6,2) * t229 + rSges(6,3) * t245;
t410 = -t140 - t170;
t141 = rSges(6,1) * t232 + rSges(6,2) * t231 + rSges(6,3) * t247;
t171 = t248 * pkin(4) + t247 * qJ(5) - qJD(5) * t358;
t409 = -t141 - t171;
t146 = rSges(7,1) * t224 + rSges(7,2) * t223 + rSges(7,3) * t290;
t408 = pkin(5) * t418 + pkin(10) * t290 + t291 * t427 + t146;
t238 = t297 * pkin(4) - qJ(5) * t358;
t406 = -t171 * t460 + t318 * t238;
t156 = t324 * t170;
t280 = -pkin(3) * t307 + pkin(9) * t308;
t255 = t324 * t280;
t405 = t156 + t255;
t168 = rSges(7,1) * t242 + rSges(7,2) * t241 - rSges(7,3) * t359;
t404 = pkin(5) * t421 - pkin(10) * t359 + t295 * t427 + t168;
t169 = rSges(7,1) * t244 + rSges(7,2) * t243 - rSges(7,3) * t358;
t403 = pkin(5) * t420 - pkin(10) * t358 + t297 * t427 + t169;
t281 = -pkin(3) * t309 + pkin(9) * t310;
t268 = t328 * t281;
t402 = t328 * t171 + t268;
t178 = rSges(6,1) * t250 + rSges(6,2) * t249 - rSges(6,3) * t359;
t401 = -t178 - t237;
t179 = rSges(6,1) * t252 + rSges(6,2) * t251 - rSges(6,3) * t358;
t400 = -t179 - t238;
t195 = rSges(6,1) * t270 + rSges(6,2) * t269 + rSges(6,3) * t290;
t204 = t291 * pkin(4) + t290 * qJ(5) - qJD(5) * t357;
t399 = -t195 - t204;
t284 = t316 * pkin(4) - qJ(5) * t357;
t398 = t311 * t204 + t308 * t284;
t305 = pkin(3) * t317 + pkin(9) * t318;
t285 = t323 * t305;
t397 = t323 * t204 + t285;
t203 = rSges(7,1) * t287 + rSges(7,2) * t286 - rSges(7,3) * t357;
t396 = -pkin(5) * t417 - pkin(10) * t357 + t316 * t427 + t203;
t217 = rSges(6,1) * t293 + rSges(6,2) * t292 - rSges(6,3) * t357;
t395 = -t217 - t284;
t282 = pkin(3) * t312 + pkin(9) * t311;
t266 = t324 * t282;
t394 = t324 * t237 + t266;
t283 = pkin(3) * t314 + pkin(9) * t313;
t271 = t328 * t283;
t393 = t328 * t238 + t271;
t392 = t281 * t445 + 0.2e1 * t255;
t306 = pkin(3) * t322 - pkin(9) * t460;
t288 = t323 * t306;
t391 = t323 * t284 + t288;
t388 = 0.2e1 * t342;
t386 = -t170 - t412;
t385 = -t171 - t411;
t384 = -t204 - t408;
t383 = -t237 - t404;
t382 = -t238 - t403;
t381 = -t284 - t396;
t376 = 2 * m(5);
t375 = 0.2e1 * m(6);
t373 = 0.2e1 * m(7);
t19 = t308 * t382 + t310 * t404 + t311 * t385 + t313 * t412 + t407;
t49 = t313 * t140 + t310 * t178 + t308 * t400 + t311 * t409 + t407;
t370 = m(6) * t49 + m(7) * t19;
t23 = t308 * t396 + t311 * t408 + t318 * t383 - t386 * t460 + t398;
t54 = t311 * t195 + t308 * t217 + t318 * t401 - t410 * t460 + t398;
t369 = m(6) * t54 + m(7) * t23;
t24 = t310 * t381 + t313 * t384 + t318 * t403 - t411 * t460 + t406;
t55 = -t141 * t460 + t318 * t179 + t310 * t395 + t313 * t399 + t406;
t368 = m(6) * t55 + m(7) * t24;
t35 = t412 * t324 + (-t281 + t385) * t323 + t405;
t66 = t324 * t140 + (-t281 + t409) * t323 + t405;
t367 = m(6) * t66 + m(7) * t35;
t47 = t408 * t323 + (-t280 + t386) * t328 + t397;
t72 = t323 * t195 + (-t280 + t410) * t328 + t397;
t366 = m(6) * t72 + m(7) * t47;
t48 = t411 * t328 + (-t305 + t384) * t324 + t402;
t73 = t328 * t141 + (-t305 + t399) * t324 + t402;
t365 = m(6) * t73 + m(7) * t48;
t53 = -t118 * t358 + t119 * t359 + t168 * t247 - t169 * t245;
t190 = rSges(5,1) * t246 - rSges(5,2) * t245 + rSges(5,3) * t308;
t350 = m(5) * t190 + m(6) * t140 + m(7) * t412;
t191 = rSges(5,1) * t248 - rSges(5,2) * t247 + rSges(5,3) * t310;
t349 = -m(5) * t191 - m(6) * t141 - m(7) * t411;
t304 = rSges(4,1) * t317 - rSges(4,2) * t318;
t303 = Icges(4,1) * t317 - Icges(4,4) * t318;
t302 = Icges(4,4) * t317 - Icges(4,2) * t318;
t301 = Icges(4,5) * t317 - Icges(4,6) * t318;
t300 = rSges(4,1) * t322 + rSges(4,2) * t460 + rSges(4,3) * t328;
t299 = Icges(4,1) * t322 + Icges(4,4) * t460 + Icges(4,5) * t328;
t298 = Icges(4,4) * t322 + Icges(4,2) * t460 + Icges(4,6) * t328;
t279 = -rSges(4,1) * t309 - rSges(4,2) * t310;
t278 = -rSges(4,1) * t307 - rSges(4,2) * t308;
t277 = -Icges(4,1) * t309 - Icges(4,4) * t310;
t276 = -Icges(4,1) * t307 - Icges(4,4) * t308;
t275 = -Icges(4,4) * t309 - Icges(4,2) * t310;
t274 = -Icges(4,4) * t307 - Icges(4,2) * t308;
t273 = -Icges(4,5) * t309 - Icges(4,6) * t310;
t272 = -Icges(4,5) * t307 - Icges(4,6) * t308;
t265 = rSges(5,1) * t316 + rSges(5,2) * t357 - rSges(5,3) * t460;
t261 = rSges(4,1) * t314 - rSges(4,2) * t313 + rSges(4,3) * t324;
t260 = rSges(4,1) * t312 - rSges(4,2) * t311 + rSges(4,3) * t323;
t259 = Icges(4,1) * t314 - Icges(4,4) * t313 + Icges(4,5) * t324;
t258 = Icges(4,1) * t312 - Icges(4,4) * t311 + Icges(4,5) * t323;
t257 = Icges(4,4) * t314 - Icges(4,2) * t313 + Icges(4,6) * t324;
t256 = Icges(4,4) * t312 - Icges(4,2) * t311 + Icges(4,6) * t323;
t240 = t311 * t284;
t234 = t460 * t238;
t228 = rSges(5,1) * t291 - rSges(5,2) * t290 + rSges(5,3) * t318;
t222 = t313 * t237;
t221 = t279 * t328 - t304 * t324;
t220 = -t278 * t328 + t304 * t323;
t219 = rSges(5,1) * t297 + rSges(5,2) * t358 + rSges(5,3) * t313;
t218 = rSges(5,1) * t295 + rSges(5,2) * t359 + rSges(5,3) * t311;
t196 = t278 * t324 - t279 * t323;
t149 = -t219 * t460 - t265 * t313;
t148 = t218 * t460 + t265 * t311;
t142 = t218 * t313 - t219 * t311;
t132 = t328 * t219 + t271 + (-t265 - t306) * t324;
t131 = t323 * t265 + t288 + (-t218 - t282) * t328;
t126 = -t169 * t357 + t203 * t358;
t125 = t168 * t357 - t203 * t359;
t124 = t324 * t218 + t266 + (-t219 - t283) * t323;
t121 = t328 * t191 + t268 + (-t228 - t305) * t324;
t120 = t323 * t228 + t285 + (-t190 - t280) * t328;
t105 = -t168 * t358 + t169 * t359;
t104 = t324 * t190 + t255 + (-t191 - t281) * t323;
t103 = -t179 * t460 + t313 * t395 - t234;
t102 = t311 * t217 - t401 * t460 + t240;
t97 = -t191 * t460 + t219 * t318 - t228 * t313 - t265 * t310;
t96 = t190 * t460 - t218 * t318 + t228 * t311 + t265 * t308;
t95 = t328 * t179 + (-t306 + t395) * t324 + t393;
t94 = t323 * t217 + (-t282 + t401) * t328 + t391;
t93 = t313 * t178 + t311 * t400 + t222;
t90 = t190 * t313 - t191 * t311 + t218 * t310 - t219 * t308;
t86 = t324 * t178 + (-t283 + t400) * t323 + t394;
t75 = t313 * t381 - t403 * t460 - t234;
t74 = t311 * t396 - t383 * t460 + t240;
t71 = t403 * t328 + (-t306 + t381) * t324 + t393;
t70 = t396 * t323 + (-t282 + t383) * t328 + t391;
t67 = t311 * t382 + t313 * t404 + t222;
t60 = t404 * t324 + (-t283 + t382) * t323 + t394;
t59 = -t119 * t357 + t146 * t358 + t169 * t290 - t203 * t247;
t58 = t118 * t357 - t146 * t359 - t168 * t290 + t203 * t245;
t40 = t106 * t328 + t323 * t87 + t324 * t88;
t39 = -t106 * t460 + t311 * t87 + t313 * t88;
t38 = -t106 * t357 - t358 * t88 - t359 * t87;
t37 = t323 * t80 + t324 * t81 + t328 * t99;
t36 = t323 * t78 + t324 * t79 + t328 * t98;
t34 = t311 * t80 + t313 * t81 - t460 * t99;
t33 = t311 * t78 + t313 * t79 - t460 * t98;
t32 = -t357 * t99 - t358 * t81 - t359 * t80;
t31 = -t357 * t98 - t358 * t79 - t359 * t78;
t22 = t323 * t68 + t324 * t69 + t328 * t89;
t21 = t323 * t63 + t324 * t64 + t328 * t77;
t20 = t323 * t61 + t324 * t62 + t328 * t76;
t15 = t323 * t50 + t324 * t51 + t328 * t65;
t14 = t323 * t43 + t324 * t44 + t328 * t57;
t13 = t323 * t41 + t324 * t42 + t328 * t56;
t12 = t29 * t323 + t30 * t324 + t328 * t52;
t10 = t27 * t323 + t28 * t324 + t328 * t46;
t9 = t25 * t323 + t26 * t324 + t328 * t45;
t5 = t106 * t290 + t245 * t87 + t247 * t88 - t29 * t359 - t30 * t358 - t357 * t52;
t2 = t245 * t80 + t247 * t81 - t27 * t359 - t28 * t358 + t290 * t99 - t357 * t46;
t1 = t245 * t78 + t247 * t79 - t25 * t359 - t26 * t358 + t290 * t98 - t357 * t45;
t7 = [0; 0; 0; m(5) * t392 / 0.2e1 + (m(4) * t278 + t350) * t324 + (-m(4) * t279 + t349) * t323 + t387 * (t171 * t445 + 0.2e1 * t156 + t392); (m(4) * t196 + m(5) * t104 + t367) * t345 + ((-m(4) * t221 - m(5) * t121 - t365) * t344 + (m(4) * t220 + m(5) * t120 + t366) * t341) * t342; (t35 * t60 + t47 * t70 + t48 * t71) * t373 + (t66 * t86 + t72 * t94 + t73 * t95) * t375 + (t104 * t124 + t120 * t131 + t121 * t132) * t376 + (t12 + t15 + t22 + (-t220 * t260 + t221 * t261) * t459 + (-t298 * t318 + t299 * t317 + t301 * t328 + t302 * t460 + t303 * t322) * t328) * t328 + (t10 + t21 + t14 + (t196 * t260 - t221 * t300) * t459 + (-t257 * t310 - t259 * t309 + t273 * t324 - t275 * t313 + t277 * t314) * t324 + (-t257 * t318 + t259 * t317 + t273 * t328 + t275 * t460 + t277 * t322 - t298 * t310 - t299 * t309 + t301 * t324 - t302 * t313 + t303 * t314) * t328) * t324 + (t9 + t13 + t20 + (-t196 * t261 + t220 * t300) * t459 + (-t256 * t308 - t258 * t307 + t272 * t323 - t274 * t311 + t276 * t312) * t323 + (-t256 * t318 + t258 * t317 + t272 * t328 + t274 * t460 + t276 * t322 - t298 * t308 - t299 * t307 + t301 * t323 - t302 * t311 + t303 * t312) * t328 + (-t256 * t310 - t257 * t308 - t258 * t309 - t259 * t307 + t272 * t324 + t273 * t323 - t274 * t313 - t275 * t311 + t276 * t314 + t277 * t312) * t324) * t323; t350 * t313 + t349 * t311 + (m(5) * t218 + m(6) * t178 + m(7) * t404) * t310 + (-m(5) * t219 - m(6) * t179 - m(7) * t403) * t308 + (-t311 * t171 - t308 * t238 + t407) * t461; (m(5) * t90 + t370) * t345 + ((-m(5) * t97 - t368) * t344 + (m(5) * t96 + t369) * t341) * t342; (t19 * t60 + t23 * t70 + t24 * t71 + t35 * t67 + t47 * t74 + t48 * t75) * m(7) - (t22 / 0.2e1 + t15 / 0.2e1 + t12 / 0.2e1) * t460 + (t21 / 0.2e1 + t14 / 0.2e1 + t10 / 0.2e1) * t313 + (t20 / 0.2e1 + t13 / 0.2e1 + t9 / 0.2e1) * t311 + m(5) * (t104 * t142 + t120 * t148 + t121 * t149 + t124 * t90 + t131 * t96 + t132 * t97) + m(6) * (t102 * t72 + t103 * t73 + t49 * t86 + t54 * t94 + t55 * t95 + t66 * t93) + (t323 * t458 + t324 * t457 + t328 * t452 + t36) * t436 + (t323 * t456 + t324 * t455 + t328 * t451 + t37) * t435 + (t323 * t454 + t324 * t453 + t328 * t450 + t40) * t433 + t448 * t432 + t447 * t431 + t446 * t430; (t19 * t67 + t23 * t74 + t24 * t75) * t373 + t318 * t39 + (t102 * t54 + t103 * t55 + t49 * t93) * t375 + (t142 * t90 + t148 * t96 + t149 * t97) * t376 - (t446 + t462) * t460 + (t318 * t453 + t447) * t313 + (t318 * t454 + t448) * t311 + (t311 * t456 + t313 * t455 - t451 * t460 + t34) * t310 + (t311 * t458 + t313 * t457 - t452 * t460 + t33) * t308; t290 * t461; t387 * (t290 * t444 + (-t245 * t344 + t247 * t341) * t388); -t367 * t357 - t366 * t358 - t365 * t359 + (m(6) * t86 + m(7) * t60) * t290 + (m(6) * t94 + m(7) * t70) * t247 + (m(6) * t95 + m(7) * t71) * t245; -t370 * t357 - t369 * t358 - t368 * t359 + (m(6) * t93 + m(7) * t67) * t290 + (m(6) * t102 + m(7) * t74) * t247 + (m(6) * t103 + m(7) * t75) * t245; 0.4e1 * t387 * (-t245 * t359 - t247 * t358 - t290 * t357); t53 * m(7); (t53 * t444 + (t341 * t58 - t344 * t59) * t388) * t442; (t105 * t35 + t125 * t47 + t126 * t48 + t53 * t60 + t58 * t70 + t59 * t71) * m(7) + t2 * t431 + t1 * t432 + t40 * t439 + t12 * t434 + t5 * t430 + t36 * t441 + t9 * t438 + t37 * t440 + t10 * t437; (t105 * t19 + t125 * t23 + t126 * t24 + t53 * t67 + t58 * t74 + t59 * t75) * m(7) + t39 * t439 + t6 * t434 + t31 * t436 + t311 * t1 / 0.2e1 + t38 * t433 - t460 * t5 / 0.2e1 + t32 * t435 + t313 * t2 / 0.2e1 + t33 * t441 + t3 * t438 + t34 * t440 + t4 * t437; (t105 * t290 + t125 * t247 + t126 * t245 - t357 * t53 - t358 * t58 - t359 * t59) * m(7); t245 * t31 - t359 * t1 + t290 * t38 - t357 * t5 + t247 * t32 - t358 * t2 + (t105 * t53 + t125 * t58 + t126 * t59) * t373;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t7(1) t7(2) t7(4) t7(7) t7(11) t7(16); t7(2) t7(3) t7(5) t7(8) t7(12) t7(17); t7(4) t7(5) t7(6) t7(9) t7(13) t7(18); t7(7) t7(8) t7(9) t7(10) t7(14) t7(19); t7(11) t7(12) t7(13) t7(14) t7(15) t7(20); t7(16) t7(17) t7(18) t7(19) t7(20) t7(21);];
Mq  = res;
