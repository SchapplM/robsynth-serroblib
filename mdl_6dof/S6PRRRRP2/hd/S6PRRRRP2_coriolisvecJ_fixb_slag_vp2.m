% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRRRRP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:28:15
% EndTime: 2018-11-23 15:28:24
% DurationCPUTime: 9.46s
% Computational Cost: add. (8028->555), mult. (20187->732), div. (0->0), fcn. (14576->10), ass. (0->265)
t414 = Ifges(6,1) + Ifges(7,1);
t403 = Ifges(6,5) + Ifges(7,4);
t413 = Ifges(6,6) - Ifges(7,6);
t238 = sin(qJ(4));
t239 = sin(qJ(3));
t242 = cos(qJ(4));
t243 = cos(qJ(3));
t207 = t238 * t239 - t242 * t243;
t309 = qJD(3) + qJD(4);
t170 = t309 * t207;
t208 = t238 * t243 + t239 * t242;
t171 = t309 * t208;
t240 = sin(qJ(2));
t235 = sin(pkin(6));
t319 = qJD(1) * t235;
t298 = t240 * t319;
t313 = qJD(3) * t239;
t308 = pkin(3) * t313;
t415 = pkin(4) * t171 + pkin(10) * t170 - t298 + t308;
t376 = -pkin(9) - pkin(8);
t299 = qJD(3) * t376;
t211 = t239 * t299;
t212 = t243 * t299;
t244 = cos(qJ(2));
t297 = t244 * t319;
t222 = t376 * t239;
t223 = t376 * t243;
t384 = t242 * t222 + t223 * t238;
t405 = qJD(4) * t384 + t207 * t297 + t211 * t242 + t212 * t238;
t237 = sin(qJ(5));
t241 = cos(qJ(5));
t232 = -pkin(3) * t243 - pkin(2);
t163 = pkin(4) * t207 - pkin(10) * t208 + t232;
t184 = t222 * t238 - t223 * t242;
t385 = t237 * t163 + t241 * t184;
t410 = -qJD(5) * t385 - t237 * t405 + t241 * t415;
t311 = qJD(5) * t241;
t312 = qJD(5) * t237;
t408 = t163 * t311 - t184 * t312 + t237 * t415 + t241 * t405;
t397 = -t237 * t413 + t241 * t403;
t343 = Ifges(7,5) * t237;
t345 = Ifges(6,4) * t237;
t394 = t241 * t414 + t343 - t345;
t412 = Ifges(5,6) * t309 / 0.2e1;
t411 = -pkin(5) * t171 - t410;
t409 = qJ(6) * t171 + qJD(6) * t207 + t408;
t110 = qJD(4) * t184 + t211 * t238 - t242 * t212;
t179 = t208 * t297;
t407 = t110 - t179;
t202 = t207 * qJD(2);
t203 = t208 * qJD(2);
t161 = mrSges(5,1) * t202 + mrSges(5,2) * t203;
t194 = qJD(2) * t232 - t297;
t406 = -m(5) * t194 - t161;
t347 = Ifges(5,4) * t203;
t404 = -Ifges(5,2) * t202 / 0.2e1 + t347 / 0.2e1;
t159 = t171 * qJD(2);
t158 = t170 * qJD(2);
t181 = t237 * t203 - t241 * t309;
t98 = -qJD(5) * t181 - t241 * t158;
t182 = t241 * t203 + t237 * t309;
t99 = qJD(5) * t182 - t237 * t158;
t402 = (-Ifges(6,4) + Ifges(7,5)) * t99 + t414 * t98 + t403 * t159;
t178 = Ifges(6,4) * t181;
t196 = qJD(5) + t202;
t339 = t181 * Ifges(7,5);
t389 = t414 * t182 + t403 * t196 - t178 + t339;
t401 = Ifges(5,5) * t309;
t399 = t194 * mrSges(5,2);
t282 = mrSges(7,1) * t237 - mrSges(7,3) * t241;
t284 = mrSges(6,1) * t237 + mrSges(6,2) * t241;
t214 = qJD(2) * pkin(8) + t298;
t288 = pkin(9) * qJD(2) + t214;
t236 = cos(pkin(6));
t318 = qJD(1) * t239;
t294 = t236 * t318;
t176 = t243 * t288 + t294;
t165 = t238 * t176;
t327 = t236 * t243;
t226 = qJD(1) * t327;
t265 = t288 * t239;
t175 = t226 - t265;
t167 = qJD(3) * pkin(3) + t175;
t101 = t242 * t167 - t165;
t89 = -pkin(4) * t309 - t101;
t40 = t181 * pkin(5) - t182 * qJ(6) + t89;
t398 = t40 * t282 + t89 * t284;
t396 = t403 * t237 + t413 * t241;
t342 = Ifges(7,5) * t241;
t344 = Ifges(6,4) * t241;
t395 = t414 * t237 - t342 + t344;
t121 = pkin(4) * t202 - pkin(10) * t203 + t194;
t317 = qJD(2) * t235;
t290 = qJD(1) * t317;
t287 = t244 * t290;
t320 = qJD(3) * t226 + t243 * t287;
t129 = -qJD(3) * t265 + t320;
t383 = -t176 * qJD(3) - t239 * t287;
t28 = qJD(4) * t101 + t242 * t129 + t238 * t383;
t197 = qJD(2) * t308 + t240 * t290;
t62 = pkin(4) * t159 + pkin(10) * t158 + t197;
t166 = t242 * t176;
t102 = t238 * t167 + t166;
t90 = pkin(10) * t309 + t102;
t6 = t121 * t311 + t237 * t62 + t241 * t28 - t312 * t90;
t37 = t121 * t237 + t241 * t90;
t7 = -qJD(5) * t37 - t237 * t28 + t241 * t62;
t393 = -t7 * t237 + t241 * t6;
t2 = qJ(6) * t159 + qJD(6) * t196 + t6;
t36 = t121 * t241 - t237 * t90;
t386 = qJD(6) - t36;
t31 = -pkin(5) * t196 + t386;
t4 = -pkin(5) * t159 - t7;
t392 = t2 * t241 + t237 * t4 + t31 * t311;
t271 = Ifges(7,3) * t237 + t342;
t277 = -Ifges(6,2) * t237 + t344;
t364 = t237 / 0.2e1;
t365 = -t237 / 0.2e1;
t370 = t182 / 0.2e1;
t372 = t181 / 0.2e1;
t373 = -t181 / 0.2e1;
t177 = Ifges(7,5) * t182;
t82 = Ifges(7,6) * t196 + Ifges(7,3) * t181 + t177;
t346 = Ifges(6,4) * t182;
t85 = -Ifges(6,2) * t181 + Ifges(6,6) * t196 + t346;
t391 = (-t237 * t37 - t241 * t36) * mrSges(6,3) + t271 * t372 + t277 * t373 + t85 * t365 + t82 * t364 + t398 + t394 * t370 + t397 * t196 / 0.2e1;
t32 = qJ(6) * t196 + t37;
t390 = -t194 * mrSges(5,1) - t36 * mrSges(6,1) + t31 * mrSges(7,1) + t37 * mrSges(6,2) - t32 * mrSges(7,3) + t404 + t412;
t293 = Ifges(4,5) * qJD(3) / 0.2e1;
t29 = qJD(4) * t102 + t129 * t238 - t242 * t383;
t35 = mrSges(6,1) * t99 + mrSges(6,2) * t98;
t388 = m(6) * t29 + t35;
t195 = pkin(5) * t312 - qJ(6) * t311 - qJD(6) * t237;
t332 = t202 * t241;
t333 = t202 * t237;
t286 = -pkin(5) * t333 + qJ(6) * t332;
t387 = -t102 - t286 + t195;
t381 = Ifges(5,5) / 0.2e1;
t380 = -Ifges(5,6) / 0.2e1;
t379 = t98 / 0.2e1;
t378 = -t99 / 0.2e1;
t377 = t99 / 0.2e1;
t374 = t159 / 0.2e1;
t371 = -t182 / 0.2e1;
t369 = -t196 / 0.2e1;
t367 = t202 / 0.2e1;
t363 = -t241 / 0.2e1;
t362 = t241 / 0.2e1;
t360 = pkin(3) * t238;
t359 = pkin(3) * t242;
t358 = pkin(5) * t203;
t47 = -mrSges(7,2) * t99 + mrSges(7,3) * t159;
t50 = -mrSges(6,2) * t159 - mrSges(6,3) * t99;
t353 = t47 + t50;
t48 = mrSges(6,1) * t159 - mrSges(6,3) * t98;
t49 = -t159 * mrSges(7,1) + t98 * mrSges(7,2);
t352 = t49 - t48;
t351 = mrSges(6,3) * t181;
t350 = mrSges(6,3) * t182;
t349 = Ifges(5,1) * t203;
t348 = Ifges(4,4) * t239;
t193 = Ifges(5,4) * t202;
t329 = t235 * t240;
t198 = -t239 * t329 + t327;
t199 = t236 * t239 + t243 * t329;
t261 = t242 * t198 - t199 * t238;
t340 = t261 * t29;
t338 = t384 * t29;
t337 = t202 * mrSges(5,3);
t336 = t203 * mrSges(5,3);
t162 = pkin(4) * t203 + pkin(10) * t202;
t46 = t241 * t101 + t237 * t162;
t334 = Ifges(4,6) * qJD(3);
t331 = t214 * t243;
t328 = t235 * t244;
t326 = t237 * t242;
t325 = t241 * t242;
t106 = t175 * t242 - t165;
t316 = qJD(2) * t239;
t138 = pkin(3) * t316 + t162;
t44 = t241 * t106 + t237 * t138;
t323 = -mrSges(5,1) * t309 + mrSges(6,1) * t181 + mrSges(6,2) * t182 + t336;
t123 = -mrSges(7,2) * t181 + mrSges(7,3) * t196;
t124 = -mrSges(6,2) * t196 - t351;
t322 = t123 + t124;
t125 = mrSges(6,1) * t196 - t350;
t126 = -mrSges(7,1) * t196 + mrSges(7,2) * t182;
t321 = -t125 + t126;
t315 = qJD(2) * t240;
t314 = qJD(2) * t243;
t310 = qJD(2) * qJD(3);
t307 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t306 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t305 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t303 = t237 * t328;
t114 = mrSges(7,1) * t181 - mrSges(7,3) * t182;
t300 = t114 + t323;
t296 = t235 * t315;
t295 = t244 * t317;
t292 = -t334 / 0.2e1;
t105 = t175 * t238 + t166;
t285 = mrSges(6,1) * t241 - mrSges(6,2) * t237;
t283 = mrSges(7,1) * t241 + mrSges(7,3) * t237;
t276 = Ifges(6,2) * t241 + t345;
t270 = -Ifges(7,3) * t241 + t343;
t269 = pkin(5) * t241 + qJ(6) * t237;
t268 = pkin(5) * t237 - qJ(6) * t241;
t45 = -t101 * t237 + t162 * t241;
t43 = -t106 * t237 + t138 * t241;
t147 = -t214 * t313 + t320;
t148 = -qJD(3) * t331 + (-qJD(3) * t236 - t295) * t318;
t263 = t147 * t243 - t148 * t239;
t103 = t163 * t241 - t184 * t237;
t144 = t198 * t238 + t199 * t242;
t216 = -pkin(4) - t269;
t119 = t144 * t237 + t241 * t328;
t186 = t294 + t331;
t215 = -qJD(2) * pkin(2) - t297;
t257 = t186 * mrSges(4,3) + t334 / 0.2e1 + (t243 * Ifges(4,2) + t348) * qJD(2) / 0.2e1 - t215 * mrSges(4,1);
t185 = -t214 * t239 + t226;
t233 = Ifges(4,4) * t314;
t256 = t215 * mrSges(4,2) + Ifges(4,1) * t316 / 0.2e1 + t233 / 0.2e1 + t293 - t185 * mrSges(4,3);
t249 = t7 * mrSges(6,1) - t4 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3);
t247 = t353 * t241 + t352 * t237 + (-t237 * t322 + t241 * t321) * qJD(5) + m(7) * (-t312 * t32 + t392) + m(6) * (-t311 * t36 - t312 * t37 + t393);
t140 = -t193 + t349 + t401;
t19 = t98 * Ifges(7,5) + t159 * Ifges(7,6) + t99 * Ifges(7,3);
t20 = t98 * Ifges(6,4) - t99 * Ifges(6,2) + t159 * Ifges(6,6);
t83 = t182 * Ifges(6,5) - t181 * Ifges(6,6) + t196 * Ifges(6,3);
t84 = t182 * Ifges(7,4) + t196 * Ifges(7,2) + t181 * Ifges(7,6);
t9 = pkin(5) * t99 - qJ(6) * t98 - qJD(6) * t182 + t29;
t246 = -t28 * mrSges(5,2) - Ifges(5,5) * t158 - Ifges(5,6) * t159 - t101 * t337 + t19 * t363 + t20 * t362 + t270 * t377 + t276 * t378 - t9 * t283 + t395 * t379 + t396 * t374 + (-t193 + t140) * t367 + t402 * t364 + (-t85 / 0.2e1 + t82 / 0.2e1) * t333 + (-mrSges(5,1) - t285) * t29 + (-t277 * t372 - t271 * t373 + t401 / 0.2e1 + t399 - t394 * t371 - t397 * t369 + t398) * t202 + (Ifges(6,6) * t372 + Ifges(7,6) * t373 + t412 - Ifges(5,2) * t367 + t403 * t371 + (Ifges(6,3) + Ifges(7,2)) * t369 + t390) * t203 + (-t332 * t36 - t333 * t37 + t393) * mrSges(6,3) - (-Ifges(5,1) * t202 - t347 + t83 + t84) * t203 / 0.2e1 + (t31 * t332 + (-t312 - t333) * t32 + t392) * mrSges(7,2) + t391 * qJD(5) + (t332 / 0.2e1 + t311 / 0.2e1) * t389;
t245 = qJD(2) ^ 2;
t221 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t314;
t220 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t316;
t206 = t216 - t359;
t205 = (mrSges(4,1) * t239 + mrSges(4,2) * t243) * t310;
t192 = t203 * qJ(6);
t189 = qJD(4) * t360 + t195;
t187 = -mrSges(5,2) * t309 - t337;
t173 = -qJD(3) * t199 - t239 * t295;
t172 = qJD(3) * t198 + t243 * t295;
t156 = Ifges(7,2) * t159;
t154 = Ifges(6,3) * t159;
t120 = t144 * t241 - t303;
t113 = pkin(5) * t182 + qJ(6) * t181;
t108 = t208 * t268 - t384;
t96 = Ifges(7,4) * t98;
t95 = Ifges(6,5) * t98;
t94 = Ifges(6,6) * t99;
t93 = Ifges(7,6) * t99;
t81 = mrSges(5,1) * t159 - mrSges(5,2) * t158;
t74 = -pkin(5) * t207 - t103;
t69 = qJ(6) * t207 + t385;
t54 = t105 + t286;
t53 = qJD(4) * t144 + t172 * t238 - t242 * t173;
t52 = qJD(4) * t261 + t172 * t242 + t173 * t238;
t42 = -t45 - t358;
t41 = t192 + t46;
t39 = -t43 - t358;
t38 = t192 + t44;
t34 = mrSges(7,1) * t99 - mrSges(7,3) * t98;
t25 = -qJD(5) * t303 + t144 * t311 + t237 * t52 - t241 * t296;
t24 = -qJD(5) * t119 + t237 * t296 + t241 * t52;
t18 = -t268 * t170 + (qJD(5) * t269 - qJD(6) * t241) * t208 + t110;
t1 = [-t144 * t159 * mrSges(5,3) + t172 * t221 + t173 * t220 + t52 * t187 + t321 * t25 + t322 * t24 + t353 * t120 + t352 * t119 + (-t198 * t243 - t199 * t239) * mrSges(4,3) * t310 + t300 * t53 - (-mrSges(5,3) * t158 + t34 + t35) * t261 + ((-mrSges(3,2) * t245 - t205 - t81) * t244 + (-mrSges(3,1) * t245 + (t161 + qJD(2) * (-mrSges(4,1) * t243 + mrSges(4,2) * t239)) * qJD(2)) * t240) * t235 + m(5) * (-t101 * t53 + t102 * t52 - t340 + t144 * t28 + (t194 * t315 - t197 * t244) * t235) + m(4) * (t147 * t199 + t148 * t198 + t172 * t186 + t173 * t185 + (t215 - t297) * t296) + m(7) * (t119 * t4 + t120 * t2 + t24 * t32 + t25 * t31 - t261 * t9 + t40 * t53) + m(6) * (-t119 * t7 + t120 * t6 + t24 * t37 - t25 * t36 + t53 * t89 - t340); t385 * t50 + t108 * t34 + t18 * t114 + t103 * t48 + t69 * t47 + t74 * t49 + (-t101 * t407 + t102 * t405 + t184 * t28 + t194 * t308 + t197 * t232 - t338) * m(5) + t263 * mrSges(4,3) + (-t170 * t381 + t171 * t380 + (-pkin(8) * t220 + t256 + t293) * t243 + (-pkin(8) * t221 + pkin(3) * t161 + t292 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t314 - t257) * t239) * qJD(3) - t384 * t35 + (t158 * t384 - t159 * t184) * mrSges(5,3) + (t19 * t364 + t20 * t365 + t197 * mrSges(5,2) + t9 * t282 + t277 * t378 + t271 * t377 - Ifges(5,1) * t158 - Ifges(5,4) * t159 + (mrSges(5,3) + t284) * t29 + (-t237 * t6 - t241 * t7) * mrSges(6,3) + (-t2 * t237 + t241 * t4) * mrSges(7,2) + (t85 * t363 + t270 * t373 + t276 * t372 + t89 * t285 + t40 * t283 + (t237 * t36 - t241 * t37) * mrSges(6,3) + (-t237 * t31 - t241 * t32) * mrSges(7,2) + t395 * t371 + t396 * t369 + t389 * t365) * qJD(5) + t394 * t379 + t397 * t374 + (qJD(5) * t82 + t402) * t362) * t208 + (-t28 * mrSges(5,3) + t95 / 0.2e1 - t94 / 0.2e1 + t154 / 0.2e1 + t96 / 0.2e1 + t156 / 0.2e1 + t93 / 0.2e1 + Ifges(5,4) * t158 + t197 * mrSges(5,1) + t305 * t99 + t307 * t98 + (Ifges(5,2) + t306) * t159 + t249) * t207 + (t83 / 0.2e1 + t84 / 0.2e1 + qJD(4) * t380 - t102 * mrSges(5,3) + t306 * t196 + t307 * t182 + t305 * t181 - t390 - t404) * t171 + m(4) * ((-t185 * t243 - t186 * t239) * qJD(3) + t263) * pkin(8) + t408 * t124 + t409 * t123 + (t103 * t7 + t36 * t410 + t37 * t408 + t385 * t6 + t407 * t89 - t338) * m(6) + t410 * t125 + t411 * t126 + (t108 * t9 + t2 * t69 + t4 * t74 + (-t179 + t18) * t40 + t409 * t32 + t411 * t31) * m(7) - pkin(2) * t205 - (t140 / 0.2e1 + t391 - t101 * mrSges(5,3) + t399 + t349 / 0.2e1 - t193 / 0.2e1 + (-t237 * t32 + t241 * t31) * mrSges(7,2) + t389 * t362 + qJD(4) * t381) * t170 + t232 * t81 + t405 * t187 + ((t220 * t239 - t221 * t243) * t244 - (pkin(2) * t315 + (-t185 * t239 + t186 * t243) * t244) * m(4) + (-t215 * m(4) + t406) * t240) * t319 + (-0.3e1 / 0.2e1 * t239 ^ 2 + 0.3e1 / 0.2e1 * t243 ^ 2) * Ifges(4,4) * t310 - t300 * t179 + t323 * t110; -t147 * mrSges(4,2) + t148 * mrSges(4,1) - t38 * t123 - t44 * t124 - t43 * t125 - t39 * t126 + t246 - m(6) * (t105 * t89 + t36 * t43 + t37 * t44) - m(7) * (t31 * t39 + t32 * t38 + t40 * t54) - t323 * t105 + m(7) * (t189 * t40 + t206 * t9) + t102 * t336 - t106 * t187 + t206 * t34 + t186 * t220 - t185 * t221 - m(5) * (-t101 * t105 + t102 * t106) + (m(5) * (t238 * t28 - t242 * t29) + (t158 * t242 - t159 * t238) * mrSges(5,3) + (t323 * t238 + (t237 * t321 + t241 * t322 + t187) * t242 + m(6) * (t238 * t89 + t325 * t37 - t326 * t36) + m(7) * (t31 * t326 + t32 * t325) + m(5) * (-t101 * t238 + t102 * t242)) * qJD(4)) * pkin(3) + t247 * (pkin(10) + t360) + ((-t233 / 0.2e1 + t293 - t256) * t243 + (t292 + (t348 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t243) * qJD(2) + t406 * pkin(3) + t257) * t239) * qJD(2) + (-t54 + t189) * t114 + t388 * (-pkin(4) - t359); t387 * t114 - m(6) * (t102 * t89 + t36 * t45 + t37 * t46) - t41 * t123 - t46 * t124 - t45 * t125 - t42 * t126 + t246 + (-t323 + t336) * t102 - t101 * t187 + t216 * t34 + t247 * pkin(10) - t388 * pkin(4) + (t216 * t9 - t31 * t42 - t32 * t41 + t387 * t40) * m(7); t156 + t154 + qJD(6) * t123 - t113 * t114 + qJ(6) * t47 - pkin(5) * t49 + t249 + t96 + t95 - t94 + (t181 * t31 + t182 * t32) * mrSges(7,2) + t93 + (-t322 - t351) * t36 + (-t321 + t350) * t37 - t40 * (mrSges(7,1) * t182 + mrSges(7,3) * t181) - t89 * (mrSges(6,1) * t182 - mrSges(6,2) * t181) + t85 * t370 + (Ifges(7,3) * t182 - t339) * t373 + (-t403 * t181 - t413 * t182) * t369 + (-pkin(5) * t4 + qJ(6) * t2 - t113 * t40 - t31 * t37 + t32 * t386) * m(7) + (-Ifges(6,2) * t182 - t178 + t389) * t372 + (-t414 * t181 + t177 - t346 + t82) * t371; t182 * t114 - t196 * t123 + 0.2e1 * (t4 / 0.2e1 + t40 * t370 + t32 * t369) * m(7) + t49;];
tauc  = t1(:);
