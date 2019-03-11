% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRRRP3
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRRP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:06:23
% EndTime: 2019-03-09 00:06:57
% DurationCPUTime: 17.39s
% Computational Cost: add. (8027->608), mult. (20249->828), div. (0->0), fcn. (14444->10), ass. (0->276)
t236 = sin(qJ(2));
t231 = sin(pkin(6));
t298 = qJD(1) * t231;
t279 = t236 * t298;
t206 = qJD(2) * pkin(8) + t279;
t235 = sin(qJ(3));
t239 = cos(qJ(3));
t232 = cos(pkin(6));
t297 = qJD(1) * t232;
t162 = -t206 * t235 + t239 * t297;
t148 = -qJD(3) * pkin(3) - t162;
t234 = sin(qJ(4));
t238 = cos(qJ(4));
t292 = qJD(3) * t238;
t296 = qJD(2) * t235;
t193 = -t234 * t296 + t292;
t112 = -pkin(4) * t193 + t148;
t237 = cos(qJ(5));
t294 = qJD(2) * t239;
t222 = qJD(4) - t294;
t194 = qJD(3) * t234 + t238 * t296;
t163 = t206 * t239 + t235 * t297;
t149 = qJD(3) * pkin(9) + t163;
t208 = -pkin(3) * t239 - pkin(9) * t235 - pkin(2);
t240 = cos(qJ(2));
t278 = t240 * t298;
t165 = qJD(2) * t208 - t278;
t94 = -t149 * t234 + t165 * t238;
t73 = -pkin(10) * t194 + t94;
t64 = pkin(4) * t222 + t73;
t233 = sin(qJ(5));
t95 = t149 * t238 + t165 * t234;
t74 = pkin(10) * t193 + t95;
t68 = t233 * t74;
t27 = t237 * t64 - t68;
t131 = t193 * t233 + t194 * t237;
t384 = qJ(6) * t131;
t18 = t27 - t384;
t213 = qJD(5) + t222;
t13 = pkin(5) * t213 + t18;
t70 = t237 * t74;
t28 = t233 * t64 + t70;
t265 = t193 * t237 - t194 * t233;
t368 = qJ(6) * t265;
t19 = t28 + t368;
t286 = qJD(2) * qJD(3);
t267 = t235 * t286;
t219 = Ifges(7,3) * t267;
t220 = Ifges(6,3) * t267;
t341 = -t213 / 0.2e1;
t351 = t131 / 0.2e1;
t352 = -t131 / 0.2e1;
t285 = qJD(3) * qJD(4);
t290 = qJD(4) * t234;
t291 = qJD(3) * t239;
t155 = t238 * t285 + (-t235 * t290 + t238 * t291) * qJD(2);
t289 = qJD(4) * t238;
t383 = t234 * t291 + t235 * t289;
t156 = -qJD(2) * t383 - t234 * t285;
t51 = qJD(5) * t265 + t155 * t237 + t156 * t233;
t308 = t231 * t240;
t276 = qJD(2) * t308;
t246 = qJD(1) * (qJD(3) * t232 + t276);
t293 = qJD(3) * t235;
t119 = -t206 * t293 + t239 * t246;
t263 = pkin(3) * t235 - pkin(9) * t239;
t204 = t263 * qJD(3);
t159 = (t204 + t279) * qJD(2);
t37 = -qJD(4) * t95 - t119 * t234 + t159 * t238;
t26 = pkin(4) * t267 - pkin(10) * t155 + t37;
t36 = t119 * t238 - t149 * t290 + t159 * t234 + t165 * t289;
t31 = pkin(10) * t156 + t36;
t6 = -qJD(5) * t28 - t233 * t31 + t237 * t26;
t2 = pkin(5) * t267 - qJ(6) * t51 - qJD(6) * t131 + t6;
t287 = qJD(5) * t237;
t288 = qJD(5) * t233;
t5 = t233 * t26 + t237 * t31 + t287 * t64 - t288 * t74;
t52 = -qJD(5) * t131 - t155 * t233 + t156 * t237;
t3 = qJ(6) * t52 + qJD(6) * t265 + t5;
t365 = mrSges(6,1) * t6 + mrSges(7,1) * t2 - mrSges(6,2) * t5 - mrSges(7,2) * t3;
t378 = -t265 / 0.2e1;
t395 = Ifges(6,5) + Ifges(7,5);
t397 = Ifges(6,1) + Ifges(7,1);
t396 = Ifges(6,4) + Ifges(7,4);
t409 = t396 * t265;
t386 = t131 * t397 + t213 * t395 + t409;
t393 = Ifges(6,6) + Ifges(7,6);
t394 = Ifges(6,2) + Ifges(7,2);
t406 = t131 * t396;
t387 = t213 * t393 + t265 * t394 + t406;
t47 = Ifges(7,6) * t52;
t48 = Ifges(6,6) * t52;
t49 = Ifges(7,5) * t51;
t50 = Ifges(6,5) * t51;
t65 = -pkin(5) * t265 + qJD(6) + t112;
t413 = t219 + t220 + t47 + t48 + t49 + t50 + t365 + (-t131 * t393 + t265 * t395) * t341 + (t13 * t265 + t131 * t19) * mrSges(7,3) + (t131 * t28 + t265 * t27) * mrSges(6,3) - t112 * (mrSges(6,1) * t131 + mrSges(6,2) * t265) - t65 * (mrSges(7,1) * t131 + mrSges(7,2) * t265) + t387 * t351 + (-t131 * t394 + t386 + t409) * t378 + (t265 * t397 - t406) * t352;
t201 = t263 * qJD(2);
t109 = -t162 * t234 + t201 * t238;
t306 = t238 * t239;
t252 = pkin(4) * t235 - pkin(10) * t306;
t357 = -pkin(10) - pkin(9);
t280 = qJD(4) * t357;
t412 = -qJD(2) * t252 + t238 * t280 - t109;
t110 = t162 * t238 + t201 * t234;
t275 = t234 * t294;
t411 = -pkin(10) * t275 - t234 * t280 + t110;
t410 = -t394 * t52 / 0.2e1 - t396 * t51 / 0.2e1 - t393 * t267 / 0.2e1;
t253 = t233 * t234 - t237 * t238;
t366 = qJD(4) + qJD(5);
t136 = t366 * t253;
t248 = t253 * t239;
t169 = qJD(2) * t248;
t408 = t136 - t169;
t196 = t233 * t238 + t234 * t237;
t137 = t366 * t196;
t249 = t196 * t239;
t168 = qJD(2) * t249;
t407 = t137 - t168;
t405 = -t294 / 0.2e1;
t192 = t238 * t208;
t332 = pkin(10) * t235;
t334 = pkin(8) * t234;
t132 = -t238 * t332 + t192 + (-pkin(4) - t334) * t239;
t224 = pkin(8) * t306;
t167 = t208 * t234 + t224;
t307 = t234 * t235;
t143 = -pkin(10) * t307 + t167;
t299 = t204 * t238 + t293 * t334;
t77 = t252 * qJD(3) + (-t224 + (-t208 + t332) * t234) * qJD(4) + t299;
t103 = t234 * t204 + t208 * t289 + (-t235 * t292 - t239 * t290) * pkin(8);
t83 = -pkin(10) * t383 + t103;
t11 = t132 * t287 - t143 * t288 + t233 * t77 + t237 * t83;
t305 = t239 * t240;
t157 = (-t234 * t305 + t236 * t238) * t298;
t158 = (t234 * t236 + t238 * t305) * t298;
t97 = t157 * t233 + t158 * t237;
t404 = t11 - t97;
t79 = t132 * t233 + t143 * t237;
t12 = -qJD(5) * t79 - t233 * t83 + t237 * t77;
t96 = t157 * t237 - t158 * t233;
t403 = t12 - t96;
t401 = t267 * t395 + t396 * t52 + t397 * t51;
t211 = t357 * t234;
t212 = t357 * t238;
t154 = t211 * t233 - t212 * t237;
t391 = -qJD(5) * t154 + t233 * t411 + t237 * t412;
t390 = t211 * t287 + t212 * t288 + t233 * t412 - t237 * t411;
t399 = -Ifges(4,1) / 0.2e1;
t339 = -t222 / 0.2e1;
t398 = Ifges(4,4) * t405;
t392 = Ifges(7,3) + Ifges(6,3);
t389 = -qJ(6) * t407 - qJD(6) * t253 + t390;
t388 = -pkin(5) * t296 + qJ(6) * t408 - qJD(6) * t196 + t391;
t133 = pkin(4) * t275 + t163;
t385 = pkin(4) * t290 + pkin(5) * t407 - t133;
t207 = -qJD(2) * pkin(2) - t278;
t325 = Ifges(5,4) * t194;
t116 = Ifges(5,2) * t193 + Ifges(5,6) * t222 + t325;
t189 = Ifges(5,4) * t193;
t117 = Ifges(5,1) * t194 + Ifges(5,5) * t222 + t189;
t254 = t234 * t95 + t238 * t94;
t256 = Ifges(5,5) * t238 - Ifges(5,6) * t234;
t323 = Ifges(5,4) * t238;
t258 = -Ifges(5,2) * t234 + t323;
t324 = Ifges(5,4) * t234;
t260 = Ifges(5,1) * t238 - t324;
t261 = mrSges(5,1) * t234 + mrSges(5,2) * t238;
t336 = t238 / 0.2e1;
t337 = -t234 / 0.2e1;
t338 = t222 / 0.2e1;
t344 = t194 / 0.2e1;
t242 = -t254 * mrSges(5,3) + t116 * t337 + t117 * t336 + t148 * t261 + t193 * t258 / 0.2e1 + t260 * t344 + t256 * t338;
t313 = Ifges(4,5) * qJD(3);
t382 = t162 * mrSges(4,3) + t296 * t399 + t398 - t313 / 0.2e1 - t207 * mrSges(4,2) - t242;
t377 = Ifges(5,3) * t339;
t329 = qJD(2) / 0.2e1;
t76 = -mrSges(6,1) * t265 + mrSges(6,2) * t131;
t268 = m(6) * t112 + t76;
t172 = t253 * t235;
t367 = -t234 * t37 + t238 * t36;
t363 = -mrSges(5,1) * t37 + mrSges(5,2) * t36 - Ifges(5,5) * t155 - Ifges(5,6) * t156;
t283 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t284 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t312 = Ifges(4,6) * qJD(3);
t326 = Ifges(4,4) * t235;
t362 = -t284 * t131 - (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t213 - t283 * t265 + t163 * mrSges(4,3) + t19 * mrSges(7,2) + t28 * mrSges(6,2) + t95 * mrSges(5,2) + t377 - t194 * Ifges(5,5) - t193 * Ifges(5,6) + t312 / 0.2e1 + (t239 * Ifges(4,2) + t326) * t329 - t13 * mrSges(7,1) - t207 * mrSges(4,1) - t27 * mrSges(6,1) - t94 * mrSges(5,1) + t393 * t378 + t395 * t352 + t392 * t341;
t361 = t51 / 0.2e1;
t360 = t52 / 0.2e1;
t354 = t265 / 0.2e1;
t350 = t155 / 0.2e1;
t349 = t156 / 0.2e1;
t346 = -t193 / 0.2e1;
t345 = -t194 / 0.2e1;
t340 = t213 / 0.2e1;
t333 = pkin(8) * t239;
t44 = -mrSges(7,2) * t267 + mrSges(7,3) * t52;
t45 = -mrSges(6,2) * t267 + mrSges(6,3) * t52;
t327 = t44 + t45;
t33 = t237 * t73 - t68;
t120 = t206 * t291 + t235 * t246;
t309 = t231 * t236;
t174 = -t232 * t239 + t235 * t309;
t311 = t120 * t174;
t105 = -mrSges(7,2) * t213 + mrSges(7,3) * t265;
t106 = -mrSges(6,2) * t213 + mrSges(6,3) * t265;
t304 = t105 + t106;
t107 = mrSges(7,1) * t213 - mrSges(7,3) * t131;
t108 = mrSges(6,1) * t213 - mrSges(6,3) * t131;
t303 = t107 + t108;
t300 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t193 + mrSges(5,2) * t194 + mrSges(4,3) * t296;
t205 = pkin(4) * t307 + pkin(8) * t235;
t295 = qJD(2) * t236;
t281 = mrSges(4,3) * t294;
t164 = pkin(4) * t383 + pkin(8) * t291;
t227 = -pkin(4) * t238 - pkin(3);
t277 = t231 * t295;
t272 = t313 / 0.2e1;
t271 = -t312 / 0.2e1;
t20 = -mrSges(7,1) * t52 + mrSges(7,2) * t51;
t266 = 0.3e1 / 0.2e1 * t294;
t32 = -t233 * t73 - t70;
t78 = t132 * t237 - t143 * t233;
t153 = t211 * t237 + t212 * t233;
t262 = mrSges(5,1) * t238 - mrSges(5,2) * t234;
t259 = Ifges(5,1) * t234 + t323;
t257 = Ifges(5,2) * t238 + t324;
t255 = Ifges(5,5) * t234 + Ifges(5,6) * t238;
t175 = t232 * t235 + t239 * t309;
t141 = -t175 * t234 - t238 * t308;
t251 = -t175 * t238 + t234 * t308;
t81 = t141 * t237 + t233 * t251;
t82 = t141 * t233 - t237 * t251;
t250 = -m(4) * t162 + m(5) * t148 + t300;
t84 = -pkin(4) * t156 + t120;
t241 = qJD(2) ^ 2;
t226 = pkin(4) * t237 + pkin(5);
t221 = Ifges(5,3) * t267;
t210 = -qJD(3) * mrSges(4,2) + t281;
t200 = (-mrSges(4,1) * t239 + mrSges(4,2) * t235) * qJD(2);
t185 = (mrSges(4,1) * t235 + mrSges(4,2) * t239) * t286;
t171 = t196 * t235;
t170 = pkin(5) * t253 + t227;
t166 = -t234 * t333 + t192;
t161 = mrSges(5,1) * t222 - mrSges(5,3) * t194;
t160 = -mrSges(5,2) * t222 + mrSges(5,3) * t193;
t140 = -qJD(3) * t174 + t239 * t276;
t139 = qJD(3) * t175 + t235 * t276;
t134 = pkin(5) * t171 + t205;
t127 = -mrSges(5,2) * t267 + mrSges(5,3) * t156;
t126 = mrSges(5,1) * t267 - mrSges(5,3) * t155;
t114 = -qJ(6) * t253 + t154;
t113 = -qJ(6) * t196 + t153;
t104 = -qJD(4) * t167 + t299;
t102 = pkin(4) * t194 + pkin(5) * t131;
t100 = -mrSges(5,1) * t156 + mrSges(5,2) * t155;
t88 = Ifges(5,1) * t155 + Ifges(5,4) * t156 + Ifges(5,5) * t267;
t87 = Ifges(5,4) * t155 + Ifges(5,2) * t156 + Ifges(5,6) * t267;
t86 = -qJD(3) * t249 + t172 * t366;
t85 = -qJD(3) * t248 - t137 * t235;
t75 = -mrSges(7,1) * t265 + mrSges(7,2) * t131;
t67 = qJD(4) * t141 + t140 * t238 + t234 * t277;
t66 = qJD(4) * t251 - t140 * t234 + t238 * t277;
t55 = -pkin(5) * t86 + t164;
t54 = -qJ(6) * t171 + t79;
t53 = -pkin(5) * t239 + qJ(6) * t172 + t78;
t43 = mrSges(6,1) * t267 - mrSges(6,3) * t51;
t42 = mrSges(7,1) * t267 - mrSges(7,3) * t51;
t29 = -pkin(5) * t52 + t84;
t23 = t33 - t384;
t22 = t32 - t368;
t21 = -mrSges(6,1) * t52 + mrSges(6,2) * t51;
t10 = -qJD(5) * t82 - t233 * t67 + t237 * t66;
t9 = qJD(5) * t81 + t233 * t66 + t237 * t67;
t8 = qJ(6) * t86 - qJD(6) * t171 + t11;
t7 = pkin(5) * t293 - qJ(6) * t85 + qJD(6) * t172 + t12;
t1 = [-t175 * mrSges(4,3) * t267 + t141 * t126 - t251 * t127 + t140 * t210 + t67 * t160 + t66 * t161 + t304 * t9 + t327 * t82 + (t43 + t42) * t81 + t303 * t10 + ((-mrSges(3,2) * t241 - t185) * t240 + (-mrSges(3,1) * t241 + qJD(2) * t200) * t236) * t231 + (qJD(3) * t281 + t100 + t20 + t21) * t174 + (t75 + t76 + t300) * t139 + m(5) * (t139 * t148 + t141 * t37 - t251 * t36 + t66 * t94 + t67 * t95 + t311) + m(6) * (t10 * t27 + t112 * t139 + t174 * t84 + t28 * t9 + t5 * t82 + t6 * t81) + m(7) * (t10 * t13 + t139 * t65 + t174 * t29 + t19 * t9 + t2 * t81 + t3 * t82) + m(4) * (t119 * t175 + t311 - t139 * t162 + t140 * t163 + (t207 - t278) * t277); (-t171 * t5 + t172 * t6 - t27 * t85 + t28 * t86) * mrSges(6,3) + t29 * (mrSges(7,1) * t171 - mrSges(7,2) * t172) + t84 * (mrSges(6,1) * t171 - mrSges(6,2) * t172) + (-t13 * t85 - t171 * t3 + t172 * t2 + t19 * t86) * mrSges(7,3) + m(5) * (t103 * t95 + t104 * t94 + t166 * t37 + t167 * t36) - m(5) * (t157 * t94 + t158 * t95) - m(7) * (t13 * t96 + t19 * t97) + (-t171 * t396 - t172 * t397) * t361 + m(7) * (t13 * t7 + t134 * t29 + t19 * t8 + t2 * t53 + t3 * t54 + t55 * t65) + (Ifges(4,4) * t266 + t250 * pkin(8) + t272 - t382) * t291 + (t119 * mrSges(4,3) - t49 / 0.2e1 - t47 / 0.2e1 - t219 / 0.2e1 - t50 / 0.2e1 - t48 / 0.2e1 - t220 / 0.2e1 - t221 / 0.2e1 - t283 * t52 - t284 * t51 + t363 - t365) * t239 + (t260 * t350 + t258 * t349 + t87 * t337 + t88 * t336 + (-t234 * t36 - t238 * t37) * mrSges(5,3) + (mrSges(4,3) + t261) * t120 + (mrSges(4,2) * t295 + (-m(7) * t65 - t250 - t268 - t75) * t240) * t298 + (t117 * t337 - t238 * t116 / 0.2e1 + t148 * t262 + t257 * t346 + t259 * t345 + t255 * t339 + (t234 * t94 - t238 * t95) * mrSges(5,3)) * qJD(4) + (Ifges(4,1) * t266 - 0.3e1 / 0.2e1 * Ifges(4,2) * t294 + (t405 + t338) * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(4,4) * t296 + t271 - t362 + (-t171 * t393 - t172 * t395 + t235 * t256 - t239 * t392) * t329) * qJD(3) + (t100 + (m(4) + m(5)) * t120 + (-m(4) * t163 - t210) * qJD(3)) * pkin(8)) * t235 + (t104 - t157) * t161 - t401 * t172 / 0.2e1 + t403 * t108 + (t112 * t164 + t205 * t84 + t27 * t403 + t28 * t404 + t5 * t79 + t6 * t78) * m(6) + t404 * t106 + (t393 * t86 + t395 * t85) * t340 + (t394 * t86 + t396 * t85) * t354 + (-t171 * t394 - t172 * t396) * t360 + (t396 * t86 + t397 * t85) * t351 + t386 * t85 / 0.2e1 + t387 * t86 / 0.2e1 + (t103 - t158) * t160 + t171 * t410 + m(4) * (-pkin(2) * qJD(1) * t277 + t119 * t333) + (-t210 * t305 + (-mrSges(4,1) * t294 - t200) * t236) * t298 + (t7 - t96) * t107 - m(4) * (t163 * t239 * t278 + t207 * t279) + t53 * t42 + t54 * t44 + t55 * t75 + t78 * t43 + t79 * t45 + t65 * (-mrSges(7,1) * t86 + mrSges(7,2) * t85) + t112 * (-mrSges(6,1) * t86 + mrSges(6,2) * t85) + t134 * t20 + t164 * t76 + t166 * t126 + t167 * t127 - pkin(2) * t185 + t205 * t21 + (t8 - t97) * t105; -m(5) * (t109 * t94 + t110 * t95 + t148 * t163) + m(5) * (-pkin(3) * t120 + pkin(9) * t367) + t367 * mrSges(5,3) + (-t126 * t234 + t127 * t238) * pkin(9) + t84 * (mrSges(6,1) * t253 + mrSges(6,2) * t196) + t29 * (mrSges(7,1) * t253 + mrSges(7,2) * t196) + t401 * t196 / 0.2e1 + t385 * t75 + t253 * t410 - t300 * t163 + (t268 * t234 * pkin(4) + (-m(5) * t254 - t160 * t234 - t161 * t238) * pkin(9) + t242) * qJD(4) + t386 * (-t136 / 0.2e1 + t169 / 0.2e1) + t387 * (-t137 / 0.2e1 + t168 / 0.2e1) + t388 * t107 + t389 * t105 + (t113 * t2 + t114 * t3 + t13 * t388 + t170 * t29 + t19 * t389 + t385 * t65) * m(7) + t390 * t106 + t391 * t108 + (-t112 * t133 + t153 * t6 + t154 * t5 + t227 * t84 + t27 * t391 + t28 * t390) * m(6) + (-t136 * t395 - t137 * t393) * t340 + (-t168 * t393 - t169 * t395) * t341 + (-t136 * t396 - t137 * t394) * t354 + (-t168 * t394 - t169 * t396) * t378 + (t196 * t396 - t253 * t394) * t360 + (-t136 * t397 - t137 * t396) * t351 + (-t168 * t396 - t169 * t397) * t352 + (t196 * t397 - t253 * t396) * t361 + ((t272 + t398 + t382) * t239 + ((t326 / 0.2e1 + (Ifges(4,2) / 0.2e1 + t399) * t239) * qJD(2) + t377 + t271 + t362) * t235 + (t196 * t395 - t253 * t393 + t255) * t293 / 0.2e1) * qJD(2) - pkin(3) * t100 + t113 * t42 + t114 * t44 - t119 * mrSges(4,2) + (-t196 * t6 - t253 * t5 + t27 * t408 - t28 * t407) * mrSges(6,3) + (t13 * t408 - t19 * t407 - t196 * t2 - t253 * t3) * mrSges(7,3) + (mrSges(6,1) * t407 - mrSges(6,2) * t408) * t112 + (mrSges(7,1) * t407 - mrSges(7,2) * t408) * t65 + t87 * t336 + t257 * t349 + t259 * t350 - t133 * t76 + t153 * t43 + t154 * t45 - t110 * t160 - t109 * t161 + t170 * t20 - t162 * t210 + t227 * t21 + t234 * t88 / 0.2e1 + (-mrSges(4,1) - t262) * t120; (-Ifges(5,2) * t194 + t117 + t189) * t346 + (t193 * t94 + t194 * t95) * mrSges(5,3) - m(6) * (t27 * t32 + t28 * t33) - t363 + t221 + (t237 * t43 + t327 * t233 + (-t233 * t303 + t237 * t304) * qJD(5) + m(6) * (t233 * t5 + t237 * t6 - t27 * t288 + t28 * t287) - t268 * t194) * pkin(4) + (t2 * t226 - t102 * t65 - t13 * t22 - t19 * t23 + (-t13 * t288 + t19 * t287 + t233 * t3) * pkin(4)) * m(7) - t102 * t75 - t23 * t105 - t33 * t106 - t22 * t107 - t32 * t108 + (Ifges(5,5) * t193 - Ifges(5,6) * t194) * t339 + t116 * t344 + (Ifges(5,1) * t193 - t325) * t345 - t94 * t160 + t95 * t161 - t148 * (mrSges(5,1) * t194 + mrSges(5,2) * t193) + t226 * t42 + t413; (-t131 * t75 + t42) * pkin(5) + (-(-t13 + t18) * t19 + (-t131 * t65 + t2) * pkin(5)) * m(7) - t18 * t105 - t27 * t106 + t19 * t107 + t28 * t108 + t413; -t265 * t105 + t131 * t107 + 0.2e1 * (t29 / 0.2e1 + t19 * t378 + t13 * t351) * m(7) + t20;];
tauc  = t1(:);
