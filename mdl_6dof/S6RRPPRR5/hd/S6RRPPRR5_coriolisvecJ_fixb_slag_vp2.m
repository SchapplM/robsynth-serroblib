% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2018-11-23 16:52
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:51:47
% EndTime: 2018-11-23 16:51:57
% DurationCPUTime: 9.62s
% Computational Cost: add. (6953->639), mult. (17987->849), div. (0->0), fcn. (12340->8), ass. (0->272)
t222 = -pkin(2) - pkin(3);
t212 = pkin(4) - t222;
t217 = sin(qJ(5));
t220 = cos(qJ(5));
t162 = pkin(5) * t220 + pkin(10) * t217 + t212;
t215 = qJ(3) - pkin(9);
t218 = sin(qJ(2));
t213 = sin(pkin(6));
t282 = qJD(1) * t213;
t276 = t218 * t282;
t279 = qJD(5) * t217;
t221 = cos(qJ(2));
t214 = cos(pkin(6));
t325 = pkin(1) * t214;
t204 = t221 * t325;
t191 = qJD(1) * t204;
t153 = -pkin(8) * t276 + t191;
t177 = qJ(4) * t276;
t123 = t177 + t153;
t275 = t221 * t282;
t180 = qJ(3) * t275;
t236 = t213 * (-pkin(9) * t221 - t212 * t218);
t82 = qJD(1) * t236 + t180;
t51 = t220 * t123 + t217 * t82;
t375 = -pkin(10) * t276 - qJD(3) * t220 - qJD(6) * t162 + t215 * t279 + t51;
t179 = qJ(4) * t275;
t203 = t218 * t325;
t262 = -pkin(5) * t217 + pkin(10) * t220;
t293 = t215 * t220;
t294 = t213 * t221;
t374 = t262 * qJD(5) - qJD(6) * t293 - t179 - (-t203 + (-pkin(8) - t262) * t294) * qJD(1);
t373 = Ifges(5,4) + Ifges(4,5);
t216 = sin(qJ(6));
t219 = cos(qJ(6));
t172 = qJD(5) + t275;
t129 = -pkin(1) * t282 - pkin(2) * t275 - qJ(3) * t276;
t103 = pkin(3) * t275 + qJD(4) - t129;
t239 = (pkin(4) * t221 - pkin(9) * t218) * t213;
t70 = qJD(1) * t239 + t103;
t198 = qJD(1) * t214 + qJD(2);
t283 = pkin(8) * t294 + t203;
t154 = t283 * qJD(1);
t125 = -t179 + t154;
t176 = t198 * qJ(3);
t97 = t176 + t125;
t81 = -pkin(9) * t198 + t97;
t34 = t217 * t70 + t220 * t81;
t29 = pkin(10) * t172 + t34;
t136 = -t198 * t217 + t220 * t276;
t138 = t198 * t220 + t217 * t276;
t364 = qJD(3) - t153;
t113 = -t198 * pkin(2) + t364;
t79 = -t198 * pkin(3) + t113 - t177;
t67 = pkin(4) * t198 - t79;
t35 = pkin(5) * t138 - pkin(10) * t136 + t67;
t10 = -t216 * t29 + t219 * t35;
t11 = t216 * t35 + t219 * t29;
t250 = t10 * t219 + t11 * t216;
t257 = mrSges(7,1) * t216 + mrSges(7,2) * t219;
t33 = -t217 * t81 + t220 * t70;
t28 = -pkin(5) * t172 - t33;
t133 = qJD(6) + t138;
t90 = t136 * t219 + t172 * t216;
t326 = Ifges(7,4) * t90;
t89 = -t136 * t216 + t172 * t219;
t31 = Ifges(7,2) * t89 + Ifges(7,6) * t133 + t326;
t85 = Ifges(7,4) * t89;
t32 = Ifges(7,1) * t90 + Ifges(7,5) * t133 + t85;
t327 = t219 / 0.2e1;
t328 = t216 / 0.2e1;
t372 = t250 * mrSges(7,3) - t28 * t257 + t31 * t328 - t32 * t327;
t370 = -mrSges(4,1) - mrSges(3,1);
t369 = mrSges(4,2) + mrSges(3,3);
t368 = Ifges(5,6) + Ifges(3,6);
t367 = t216 * t375 + t374 * t219;
t366 = t374 * t216 - t219 * t375;
t352 = t283 * qJD(2);
t140 = qJD(1) * t352;
t295 = t213 * t218;
t193 = qJD(4) * t295;
t271 = qJD(2) * t282;
t268 = t221 * t271;
t87 = -qJ(4) * t268 - qJD(1) * t193 + t140;
t252 = Ifges(7,5) * t219 - Ifges(7,6) * t216;
t310 = Ifges(7,4) * t219;
t254 = -Ifges(7,2) * t216 + t310;
t311 = Ifges(7,4) * t216;
t256 = Ifges(7,1) * t219 - t311;
t330 = t133 / 0.2e1;
t338 = t90 / 0.2e1;
t340 = t89 / 0.2e1;
t365 = t252 * t330 + t254 * t340 + t256 * t338 - t372;
t280 = qJD(2) * t221;
t274 = t213 * t280;
t102 = -qJ(4) * t274 - t193 + t352;
t269 = t218 * t271;
t273 = t220 * t280;
t278 = qJD(5) * t220;
t95 = -t198 * t278 + (-t218 * t279 + t273) * t282;
t39 = qJD(6) * t89 - t216 * t269 + t219 * t95;
t234 = (t217 * t280 + t218 * t278) * t213;
t96 = qJD(1) * t234 - t198 * t279;
t25 = mrSges(7,1) * t96 - mrSges(7,3) * t39;
t40 = -qJD(6) * t90 - t216 * t95 - t219 * t269;
t26 = -mrSges(7,2) * t96 + mrSges(7,3) * t40;
t56 = -mrSges(7,2) * t133 + mrSges(7,3) * t89;
t57 = mrSges(7,1) * t133 - mrSges(7,3) * t90;
t363 = -t216 * t25 + t219 * t26 + (-m(7) * t250 - t216 * t56 - t219 * t57) * qJD(6);
t362 = Ifges(6,6) / 0.2e1;
t336 = t96 / 0.2e1;
t361 = -t198 / 0.2e1;
t360 = t198 / 0.2e1;
t357 = -t95 * Ifges(6,4) / 0.2e1;
t141 = t214 * qJ(3) + t283;
t114 = -qJ(4) * t294 + t141;
t104 = -pkin(9) * t214 + t114;
t142 = -t213 * pkin(1) - pkin(2) * t294 - qJ(3) * t295;
t115 = pkin(3) * t294 - t142;
t84 = t239 + t115;
t355 = t220 * t104 + t217 * t84;
t132 = Ifges(6,4) * t138;
t300 = t172 * Ifges(6,5);
t302 = t136 * Ifges(6,1);
t64 = -t132 + t300 + t302;
t354 = t33 * mrSges(6,3) - t64 / 0.2e1 - t300 / 0.2e1 - t67 * mrSges(6,2) + t132 / 0.2e1;
t296 = qJD(5) * t34;
t230 = qJD(2) * t236;
t194 = qJD(3) * t295;
t286 = qJ(3) * t268 + qJD(1) * t194;
t59 = qJD(1) * t230 + t286;
t171 = qJD(2) * t191;
t173 = t198 * qJD(3);
t281 = qJD(2) * t218;
t231 = (-qJD(4) * t221 + (-pkin(8) + qJ(4)) * t281) * t213;
t71 = qJD(1) * t231 + t171 + t173;
t13 = -t217 * t71 + t220 * t59 - t296;
t353 = -t10 * t216 + t11 * t219;
t27 = pkin(5) * t96 - pkin(10) * t95 - t87;
t12 = t217 * t59 + t220 * t71 + t70 * t278 - t279 * t81;
t5 = -pkin(10) * t269 + t12;
t1 = qJD(6) * t10 + t216 * t27 + t219 * t5;
t2 = -qJD(6) * t11 - t216 * t5 + t219 * t27;
t261 = t1 * t219 - t2 * t216;
t351 = t13 * mrSges(6,1) - t12 * mrSges(6,2) + Ifges(6,5) * t95 - Ifges(6,6) * t96;
t350 = -t373 * t276 / 0.2e1;
t122 = t176 + t154;
t313 = Ifges(3,4) * t218;
t348 = -(Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1 + Ifges(3,6) / 0.2e1) * t198 + t129 * mrSges(4,1) + t34 * mrSges(6,2) + t97 * mrSges(5,3) + Ifges(4,6) * t360 - (Ifges(3,2) * t221 + t313) * t282 / 0.2e1 - t172 * Ifges(6,3) - t136 * Ifges(6,5) + t138 * Ifges(6,6) - t103 * mrSges(5,1) - t122 * mrSges(4,2) - t154 * mrSges(3,3) - t33 * mrSges(6,1) - t350 + t368 * t361 - (Ifges(4,3) + Ifges(5,2)) * t275 / 0.2e1;
t185 = Ifges(3,4) * t275;
t277 = Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t347 = -(Ifges(5,5) / 0.2e1 - t277) * t198 + t103 * mrSges(5,2) + t113 * mrSges(4,2) + Ifges(5,5) * t361 + Ifges(3,1) * t276 / 0.2e1 + t185 / 0.2e1 - t129 * mrSges(4,3) - t153 * mrSges(3,3) - t79 * mrSges(5,3) + (Ifges(4,4) + Ifges(3,5)) * t360 + (-t373 * t221 + (Ifges(4,1) + Ifges(5,1)) * t218) * t282 / 0.2e1;
t284 = qJ(3) * t274 + t194;
t69 = t230 + t284;
t192 = qJD(2) * t204;
t206 = t214 * qJD(3);
t88 = t192 + t206 + t231;
t18 = -qJD(5) * t355 - t217 * t88 + t220 * t69;
t346 = Ifges(6,2) / 0.2e1;
t345 = t31 / 0.2e1;
t344 = t39 / 0.2e1;
t343 = t40 / 0.2e1;
t342 = Ifges(6,2) * t336 + t269 * t362 + t357;
t341 = -t89 / 0.2e1;
t339 = -t90 / 0.2e1;
t337 = -t96 / 0.2e1;
t335 = m(7) * t28;
t334 = pkin(1) * mrSges(3,1);
t333 = pkin(1) * mrSges(3,2);
t331 = -t133 / 0.2e1;
t157 = t214 * t220 + t217 * t295;
t329 = t157 / 0.2e1;
t38 = Ifges(7,5) * t39;
t37 = Ifges(7,6) * t40;
t322 = t89 * Ifges(7,6);
t321 = t90 * Ifges(7,5);
t320 = t95 * Ifges(6,1);
t318 = t96 * Ifges(6,4);
t317 = mrSges(4,2) - mrSges(5,3);
t16 = -mrSges(7,1) * t40 + mrSges(7,2) * t39;
t73 = -mrSges(6,1) * t269 - mrSges(6,3) * t95;
t315 = t16 - t73;
t47 = -mrSges(7,1) * t89 + mrSges(7,2) * t90;
t99 = mrSges(6,1) * t172 - mrSges(6,3) * t136;
t314 = t47 - t99;
t312 = Ifges(6,4) * t136;
t305 = t133 * Ifges(7,3);
t304 = t138 * Ifges(6,2);
t299 = t172 * Ifges(6,6);
t297 = -mrSges(5,1) * t198 - mrSges(6,1) * t138 - mrSges(6,2) * t136 - mrSges(5,3) * t276;
t292 = t217 * t221;
t291 = t219 * t220;
t290 = t220 * t221;
t289 = t198 * t370 + t276 * t369;
t147 = mrSges(5,2) * t198 - mrSges(5,3) * t275;
t149 = mrSges(4,2) * t275 + mrSges(4,3) * t198;
t288 = t149 + t147;
t159 = -pkin(8) * t295 + t204;
t7 = Ifges(7,3) * t96 + t37 + t38;
t143 = -t214 * pkin(2) - t159;
t272 = t213 * t281;
t270 = t222 * t295;
t267 = 0.3e1 / 0.2e1 * Ifges(4,5) - 0.3e1 / 0.2e1 * Ifges(3,4) + 0.3e1 / 0.2e1 * Ifges(5,4);
t264 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t52 = t96 * mrSges(6,1) + t95 * mrSges(6,2);
t260 = t1 * t216 + t2 * t219;
t258 = -mrSges(7,1) * t219 + mrSges(7,2) * t216;
t255 = Ifges(7,1) * t216 + t310;
t253 = Ifges(7,2) * t219 + t311;
t251 = Ifges(7,5) * t216 + Ifges(7,6) * t219;
t42 = pkin(10) * t294 + t355;
t158 = -t214 * t217 + t220 * t295;
t105 = -t214 * pkin(3) - qJ(4) * t295 + t143;
t94 = t214 * pkin(4) - t105;
t55 = pkin(5) * t157 - pkin(10) * t158 + t94;
t20 = t216 * t55 + t219 * t42;
t19 = -t216 * t42 + t219 * t55;
t247 = t216 * t57 - t219 * t56;
t48 = -t104 * t217 + t220 * t84;
t50 = -t123 * t217 + t220 * t82;
t242 = qJD(2) * t270;
t155 = -pkin(8) * t272 + t192;
t98 = -mrSges(6,2) * t172 - mrSges(6,3) * t138;
t241 = t247 - t98;
t111 = -t158 * t216 + t219 * t294;
t112 = t158 * t219 + t216 * t294;
t17 = -t104 * t279 + t217 * t69 + t220 * t88 + t84 * t278;
t139 = -pkin(8) * t269 + t171;
t108 = t139 + t173;
t229 = -t87 * mrSges(5,1) - t139 * mrSges(3,2) + t71 * mrSges(5,2) + t108 * mrSges(4,3) + t140 * t370;
t228 = -t302 / 0.2e1 + t354;
t30 = t305 + t321 + t322;
t63 = t299 - t304 + t312;
t226 = -t10 * mrSges(7,1) + t11 * mrSges(7,2) - t67 * mrSges(6,1) + t299 / 0.2e1 - t305 / 0.2e1 - t322 / 0.2e1 - t321 / 0.2e1 - t30 / 0.2e1 + t63 / 0.2e1 + t34 * mrSges(6,3) + t312 / 0.2e1;
t225 = (-t304 / 0.2e1 + t226) * t217;
t170 = Ifges(4,4) * t268;
t169 = Ifges(3,5) * t268;
t168 = Ifges(4,6) * t269;
t167 = mrSges(5,2) * t268;
t152 = pkin(2) * t276 - t180;
t151 = (mrSges(5,1) * t221 + mrSges(5,2) * t218) * t282;
t150 = (-mrSges(4,1) * t221 - mrSges(4,3) * t218) * t282;
t148 = -mrSges(3,2) * t198 + mrSges(3,3) * t275;
t134 = t155 + t206;
t131 = -t216 * t276 + t275 * t291;
t130 = (t216 * t290 + t218 * t219) * t282;
t128 = t162 * t216 + t215 * t291;
t127 = t162 * t219 - t216 * t293;
t126 = pkin(2) * t272 - t284;
t124 = qJD(1) * t270 + t180;
t110 = -t214 * t279 + t234;
t109 = -qJD(5) * t157 + t213 * t273;
t106 = pkin(2) * t269 - t286;
t101 = t242 + t284;
t86 = qJD(1) * t242 + t286;
t77 = pkin(5) * t136 + pkin(10) * t138;
t74 = mrSges(6,2) * t269 - mrSges(6,3) * t96;
t54 = -qJD(6) * t112 - t109 * t216 - t219 * t272;
t53 = qJD(6) * t111 + t109 * t219 - t216 * t272;
t46 = -Ifges(6,5) * t269 - t318 + t320;
t43 = pkin(5) * t276 - t50;
t41 = -pkin(5) * t294 - t48;
t36 = pkin(5) * t110 - pkin(10) * t109 - t102;
t22 = t216 * t77 + t219 * t33;
t21 = -t216 * t33 + t219 * t77;
t15 = pkin(5) * t272 - t18;
t14 = -pkin(10) * t272 + t17;
t9 = t39 * Ifges(7,1) + t40 * Ifges(7,4) + t96 * Ifges(7,5);
t8 = t39 * Ifges(7,4) + t40 * Ifges(7,2) + t96 * Ifges(7,6);
t6 = pkin(5) * t269 - t13;
t4 = -qJD(6) * t20 - t14 * t216 + t219 * t36;
t3 = qJD(6) * t19 + t14 * t219 + t216 * t36;
t23 = [((t143 * mrSges(4,2) - t159 * mrSges(3,3) - t142 * mrSges(4,3) - t105 * mrSges(5,3) + (-t221 * t267 - 0.2e1 * t333) * t213 + (-Ifges(5,5) + t277) * t214) * t280 + (t142 * mrSges(4,1) - t115 * mrSges(5,1) + t114 * mrSges(5,3) - t141 * mrSges(4,2) - t283 * mrSges(3,3) - Ifges(6,5) * t158 / 0.2e1 + Ifges(6,6) * t329 + (t218 * t267 - 0.2e1 * t334) * t213 + (Ifges(4,6) / 0.2e1 - t368) * t214 + (-0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(5,1) - Ifges(6,3)) * t294) * t281) * t282 + (t169 / 0.2e1 + t170 / 0.2e1 + t168 / 0.2e1 + t229) * t214 + t95 * (Ifges(6,1) * t158 - Ifges(6,4) * t157) / 0.2e1 + t136 * (Ifges(6,1) * t109 - Ifges(6,4) * t110) / 0.2e1 + t7 * t329 + (Ifges(7,5) * t53 + Ifges(7,6) * t54 + Ifges(7,3) * t110) * t330 + (Ifges(7,5) * t112 + Ifges(7,6) * t111 + Ifges(7,3) * t157) * t336 + (Ifges(6,4) * t158 - Ifges(6,2) * t157) * t337 + (Ifges(7,1) * t53 + Ifges(7,4) * t54 + Ifges(7,5) * t110) * t338 + (Ifges(7,4) * t53 + Ifges(7,2) * t54 + Ifges(7,6) * t110) * t340 + t157 * t342 + (Ifges(7,4) * t112 + Ifges(7,2) * t111 + Ifges(7,6) * t157) * t343 + ((t86 * mrSges(5,2) - t106 * mrSges(4,3) - t87 * mrSges(5,3) + t140 * t369) * t218 + (-t106 * mrSges(4,1) + t86 * mrSges(5,1) + t108 * mrSges(4,2) + t139 * mrSges(3,3) - t71 * mrSges(5,3) + t351) * t221 + (t218 * t348 + t221 * t347) * qJD(2)) * t213 + t172 * (Ifges(6,5) * t109 - Ifges(6,6) * t110) / 0.2e1 + m(6) * (-t102 * t67 + t12 * t355 + t13 * t48 + t17 * t34 + t18 * t33 - t87 * t94) + t355 * t74 - t138 * (Ifges(6,4) * t109 - Ifges(6,2) * t110) / 0.2e1 + m(3) * (t139 * t283 - t140 * t159 - t153 * t352 + t154 * t155) + t289 * t352 + m(4) * (t106 * t142 + t108 * t141 + t113 * t352 + t122 * t134 + t126 * t129 + t140 * t143) + t297 * t102 - t87 * (mrSges(6,1) * t157 + mrSges(6,2) * t158) + t158 * t46 / 0.2e1 + t1 * (-mrSges(7,2) * t157 + mrSges(7,3) * t111) + t2 * (mrSges(7,1) * t157 - mrSges(7,3) * t112) + t155 * t148 + t88 * t147 + t134 * t149 + t126 * t150 + t101 * t151 - t110 * t63 / 0.2e1 + t111 * t8 / 0.2e1 + t6 * (-mrSges(7,1) * t111 + mrSges(7,2) * t112) + t112 * t9 / 0.2e1 + t109 * t64 / 0.2e1 + t110 * t30 / 0.2e1 + t10 * (mrSges(7,1) * t110 - mrSges(7,3) * t53) + t11 * (-mrSges(7,2) * t110 + mrSges(7,3) * t54) + t67 * (mrSges(6,1) * t110 + mrSges(6,2) * t109) + t17 * t98 + t18 * t99 + t94 * t52 + (Ifges(7,1) * t112 + Ifges(7,4) * t111 + Ifges(7,5) * t157) * t344 + t54 * t345 + (-t109 * t33 - t110 * t34 - t12 * t157 - t13 * t158) * mrSges(6,3) + m(7) * (t1 * t20 + t10 * t4 + t11 * t3 + t15 * t28 + t19 * t2 + t41 * t6) + m(5) * (t101 * t103 + t102 * t79 + t105 * t87 + t114 * t71 + t115 * t86 + t88 * t97) + t115 * t167 + t19 * t25 + t20 * t26 + t41 * t16 + t15 * t47 + t53 * t32 / 0.2e1 + t28 * (-mrSges(7,1) * t54 + mrSges(7,2) * t53) + t3 * t56 + t4 * t57 + t48 * t73; t229 - m(6) * (-t125 * t67 + t33 * t50 + t34 * t51) + (-t148 - t149) * t153 + (-t46 / 0.2e1 - t39 * t256 / 0.2e1 - t40 * t254 / 0.2e1 + t318 / 0.2e1 - t320 / 0.2e1 + t13 * mrSges(6,3) + t87 * mrSges(6,2) - t6 * t257 + t252 * t337 - t219 * t9 / 0.2e1 + t8 * t328 + t315 * t215 + t314 * qJD(3) + t260 * mrSges(7,3) + m(6) * (-qJD(3) * t33 - t13 * t215) + (mrSges(7,3) * t353 + t251 * t330 + t253 * t340 + t255 * t338 + t258 * t28 + t31 * t327 + t32 * t328) * qJD(6)) * t217 + (t10 * t131 + t11 * t130) * mrSges(7,3) + (((Ifges(6,5) * t217 / 0.2e1 + t220 * t362 - t317 * qJ(3) - t368) * qJD(2) + (t313 / 0.2e1 + t334) * t282 - t348 + t350) * t218 + (-t185 / 0.2e1 + (t333 + (Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t221) * t282 + (-mrSges(4,2) * pkin(2) - mrSges(5,3) * t222 - Ifges(5,5)) * qJD(2) + (-Ifges(4,1) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1) * t276 + t228 * t220 + t225 - t347) * t221) * t282 + t169 + t170 + t168 + (-pkin(2) * t140 + qJ(3) * t108 - t113 * t154 + t122 * t364 - t129 * t152) * m(4) + (Ifges(7,5) * t131 - Ifges(7,6) * t130) * t331 + (Ifges(7,1) * t131 - Ifges(7,4) * t130) * t339 + (Ifges(7,4) * t131 - Ifges(7,2) * t130) * t341 + t366 * t56 + (-t28 * t43 + t1 * t128 + t127 * t2 + (qJD(3) * t28 + t215 * t6) * t217 + t366 * t11 + t367 * t10) * m(7) + t367 * t57 + (t225 + (-m(6) * t34 - t98) * t215 * t217 + (t252 * t331 + t254 * t341 + t256 * t339 + t228 + (-m(6) * t33 + t314 + t335) * t215 + t372) * t220) * qJD(5) - t297 * t125 - t289 * t154 + t288 * qJD(3) + (-m(6) * t87 + t52) * t212 - t123 * t147 - t124 * t151 - t152 * t150 - t131 * t32 / 0.2e1 + t127 * t25 + t128 * t26 - t28 * (mrSges(7,1) * t130 + mrSges(7,2) * t131) - t51 * t98 - t50 * t99 + t130 * t345 + (t71 * qJ(3) - t103 * t124 - t125 * t79 + t222 * t87 + (-t123 + qJD(3)) * t97) * m(5) + (t38 / 0.2e1 + t37 / 0.2e1 + t357 - t87 * mrSges(6,1) + t7 / 0.2e1 + t342 - t12 * mrSges(6,3) + qJD(3) * t98 + t215 * t74 + m(6) * (qJD(3) * t34 + t12 * t215) + (t346 + Ifges(7,3) / 0.2e1) * t96 + t264) * t220 - t43 * t47; -t216 * t26 - t219 * t25 - t288 * t198 + t314 * t136 + t247 * qJD(6) + t241 * t138 + ((t150 - t151) * t218 + t317 * t280) * t282 - t52 + (-t133 * t353 + t136 * t28 - t260) * m(7) + (-t136 * t33 - t138 * t34 + t87) * m(6) + (-t122 * t198 + t129 * t276 + t140) * m(4) + (-t103 * t276 - t198 * t97 + t87) * m(5); -t130 * t57 + t131 * t56 + t167 - m(7) * (t10 * t130 - t11 * t131) + m(5) * t86 + (-m(7) * t6 + m(6) * (t13 + t296) - t315 + (m(7) * t353 - t241) * qJD(5)) * t220 + (t74 + t314 * qJD(5) + m(7) * (qJD(5) * t28 + t261) + m(6) * (-qJD(5) * t33 + t12) + t363) * t217 + ((-mrSges(5,1) * qJD(2) + t297) * t218 + (t217 * t314 + t220 * t98 + t147) * t221 - m(6) * (t218 * t67 - t290 * t34 + t292 * t33) - m(5) * (-t218 * t79 - t221 * t97) + t292 * t335) * t282; t351 + t226 * t136 + (-pkin(5) * t6 - t10 * t21 - t11 * t22 - t28 * t34) * m(7) + t8 * t327 + t9 * t328 + t251 * t336 + t365 * qJD(6) + (m(7) * t261 + t363) * pkin(10) - ((-Ifges(6,1) / 0.2e1 + t346) * t136 + t354 - t365) * t138 - t314 * t34 - t33 * t98 + t253 * t343 + t255 * t344 + t6 * t258 + t261 * mrSges(7,3) - Ifges(6,3) * t269 - pkin(5) * t16 - t22 * t56 - t21 * t57; -t28 * (mrSges(7,1) * t90 + mrSges(7,2) * t89) + (Ifges(7,1) * t89 - t326) * t339 + t31 * t338 + (Ifges(7,5) * t89 - Ifges(7,6) * t90) * t331 - t10 * t56 + t11 * t57 + (t10 * t89 + t11 * t90) * mrSges(7,3) + t264 + t7 + (-Ifges(7,2) * t90 + t32 + t85) * t341;];
tauc  = t23(:);
