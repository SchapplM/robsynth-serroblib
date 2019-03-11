% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPPRR7
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:17:11
% EndTime: 2019-03-09 09:17:36
% DurationCPUTime: 12.10s
% Computational Cost: add. (7379->645), mult. (18629->864), div. (0->0), fcn. (12741->8), ass. (0->282)
t217 = sin(qJ(2));
t212 = sin(pkin(6));
t296 = qJD(1) * t212;
t282 = t217 * t296;
t177 = qJD(5) + t282;
t181 = qJ(4) * t282;
t220 = cos(qJ(2));
t213 = cos(pkin(6));
t291 = t213 * qJD(1);
t289 = pkin(1) * t291;
t193 = t220 * t289;
t297 = pkin(8) * t282 - t193;
t123 = -t181 + t297;
t385 = -qJD(3) - t123;
t216 = sin(qJ(5));
t219 = cos(qJ(5));
t384 = -t385 - t177 * (pkin(5) * t216 - pkin(10) * t219);
t221 = -pkin(2) - pkin(3);
t211 = -pkin(9) + t221;
t281 = t220 * t296;
t293 = qJD(5) * t216;
t154 = pkin(8) * t281 + t217 * t289;
t125 = -qJ(4) * t281 + t154;
t183 = qJ(3) * t281;
t232 = t212 * (pkin(4) * t220 + t211 * t217);
t81 = qJD(1) * t232 + t183;
t51 = t219 * t125 + t216 * t81;
t383 = pkin(10) * t281 + t211 * t293 + t51;
t215 = sin(qJ(6));
t218 = cos(qJ(6));
t198 = qJD(2) + t291;
t276 = qJD(3) + t297;
t269 = -t181 + t276;
t67 = t198 * t211 + t269;
t129 = -pkin(1) * t296 - pkin(2) * t281 - qJ(3) * t282;
t103 = pkin(3) * t281 + qJD(4) - t129;
t238 = (pkin(4) * t217 + pkin(9) * t220) * t212;
t69 = qJD(1) * t238 + t103;
t31 = t216 * t69 + t219 * t67;
t29 = pkin(10) * t177 + t31;
t137 = -t198 * t219 + t216 * t281;
t138 = t216 * t198 + t219 * t281;
t180 = t198 * qJ(3);
t98 = -t125 - t180;
t80 = pkin(4) * t198 - t98;
t42 = -pkin(5) * t137 + pkin(10) * t138 + t80;
t12 = -t215 * t29 + t218 * t42;
t13 = t215 * t42 + t218 * t29;
t250 = t12 * t218 + t13 * t215;
t257 = mrSges(7,1) * t215 + mrSges(7,2) * t218;
t30 = -t216 * t67 + t219 * t69;
t28 = -pkin(5) * t177 - t30;
t135 = qJD(6) - t137;
t241 = t138 * t218 - t177 * t215;
t337 = Ifges(7,4) * t241;
t89 = t138 * t215 + t177 * t218;
t33 = Ifges(7,2) * t89 + Ifges(7,6) * t135 - t337;
t338 = t218 / 0.2e1;
t339 = t215 / 0.2e1;
t84 = Ifges(7,4) * t89;
t34 = -Ifges(7,1) * t241 + Ifges(7,5) * t135 + t84;
t382 = t250 * mrSges(7,3) - t28 * t257 + t33 * t339 - t34 * t338;
t381 = -Ifges(6,6) / 0.2e1;
t352 = t89 / 0.2e1;
t342 = t135 / 0.2e1;
t350 = -t241 / 0.2e1;
t380 = -mrSges(4,1) - mrSges(3,1);
t379 = mrSges(4,2) + mrSges(3,3);
t378 = Ifges(5,5) + Ifges(3,6);
t214 = qJ(3) + pkin(4);
t164 = pkin(5) * t219 + pkin(10) * t216 + t214;
t306 = t211 * t219;
t127 = t164 * t218 - t215 * t306;
t377 = qJD(6) * t127 + t215 * t384 - t383 * t218;
t302 = t218 * t219;
t128 = t164 * t215 + t211 * t302;
t376 = -qJD(6) * t128 + t383 * t215 + t218 * t384;
t252 = Ifges(7,5) * t218 - Ifges(7,6) * t215;
t321 = Ifges(7,4) * t218;
t254 = -Ifges(7,2) * t215 + t321;
t322 = Ifges(7,4) * t215;
t256 = Ifges(7,1) * t218 - t322;
t375 = t252 * t342 + t254 * t352 + t256 * t350 - t382;
t334 = pkin(1) * t220;
t290 = t213 * t334;
t194 = qJD(2) * t290;
t206 = t213 * qJD(3);
t295 = qJD(2) * t217;
t280 = t212 * t295;
t369 = (pkin(8) * t295 + qJD(4) * t220) * t212;
t87 = -qJ(4) * t280 - t194 - t206 + t369;
t176 = qJD(2) * t193;
t178 = t198 * qJD(3);
t275 = qJD(2) * t296;
t268 = t217 * t275;
t70 = -qJ(4) * t268 + qJD(1) * t369 - t176 - t178;
t278 = t219 * t295;
t292 = qJD(5) * t219;
t96 = -t198 * t292 + (t220 * t293 + t278) * t296;
t304 = t212 * t220;
t286 = t219 * t304;
t270 = qJD(5) * t286;
t97 = qJD(1) * t270 + t198 * t293 - t216 * t268;
t27 = -pkin(5) * t97 - pkin(10) * t96 - t70;
t230 = qJD(2) * t232;
t305 = t212 * t217;
t195 = qJD(3) * t305;
t267 = t220 * t275;
t300 = qJ(3) * t267 + qJD(1) * t195;
t59 = qJD(1) * t230 + t300;
t203 = t213 * t217 * pkin(1);
t102 = -qJD(4) * t305 + (t203 + (pkin(8) - qJ(4)) * t304) * qJD(2);
t86 = t102 * qJD(1);
t10 = t216 * t59 + t219 * t86 + t69 * t292 - t293 * t67;
t5 = pkin(10) * t267 + t10;
t1 = qJD(6) * t12 + t215 * t27 + t218 * t5;
t2 = -qJD(6) * t13 - t215 * t5 + t218 * t27;
t259 = t1 * t218 - t2 * t215;
t374 = qJD(5) * t28 - qJD(6) * t250 + t259;
t372 = t198 / 0.2e1;
t371 = -t296 / 0.2e1;
t370 = t296 / 0.2e1;
t134 = Ifges(6,4) * t137;
t320 = Ifges(6,5) * t177;
t324 = Ifges(6,1) * t138;
t62 = t134 + t320 - t324;
t368 = t30 * mrSges(6,3) - t62 / 0.2e1 - t320 / 0.2e1 - t80 * mrSges(6,2) - t134 / 0.2e1;
t307 = qJD(5) * t31;
t11 = -t216 * t86 + t219 * t59 - t307;
t367 = -t12 * t215 + t13 * t218;
t366 = t11 * mrSges(6,1) - t10 * mrSges(6,2) + Ifges(6,5) * t96 + Ifges(6,6) * t97;
t38 = qJD(6) * t89 + t215 * t267 + t218 * t96;
t39 = qJD(6) * t241 - t215 * t96 + t218 * t267;
t14 = -mrSges(7,1) * t39 + mrSges(7,2) * t38;
t72 = mrSges(6,1) * t267 - mrSges(6,3) * t96;
t326 = t72 - t14;
t6 = -pkin(5) * t267 - t11;
t365 = -m(7) * t6 + t326;
t364 = t177 * t381 + Ifges(6,4) * t138 / 0.2e1;
t113 = -pkin(2) * t198 + t276;
t187 = Ifges(3,4) * t281;
t288 = Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t78 = t198 * t221 + t269;
t362 = (Ifges(5,6) / 0.2e1 + t288) * t198 + t103 * mrSges(5,1) + t113 * mrSges(4,2) + t30 * mrSges(6,1) + (-Ifges(5,4) * t220 - Ifges(5,2) * t217) * t371 + (Ifges(4,1) * t217 - Ifges(4,5) * t220) * t370 + Ifges(3,1) * t282 / 0.2e1 + t187 / 0.2e1 + t177 * Ifges(6,3) - t138 * Ifges(6,5) + t137 * Ifges(6,6) - t129 * mrSges(4,3) + t297 * mrSges(3,3) - t31 * mrSges(6,2) - t78 * mrSges(5,3) + (Ifges(5,6) + Ifges(4,4) + Ifges(3,5)) * t372;
t122 = t180 + t154;
t186 = Ifges(4,5) * t282;
t359 = Ifges(4,6) / 0.2e1;
t361 = t103 * mrSges(5,2) + t129 * mrSges(4,1) + Ifges(4,6) * t372 - Ifges(4,3) * t281 / 0.2e1 + t186 / 0.2e1 + (Ifges(3,4) * t217 + t220 * Ifges(3,2)) * t371 + (-t220 * Ifges(5,1) - Ifges(5,4) * t217) * t370 - t122 * mrSges(4,2) - t154 * mrSges(3,3) - t98 * mrSges(5,3) + (t359 - Ifges(3,6) / 0.2e1 - Ifges(5,5) / 0.2e1 - t378 / 0.2e1) * t198;
t142 = -t212 * pkin(1) - pkin(2) * t304 - qJ(3) * t305;
t115 = pkin(3) * t304 - t142;
t83 = t238 + t115;
t199 = pkin(8) * t305;
t260 = -qJ(4) * t305 + t199;
t283 = -pkin(2) - t334;
t274 = -pkin(3) + t283;
t95 = (-pkin(9) + t274) * t213 + t260;
t325 = t216 * t83 + t219 * t95;
t294 = qJD(2) * t220;
t279 = t212 * t294;
t298 = qJ(3) * t279 + t195;
t68 = t230 + t298;
t18 = -qJD(5) * t325 - t102 * t216 + t219 * t68;
t358 = Ifges(7,5) * t350 + Ifges(7,6) * t352 + Ifges(7,3) * t342;
t360 = -Ifges(6,2) / 0.2e1;
t357 = t33 / 0.2e1;
t356 = t38 / 0.2e1;
t355 = t39 / 0.2e1;
t319 = Ifges(6,2) * t137;
t354 = -t319 / 0.2e1 + t364;
t353 = -t89 / 0.2e1;
t351 = t241 / 0.2e1;
t349 = -t97 / 0.2e1;
t348 = t97 / 0.2e1;
t347 = m(7) * t28;
t346 = pkin(1) * mrSges(3,1);
t345 = pkin(1) * mrSges(3,2);
t343 = -t135 / 0.2e1;
t157 = -t213 * t219 + t216 * t304;
t341 = t157 / 0.2e1;
t158 = -t213 * t216 - t286;
t340 = t158 / 0.2e1;
t37 = Ifges(7,5) * t38;
t36 = Ifges(7,6) * t39;
t331 = t96 * Ifges(6,1);
t330 = t96 * Ifges(6,4);
t329 = t97 * Ifges(6,4);
t328 = mrSges(4,2) - mrSges(5,3);
t100 = mrSges(6,1) * t177 + mrSges(6,3) * t138;
t49 = -mrSges(7,1) * t89 - mrSges(7,2) * t241;
t309 = t100 - t49;
t308 = -mrSges(5,1) * t198 + mrSges(6,1) * t137 + mrSges(6,2) * t138 + mrSges(5,3) * t281;
t303 = t216 * t217;
t301 = t198 * t380 + t282 * t379;
t299 = mrSges(5,1) * t267 + mrSges(5,2) * t268;
t160 = pkin(8) * t304 + t203;
t7 = -Ifges(7,3) * t97 + t36 + t37;
t149 = mrSges(4,2) * t281 + mrSges(4,3) * t198;
t287 = t149 - t308;
t141 = t213 * qJ(3) + t160;
t273 = t221 * t305;
t272 = t216 * t282;
t271 = t219 * t282;
t265 = 0.3e1 / 0.2e1 * Ifges(5,4) - 0.3e1 / 0.2e1 * Ifges(4,5) + 0.3e1 / 0.2e1 * Ifges(3,4);
t263 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t258 = -mrSges(7,1) * t218 + mrSges(7,2) * t215;
t255 = Ifges(7,1) * t215 + t321;
t253 = Ifges(7,2) * t218 + t322;
t251 = Ifges(7,5) * t215 + Ifges(7,6) * t218;
t25 = -mrSges(7,1) * t97 - mrSges(7,3) * t38;
t26 = mrSges(7,2) * t97 + mrSges(7,3) * t39;
t248 = -t215 * t25 + t218 * t26;
t41 = pkin(10) * t305 + t325;
t114 = qJ(4) * t304 - t141;
t104 = t213 * pkin(4) - t114;
t55 = -pkin(5) * t157 - pkin(10) * t158 + t104;
t22 = t215 * t55 + t218 * t41;
t21 = -t215 * t41 + t218 * t55;
t56 = -mrSges(7,2) * t135 + mrSges(7,3) * t89;
t57 = mrSges(7,1) * t135 + mrSges(7,3) * t241;
t247 = -t215 * t56 - t218 * t57;
t246 = t216 * t31 + t219 * t30;
t47 = -t216 * t95 + t219 * t83;
t50 = -t125 * t216 + t219 * t81;
t242 = qJD(2) * t273;
t155 = -pkin(8) * t280 + t194;
t99 = -mrSges(6,2) * t177 + mrSges(6,3) * t137;
t240 = t215 * t57 - t218 * t56 - t99;
t111 = -t158 * t215 + t218 * t305;
t112 = t158 * t218 + t215 * t305;
t17 = t219 * t102 + t216 * t68 + t83 * t292 - t293 * t95;
t139 = -pkin(8) * t268 + t176;
t156 = t160 * qJD(2);
t108 = t139 + t178;
t140 = qJD(1) * t156;
t229 = -t70 * mrSges(5,1) - t139 * mrSges(3,2) + t86 * mrSges(5,2) + t108 * mrSges(4,3) + t140 * t380;
t228 = t324 / 0.2e1 + t368;
t73 = -mrSges(6,2) * t267 + mrSges(6,3) * t97;
t226 = -qJD(5) * t309 + t247 * qJD(6) + t248 + t73;
t225 = t80 * mrSges(6,1) + t12 * mrSges(7,1) - t13 * mrSges(7,2) - t31 * mrSges(6,3) + t354 + 0.2e1 * t358 + t364;
t224 = (t319 / 0.2e1 - t225) * t216;
t175 = Ifges(4,4) * t267;
t174 = Ifges(3,5) * t267;
t173 = Ifges(4,6) * t268;
t172 = Ifges(6,3) * t267;
t159 = -t199 + t290;
t152 = (mrSges(5,1) * t217 - mrSges(5,2) * t220) * t296;
t151 = pkin(2) * t282 - t183;
t150 = (-mrSges(4,1) * t220 - mrSges(4,3) * t217) * t296;
t148 = -mrSges(3,2) * t198 + mrSges(3,3) * t281;
t144 = mrSges(5,2) * t198 - mrSges(5,3) * t282;
t143 = t213 * t283 + t199;
t136 = t155 + t206;
t133 = (t215 * t220 + t217 * t302) * t296;
t132 = t215 * t271 - t218 * t281;
t131 = t198 * t215 + t218 * t272;
t130 = t198 * t218 - t215 * t272;
t126 = pkin(2) * t280 - t298;
t124 = qJD(1) * t273 + t183;
t110 = qJD(5) * t157 + t212 * t278;
t109 = -t213 * t293 + t216 * t280 - t270;
t106 = pkin(2) * t268 - t300;
t105 = t213 * t274 + t260;
t101 = t242 + t298;
t85 = qJD(1) * t242 + t300;
t76 = -pkin(5) * t138 - pkin(10) * t137;
t54 = qJD(6) * t111 + t110 * t218 + t215 * t279;
t53 = -qJD(6) * t112 - t110 * t215 + t218 * t279;
t52 = -mrSges(6,1) * t97 + mrSges(6,2) * t96;
t46 = Ifges(6,5) * t267 + t329 + t331;
t45 = t97 * Ifges(6,2) + Ifges(6,6) * t267 + t330;
t43 = -pkin(5) * t281 - t50;
t40 = -pkin(5) * t305 - t47;
t35 = pkin(5) * t109 - pkin(10) * t110 - t87;
t20 = t215 * t76 + t218 * t30;
t19 = -t215 * t30 + t218 * t76;
t16 = -pkin(5) * t279 - t18;
t15 = pkin(10) * t279 + t17;
t9 = Ifges(7,1) * t38 + Ifges(7,4) * t39 - Ifges(7,5) * t97;
t8 = Ifges(7,4) * t38 + Ifges(7,2) * t39 - Ifges(7,6) * t97;
t4 = -qJD(6) * t22 - t15 * t215 + t218 * t35;
t3 = qJD(6) * t21 + t15 * t218 + t215 * t35;
t23 = [(t10 * t157 - t109 * t31 - t11 * t158 - t110 * t30) * mrSges(6,3) + m(7) * (t1 * t22 + t12 * t4 + t13 * t3 + t16 * t28 + t2 * t21 + t40 * t6) + m(5) * (t101 * t103 + t102 * t78 + t105 * t86 + t114 * t70 + t115 * t85 + t87 * t98) + m(4) * (t106 * t142 + t108 * t141 + t113 * t156 + t122 * t136 + t126 * t129 + t140 * t143) + ((t142 * mrSges(4,1) - t141 * mrSges(4,2) - t160 * mrSges(3,3) - t114 * mrSges(5,3) + (-t217 * t265 - 0.2e1 * t346) * t212 + (t359 - t378) * t213) * t295 + (t143 * mrSges(4,2) - t159 * mrSges(3,3) - t105 * mrSges(5,3) + Ifges(6,5) * t340 + Ifges(6,6) * t341 - t142 * mrSges(4,3) + (t220 * t265 - 0.2e1 * t345) * t212 + (Ifges(5,6) + t288) * t213 + (-0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + Ifges(6,3) / 0.2e1) * t305) * t294) * t296 + m(6) * (t10 * t325 - t104 * t70 + t11 * t47 + t17 * t31 + t18 * t30 - t80 * t87) + t325 * t73 + (t174 / 0.2e1 + t175 / 0.2e1 + t173 / 0.2e1 + t229) * t213 + t137 * (Ifges(6,4) * t110 - Ifges(6,2) * t109) / 0.2e1 + t177 * (Ifges(6,5) * t110 - Ifges(6,6) * t109) / 0.2e1 + t308 * t87 + t115 * t299 + t301 * t156 + t155 * t148 + t2 * (-mrSges(7,1) * t157 - mrSges(7,3) * t112) + t1 * (mrSges(7,2) * t157 + mrSges(7,3) * t111) - t157 * t7 / 0.2e1 - t70 * (-mrSges(6,1) * t157 + mrSges(6,2) * t158) + t136 * t149 + t126 * t150 + t101 * t152 + t102 * t144 + t6 * (-mrSges(7,1) * t111 + mrSges(7,2) * t112) + t112 * t9 / 0.2e1 + t111 * t8 / 0.2e1 + t13 * (-mrSges(7,2) * t109 + mrSges(7,3) * t53) + t12 * (mrSges(7,1) * t109 - mrSges(7,3) * t54) + t80 * (mrSges(6,1) * t109 + mrSges(6,2) * t110) + t110 * t62 / 0.2e1 + t104 * t52 + t17 * t99 + t18 * t100 + t47 * t72 + t3 * t56 + t4 * t57 - t138 * (Ifges(6,1) * t110 - Ifges(6,4) * t109) / 0.2e1 + t96 * (Ifges(6,1) * t158 + Ifges(6,4) * t157) / 0.2e1 + t46 * t340 + t45 * t341 + (Ifges(7,5) * t54 + Ifges(7,6) * t53 + Ifges(7,3) * t109) * t342 + (Ifges(6,4) * t158 + Ifges(6,2) * t157) * t348 + (Ifges(7,5) * t112 + Ifges(7,6) * t111 - Ifges(7,3) * t157) * t349 + (Ifges(7,1) * t54 + Ifges(7,4) * t53 + Ifges(7,5) * t109) * t350 + (Ifges(7,4) * t54 + Ifges(7,2) * t53 + Ifges(7,6) * t109) * t352 + t109 * t354 + (Ifges(7,4) * t112 + Ifges(7,2) * t111 - Ifges(7,6) * t157) * t355 + (Ifges(7,1) * t112 + Ifges(7,4) * t111 - Ifges(7,5) * t157) * t356 + t53 * t357 + t109 * t358 + ((-mrSges(4,1) * t106 + mrSges(4,2) * t108 - mrSges(5,2) * t85 + mrSges(3,3) * t139 + mrSges(5,3) * t70) * t220 + (t172 / 0.2e1 - t106 * mrSges(4,3) + t85 * mrSges(5,1) - t86 * mrSges(5,3) + t379 * t140 + t366) * t217 + (t217 * t361 + t220 * t362) * qJD(2)) * t212 + t21 * t25 + t22 * t26 + m(3) * (t139 * t160 - t140 * t159 + t154 * t155 + t156 * t297) + t40 * t14 + t16 * t49 + t54 * t34 / 0.2e1 + t28 * (-mrSges(7,1) * t53 + mrSges(7,2) * t54); (-t218 * t9 / 0.2e1 + t8 * t339 - t6 * t257 - t38 * t256 / 0.2e1 - t39 * t254 / 0.2e1 + t252 * t348 - t46 / 0.2e1 - t331 / 0.2e1 - t329 / 0.2e1 + t11 * mrSges(6,3) + t70 * mrSges(6,2) + (t1 * t215 + t2 * t218) * mrSges(7,3) + (-m(6) * t11 - t365) * t211 + (mrSges(7,3) * t367 + t251 * t342 + t253 * t352 + t255 * t350 + t258 * t28 + t33 * t338 + t339 * t34) * qJD(6)) * t216 + ((-t187 / 0.2e1 + (Ifges(5,6) - t221 * mrSges(5,3) - pkin(2) * mrSges(4,2) - Ifges(6,5) * t216 / 0.2e1 + t219 * t381) * qJD(2) + (t345 + (Ifges(4,5) / 0.2e1 - Ifges(5,4) / 0.2e1) * t220) * t296 - t362) * t220 + (-t186 / 0.2e1 + (t346 + (Ifges(3,4) / 0.2e1 + Ifges(5,4) / 0.2e1) * t217) * t296 + (Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1 + Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1 - Ifges(4,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t281 + (-qJ(3) * t328 - t378) * qJD(2) + t228 * t219 + t224 - t361) * t217) * t296 + t173 + t174 + t175 + (t37 / 0.2e1 + t36 / 0.2e1 - t70 * mrSges(6,1) - t330 / 0.2e1 + t7 / 0.2e1 - t45 / 0.2e1 - t10 * mrSges(6,3) + t211 * t73 + (-Ifges(7,3) / 0.2e1 + t360) * t97 + t263) * t219 - m(6) * (-t123 * t80 + t30 * t50 + t31 * t51) + t229 + (t224 + (-m(6) * t246 - t216 * t99) * t211 + (t252 * t343 + t254 * t353 + t256 * t351 + t228 + (-t309 + t347) * t211 + t382) * t219) * qJD(5) - t308 * t123 + m(6) * (qJD(3) * t80 + t10 * t306 - t214 * t70) - t301 * t154 + t287 * qJD(3) + (-t70 * qJ(3) - t103 * t124 - t125 * t78 + t221 * t86 + t385 * t98) * m(5) + t214 * t52 - t151 * t150 - t124 * t152 - t125 * t144 - t28 * (mrSges(7,1) * t132 + mrSges(7,2) * t133) - t133 * t34 / 0.2e1 + t127 * t25 + t128 * t26 - t51 * t99 - t50 * t100 + (t12 * t133 + t13 * t132) * mrSges(7,3) + (Ifges(7,5) * t133 - Ifges(7,6) * t132) * t343 + (Ifges(7,1) * t133 - Ifges(7,4) * t132) * t351 + (Ifges(7,4) * t133 - Ifges(7,2) * t132) * t353 + t132 * t357 + (-pkin(2) * t140 + qJ(3) * t108 - t113 * t154 + t122 * t276 - t129 * t151) * m(4) + t376 * t57 + t377 * t56 + (t1 * t128 + t12 * t376 + t127 * t2 + t13 * t377 - t28 * t43) * m(7) - (-t149 - t148) * t297 - t43 * t49; -t130 * t57 - t131 * t56 - t287 * t198 + ((t150 - t152) * t217 + t328 * t294) * t296 + (t240 * qJD(5) - t99 * t282 - t326) * t216 + (-t309 * t282 + t226) * t219 + ((-qJD(5) * t367 + t6) * t216 + t374 * t219 - t12 * t130 - t13 * t131 + t28 * t271) * m(7) + (t10 * t219 - t11 * t216 - t177 * t246 - t198 * t80) * m(6) + (-t103 * t282 + t198 * t98 + t86) * m(5) + (-t122 * t198 + t129 * t282 + t140) * m(4); -t132 * t57 + t133 * t56 - m(7) * (t12 * t132 - t13 * t133) + m(5) * t85 + (m(6) * (t11 + t307) + (m(7) * t367 - t240) * qJD(5) + t365) * t219 + (m(7) * t374 + m(6) * (-qJD(5) * t30 + t10) + t226) * t216 + (-t308 * t220 + (-t216 * t309 + t219 * t99 + t144) * t217 - m(6) * (-t217 * t219 * t31 - t220 * t80 + t30 * t303) - m(5) * (-t217 * t78 + t220 * t98) + t303 * t347) * t296 + t299; t366 + t225 * t138 + t172 + ((t360 + Ifges(6,1) / 0.2e1) * t138 + t368 - t375) * t137 + t375 * qJD(6) + (m(7) * t259 + t248 + (-m(7) * t250 + t247) * qJD(6)) * pkin(10) + (-pkin(5) * t6 - t12 * t19 - t13 * t20 - t28 * t31) * m(7) + t259 * mrSges(7,3) + t6 * t258 + t309 * t31 - t30 * t99 - t20 * t56 - t19 * t57 + t8 * t338 + t9 * t339 + t251 * t349 + t253 * t355 + t255 * t356 - pkin(5) * t14; -t28 * (-mrSges(7,1) * t241 + mrSges(7,2) * t89) + (Ifges(7,1) * t89 + t337) * t351 + t33 * t350 + (Ifges(7,5) * t89 + Ifges(7,6) * t241) * t343 - t12 * t56 + t13 * t57 + (t12 * t89 - t13 * t241) * mrSges(7,3) + t263 + t7 + (Ifges(7,2) * t241 + t34 + t84) * t353;];
tauc  = t23(:);
