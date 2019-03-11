% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:16:47
% EndTime: 2019-03-09 12:17:05
% DurationCPUTime: 9.14s
% Computational Cost: add. (7639->576), mult. (17592->720), div. (0->0), fcn. (10996->6), ass. (0->256)
t209 = sin(qJ(2));
t283 = qJD(1) * t209;
t194 = pkin(7) * t283;
t152 = pkin(8) * t283 - t194;
t383 = qJD(3) - t152;
t208 = sin(qJ(4));
t211 = cos(qJ(4));
t212 = cos(qJ(2));
t282 = qJD(1) * t212;
t129 = -t208 * t283 - t211 * t282;
t123 = qJD(5) - t129;
t382 = Ifges(6,1) + Ifges(7,1);
t374 = Ifges(7,4) + Ifges(6,5);
t381 = -Ifges(6,6) + Ifges(7,6);
t380 = -mrSges(3,1) - mrSges(4,1);
t195 = pkin(7) * t282;
t153 = -pkin(8) * t282 + t195;
t213 = -pkin(2) - pkin(3);
t158 = t211 * qJ(3) + t208 * t213;
t371 = qJD(4) * t158 + t211 * t153 + t208 * t383;
t207 = sin(qJ(5));
t210 = cos(qJ(5));
t240 = pkin(5) * t207 - qJ(6) * t210;
t379 = -t207 * qJD(6) + t123 * t240;
t138 = -t208 * t212 + t209 * t211;
t131 = -t208 * t282 + t211 * t283;
t273 = qJD(2) - qJD(4);
t102 = t207 * t131 + t210 * t273;
t103 = t210 * t131 - t207 * t273;
t167 = mrSges(7,1) * t210 + mrSges(7,3) * t207;
t301 = Ifges(7,5) * t207;
t169 = -Ifges(7,3) * t210 + t301;
t133 = -qJD(1) * pkin(1) - pkin(2) * t282 - qJ(3) * t283;
t110 = pkin(3) * t282 - t133;
t69 = -pkin(4) * t129 - pkin(9) * t131 + t110;
t262 = t213 * qJD(2);
t112 = t262 + t383;
t205 = qJD(2) * qJ(3);
t132 = t153 + t205;
t82 = t208 * t112 + t211 * t132;
t77 = -pkin(9) * t273 + t82;
t19 = -t207 * t77 + t210 * t69;
t347 = qJD(6) - t19;
t17 = -pkin(5) * t123 + t347;
t305 = Ifges(6,4) * t207;
t172 = Ifges(6,2) * t210 + t305;
t300 = Ifges(7,5) * t210;
t173 = Ifges(7,1) * t207 - t300;
t304 = Ifges(6,4) * t210;
t174 = Ifges(6,1) * t207 + t304;
t20 = t207 * t69 + t210 * t77;
t18 = qJ(6) * t123 + t20;
t242 = Ifges(7,3) * t207 + t300;
t250 = mrSges(6,1) * t210 - mrSges(6,2) * t207;
t306 = Ifges(6,4) * t103;
t49 = -t102 * Ifges(6,2) + t123 * Ifges(6,6) + t306;
t332 = -t49 / 0.2e1;
t342 = t382 * t210 + t301 - t305;
t343 = t381 * t207 + t374 * t210;
t245 = -Ifges(6,2) * t207 + t304;
t360 = -t245 / 0.2e1;
t248 = mrSges(7,1) * t207 - mrSges(7,3) * t210;
t249 = mrSges(6,1) * t207 + mrSges(6,2) * t210;
t81 = t211 * t112 - t208 * t132;
t76 = pkin(4) * t273 - t81;
t28 = t102 * pkin(5) - t103 * qJ(6) + t76;
t369 = t248 * t28 + t249 * t76;
t281 = qJD(2) * t209;
t329 = pkin(7) - pkin(8);
t154 = t329 * t281;
t204 = qJD(2) * qJD(3);
t116 = -qJD(1) * t154 + t204;
t176 = t329 * t212;
t155 = qJD(2) * t176;
t234 = qJD(1) * t155;
t350 = qJD(4) * t81;
t39 = t211 * t116 + t208 * t234 + t350;
t293 = qJD(4) * t82;
t40 = t116 * t208 - t211 * t234 + t293;
t100 = Ifges(7,5) * t103;
t46 = t123 * Ifges(7,6) + t102 * Ifges(7,3) + t100;
t302 = Ifges(7,5) * t102;
t50 = Ifges(7,1) * t103 + Ifges(7,4) * t123 + t302;
t101 = Ifges(6,4) * t102;
t51 = Ifges(6,1) * t103 + Ifges(6,5) * t123 - t101;
t137 = t208 * t209 + t211 * t212;
t97 = t273 * t137;
t88 = t97 * qJD(1);
t58 = -qJD(5) * t102 + t210 * t88;
t59 = qJD(5) * t103 + t207 * t88;
t7 = pkin(5) * t59 - qJ(6) * t58 - qJD(6) * t103 + t40;
t378 = ((t46 / 0.2e1 + t332 - t20 * mrSges(6,3) - t18 * mrSges(7,2)) * t207 + (t50 / 0.2e1 + t51 / 0.2e1 + t17 * mrSges(7,2) - t19 * mrSges(6,3)) * t210) * qJD(5) + (t174 / 0.2e1 + t173 / 0.2e1) * t58 + (t169 / 0.2e1 - t172 / 0.2e1) * t59 + (-t250 - mrSges(5,1)) * t40 - t39 * mrSges(5,2) - t7 * t167 + (t102 * t360 + t369) * qJD(5) + (t102 * t242 + t103 * t342 + t123 * t343) * qJD(5) / 0.2e1;
t377 = Ifges(5,1) / 0.2e1;
t376 = -Ifges(5,2) / 0.2e1;
t375 = -t282 / 0.2e1;
t373 = -t379 + t371;
t297 = -mrSges(5,1) * t273 - mrSges(6,1) * t102 - mrSges(6,2) * t103 - t131 * mrSges(5,3);
t157 = -t208 * qJ(3) + t211 * t213;
t113 = t211 * qJD(3) + qJD(4) * t157;
t95 = t152 * t211 + t153 * t208;
t372 = t113 - t95;
t370 = (t17 * t210 - t18 * t207) * mrSges(7,2) - (t19 * t210 + t20 * t207) * mrSges(6,3);
t255 = qJD(1) * t281;
t279 = qJD(2) * t212;
t344 = t138 * qJD(4) + t208 * t279;
t89 = qJD(1) * t344 - t211 * t255;
t21 = -mrSges(7,2) * t59 + mrSges(7,3) * t89;
t22 = mrSges(6,1) * t89 - mrSges(6,3) * t58;
t23 = -t89 * mrSges(7,1) + t58 * mrSges(7,2);
t24 = -mrSges(6,2) * t89 - mrSges(6,3) * t59;
t309 = mrSges(6,3) * t103;
t74 = mrSges(6,1) * t123 - t309;
t75 = -mrSges(7,1) * t123 + mrSges(7,2) * t103;
t312 = -t75 + t74;
t72 = -mrSges(7,2) * t102 + mrSges(7,3) * t123;
t310 = mrSges(6,3) * t102;
t73 = -mrSges(6,2) * t123 - t310;
t313 = t72 + t73;
t218 = (t21 + t24) * t210 + (-t22 + t23) * t207 + (-t207 * t313 - t210 * t312) * qJD(5);
t275 = qJD(5) * t210;
t276 = qJD(5) * t207;
t251 = t209 * t262;
t189 = qJ(3) * t282;
t198 = t209 * qJD(3);
t285 = qJD(1) * t198 + qJD(2) * t189;
t99 = qJD(1) * t251 + t285;
t29 = pkin(4) * t89 - pkin(9) * t88 + t99;
t3 = t207 * t29 + t210 * t39 + t69 * t275 - t276 * t77;
t4 = -qJD(5) * t20 - t207 * t39 + t210 * t29;
t221 = -t19 * t275 - t20 * t276 - t207 * t4 + t210 * t3;
t1 = qJ(6) * t89 + qJD(6) * t123 + t3;
t2 = -pkin(5) * t89 - t4;
t222 = t1 * t210 + t17 * t275 - t18 * t276 + t2 * t207;
t368 = m(6) * t221 + m(7) * t222 + t218;
t159 = -qJD(2) * pkin(2) + qJD(3) + t194;
t303 = Ifges(4,5) * t212;
t356 = -qJD(2) / 0.2e1;
t357 = -qJD(1) / 0.2e1;
t361 = Ifges(3,4) * t375;
t364 = -Ifges(3,1) / 0.2e1;
t367 = (m(4) * t159 + (mrSges(4,2) + mrSges(3,3)) * t283 + t380 * qJD(2)) * pkin(7) - t133 * mrSges(4,3) - (t209 * Ifges(4,1) - t303) * t357 - t283 * t364 - t361 + t159 * mrSges(4,2) - (Ifges(4,4) + Ifges(3,5)) * t356;
t164 = t195 + t205;
t166 = mrSges(4,2) * t282 + qJD(2) * mrSges(4,3);
t191 = Ifges(4,5) * t283;
t307 = Ifges(3,4) * t209;
t363 = Ifges(4,6) / 0.2e1;
t366 = -(m(4) * t164 - qJD(2) * mrSges(3,2) + mrSges(3,3) * t282 + t166) * pkin(7) + t133 * mrSges(4,1) + qJD(2) * t363 + Ifges(4,3) * t375 + t191 / 0.2e1 + Ifges(3,6) * t356 + (t212 * Ifges(3,2) + t307) * t357 - t164 * mrSges(4,2);
t326 = -t102 / 0.2e1;
t325 = t102 / 0.2e1;
t323 = t103 / 0.2e1;
t119 = Ifges(5,4) * t129;
t263 = t119 / 0.2e1;
t120 = Ifges(5,4) * t131;
t322 = -t120 / 0.2e1;
t362 = t123 / 0.2e1;
t359 = t129 * t376;
t358 = t131 * t377;
t354 = t50 + t51;
t351 = -t82 + t379;
t348 = (t376 + t377) * t131;
t175 = t329 * t209;
t346 = t211 * t175 - t176 * t208;
t265 = Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1;
t267 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t269 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t337 = t265 * t102 - t269 * t103 - t267 * t123 - t110 * mrSges(5,1) - t19 * mrSges(6,1) + t17 * mrSges(7,1) + t20 * mrSges(6,2) + t82 * mrSges(5,3) - t18 * mrSges(7,3) - Ifges(5,6) * t273 - Ifges(6,6) * t326 - Ifges(7,6) * t325 - t322 - t359 - (Ifges(7,2) + Ifges(6,3)) * t362 - t374 * t323;
t345 = -t337 + t322;
t280 = qJD(2) * t211;
t341 = t263 + t370;
t105 = t175 * t208 + t176 * t211;
t161 = -t212 * pkin(2) - t209 * qJ(3) - pkin(1);
t136 = t212 * pkin(3) - t161;
t80 = pkin(4) * t137 - pkin(9) * t138 + t136;
t311 = t210 * t105 + t207 * t80;
t284 = qJ(3) * t279 + t198;
t108 = t251 + t284;
t96 = -t209 * t280 + t344;
t38 = pkin(4) * t96 - pkin(9) * t97 + t108;
t63 = qJD(4) * t346 - t154 * t211 + t155 * t208;
t9 = -qJD(5) * t311 - t207 * t63 + t210 * t38;
t317 = t210 / 0.2e1;
t318 = t207 / 0.2e1;
t319 = -t207 / 0.2e1;
t338 = -t110 * mrSges(5,2) + t81 * mrSges(5,3) + Ifges(5,5) * t273 - t317 * t354 - t46 * t318 - t49 * t319 - t263 - t358 - t369;
t335 = m(5) / 0.2e1;
t334 = m(6) / 0.2e1;
t333 = m(7) / 0.2e1;
t328 = pkin(1) * mrSges(3,1);
t327 = pkin(1) * mrSges(3,2);
t324 = -t103 / 0.2e1;
t321 = -t123 / 0.2e1;
t316 = pkin(5) * t131;
t92 = pkin(4) * t131 - pkin(9) * t129;
t42 = t207 * t92 + t210 * t81;
t118 = t213 * t283 + t189;
t70 = t118 - t92;
t34 = t207 * t70 + t210 * t95;
t298 = t346 * t40;
t294 = qJ(6) * t131;
t292 = t113 * t207;
t291 = t113 * t210;
t278 = qJD(4) * t207;
t277 = qJD(4) * t210;
t271 = Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t270 = -0.3e1 / 0.2e1 * Ifges(4,5) + 0.3e1 / 0.2e1 * Ifges(3,4);
t266 = t363 - Ifges(3,6) / 0.2e1;
t264 = m(4) * pkin(7) + mrSges(4,2);
t170 = Ifges(6,5) * t207 + Ifges(6,6) * t210;
t171 = Ifges(7,4) * t207 - Ifges(7,6) * t210;
t252 = t170 / 0.2e1 + t171 / 0.2e1 - Ifges(5,6);
t241 = t210 * pkin(5) + t207 * qJ(6);
t33 = -t207 * t95 + t210 * t70;
t41 = -t207 * t81 + t210 * t92;
t44 = -t105 * t207 + t210 * t80;
t160 = -pkin(4) - t241;
t8 = -t105 * t276 + t207 * t38 + t210 * t63 + t80 * t275;
t11 = t58 * Ifges(7,5) + t89 * Ifges(7,6) + t59 * Ifges(7,3);
t12 = t58 * Ifges(6,4) - t59 * Ifges(6,2) + t89 * Ifges(6,6);
t231 = -t11 / 0.2e1 + t12 / 0.2e1 + t1 * mrSges(7,2) + t3 * mrSges(6,3);
t13 = t58 * Ifges(7,1) + t89 * Ifges(7,4) + t59 * Ifges(7,5);
t14 = t58 * Ifges(6,1) - t59 * Ifges(6,4) + t89 * Ifges(6,5);
t230 = t13 / 0.2e1 + t14 / 0.2e1 + t2 * mrSges(7,2) - t4 * mrSges(6,3);
t224 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t64 = qJD(4) * t105 - t154 * t208 - t211 * t155;
t214 = t242 * t326 + t245 * t325 + t321 * t343 + t324 * t342 + t338;
t156 = -pkin(7) * t255 + t204;
t150 = pkin(4) - t157;
t140 = (-t212 * mrSges(4,1) - mrSges(4,3) * t209) * qJD(1);
t130 = t207 * t283 + t210 * t280;
t128 = t207 * t280 - t210 * t283;
t122 = pkin(2) * t281 - t284;
t111 = -t157 - t160;
t109 = pkin(2) * t255 - t285;
t106 = mrSges(5,2) * t273 + t129 * mrSges(5,3);
t90 = -mrSges(5,1) * t129 + mrSges(5,2) * t131;
t87 = Ifges(7,2) * t89;
t86 = Ifges(6,3) * t89;
t66 = mrSges(7,1) * t102 - mrSges(7,3) * t103;
t65 = pkin(5) * t103 + qJ(6) * t102;
t61 = t138 * t240 - t346;
t57 = Ifges(7,4) * t58;
t56 = Ifges(6,5) * t58;
t55 = Ifges(6,6) * t59;
t54 = Ifges(7,6) * t59;
t35 = -pkin(5) * t137 - t44;
t32 = qJ(6) * t137 + t311;
t31 = -t41 - t316;
t30 = t294 + t42;
t26 = -t33 + t316;
t25 = -t294 + t34;
t16 = mrSges(6,1) * t59 + mrSges(6,2) * t58;
t15 = mrSges(7,1) * t59 - mrSges(7,3) * t58;
t10 = t240 * t97 + (qJD(5) * t241 - qJD(6) * t210) * t138 + t64;
t6 = -pkin(5) * t96 - t9;
t5 = qJ(6) * t96 + qJD(6) * t137 + t8;
t27 = [(Ifges(5,1) * t88 + t7 * t248 + t99 * mrSges(5,2) + t11 * t318 + t12 * t319 + (mrSges(5,3) + t249) * t40 + (-t207 * t3 - t210 * t4) * mrSges(6,3) + (-t1 * t207 + t2 * t210) * mrSges(7,2) + (t172 * t325 + t169 * t326 + t76 * t250 + t28 * t167 + t210 * t332 + (t19 * t207 - t20 * t210) * mrSges(6,3) + (-t17 * t207 - t18 * t210) * mrSges(7,2) + (t174 + t173) * t324 + (t171 + t170) * t321 + t354 * t319) * qJD(5) + (t242 / 0.2e1 + t360) * t59 + t342 * t58 / 0.2e1 + (qJD(5) * t46 + t13 + t14) * t317 + (-Ifges(5,4) + t343 / 0.2e1) * t89) * t138 + t311 * t24 + m(7) * (t1 * t32 + t10 * t28 + t17 * t6 + t18 * t5 + t2 * t35 + t61 * t7) + m(4) * (t109 * t161 + t122 * t133) + (-t109 * mrSges(4,3) + (t266 * qJD(2) + (t161 * mrSges(4,1) - t209 * t270 - 0.2e1 * t328) * qJD(1) + t366) * qJD(2)) * t209 + (-t109 * mrSges(4,1) + t264 * t156 + (t271 * qJD(2) + (-0.2e1 * t327 - t161 * mrSges(4,3) + t270 * t212 + (-0.3e1 / 0.2e1 * Ifges(4,3) - 0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1) + 0.3e1 / 0.2e1 * Ifges(4,1) + t264 * pkin(7)) * t209) * qJD(1) + t367) * qJD(2)) * t212 + (-t105 * t89 - t346 * t88) * mrSges(5,3) - t346 * t16 + m(6) * (t19 * t9 + t20 * t8 + t3 * t311 + t4 * t44 + t64 * t76 - t298) + m(5) * (t105 * t39 + t108 * t110 + t136 * t99 + t63 * t82 - t64 * t81 - t298) + t122 * t140 + t136 * (mrSges(5,1) * t89 + mrSges(5,2) * t88) + t63 * t106 + t108 * t90 + t5 * t72 + t8 * t73 + t9 * t74 + t6 * t75 + t61 * t15 + t10 * t66 + t44 * t22 + t35 * t23 + t32 * t21 + (-t214 + t358 + t341) * t97 + (t359 + t345) * t96 + (t99 * mrSges(5,1) - Ifges(5,4) * t88 + t56 / 0.2e1 - t55 / 0.2e1 + t86 / 0.2e1 + t57 / 0.2e1 + t87 / 0.2e1 + t54 / 0.2e1 - t39 * mrSges(5,3) - t265 * t59 + t269 * t58 + (Ifges(5,2) + t267) * t89 + t224) * t137 - t297 * t64; -t378 + ((t361 + (t327 + t303 / 0.2e1) * qJD(1) + (-pkin(2) * mrSges(4,2) + (-m(4) * pkin(2) + t380) * pkin(7) + t271) * qJD(2) - t367) * t212 + (-t191 / 0.2e1 + (t328 + t307 / 0.2e1) * qJD(1) + (Ifges(4,3) / 0.2e1 + t364 - Ifges(4,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t282 + (pkin(7) * mrSges(3,2) - qJ(3) * mrSges(4,2) + t266) * qJD(2) - t366) * t209) * qJD(1) + t345 * t131 + t373 * t66 + (t111 * t7 + t373 * t28 + (-t25 + t291) * t18 + (-t26 + t292) * t17) * m(7) + (-t113 * t312 - t230) * t207 + (t113 * t313 - t231) * t210 + (t150 * t40 + t371 * t76 + (t291 - t34) * t20 + (-t292 - t33) * t19) * m(6) + (-t110 * t118 - t157 * t40 + t158 * t39 - t371 * t81 + t372 * t82) * m(5) + t372 * t106 - t297 * t371 + t368 * (-pkin(9) + t158) + m(4) * (qJ(3) * t156 + qJD(3) * t164) + (-m(4) * t133 - t140) * (pkin(2) * t283 - t189) + t156 * mrSges(4,3) + qJD(3) * t166 + t150 * t16 - t118 * t90 + t111 * t15 - t25 * t72 - t34 * t73 - t33 * t74 - t26 * t75 + (t242 * t325 + t245 * t326 + t323 * t342 + t343 * t362 - t338 + t341 + t348) * t129 + (-t158 * mrSges(5,3) - t252) * t89 + (-mrSges(5,3) * t157 - Ifges(5,5)) * t88; -qJD(2) * t166 - t313 * t130 + t312 * t128 + ((t140 - t90) * t209 + t264 * t279) * qJD(1) + (-t88 * mrSges(5,3) - qJD(2) * t106 - t15 - t16 + (-t312 * t207 + t313 * t210 + t106) * qJD(4)) * t211 + (-t89 * mrSges(5,3) + t218 + t273 * (-t66 + t297)) * t208 - m(4) * (qJD(2) * t164 - t133 * t283) - m(7) * (t128 * t17 + t130 * t18) - m(6) * (-t128 * t19 + t130 * t20) + 0.2e1 * ((t17 * t278 + t18 * t277 - t7) * t333 + (-t19 * t278 + t20 * t277 - t40) * t334 + (-t40 + t293) * t335) * t211 + 0.2e1 * ((qJD(4) * t28 + t222) * t333 + (qJD(4) * t76 + t221) * t334 + (t39 - t350) * t335 + (t81 * t335 - m(7) * t28 / 0.2e1 - m(6) * t76 / 0.2e1) * qJD(2)) * t208 + (-t110 * t283 - t280 * t82) * m(5); (t160 * t7 - t17 * t31 - t18 * t30 + t28 * t351) * m(7) + t351 * t66 + t368 * pkin(9) + t230 * t207 + t231 * t210 + (-pkin(4) * t40 - t19 * t41 - t20 * t42 - t76 * t82) * m(6) + t160 * t15 - t81 * t106 + Ifges(5,5) * t88 - t30 * t72 - t42 * t73 - t41 * t74 - t31 * t75 - pkin(4) * t16 + t252 * t89 + (t214 - t348 - t119 / 0.2e1 - t370) * t129 + t297 * t82 + (t120 / 0.2e1 + t337) * t131 + t378; (Ifges(7,3) * t103 - t302) * t326 + t49 * t323 + (t102 * t17 + t103 * t18) * mrSges(7,2) + t57 + t56 - t55 + t54 + t224 - t76 * (mrSges(6,1) * t103 - mrSges(6,2) * t102) - t28 * (mrSges(7,1) * t103 + mrSges(7,3) * t102) + qJD(6) * t72 - t65 * t66 + qJ(6) * t21 - pkin(5) * t23 + t86 + t87 + (t309 + t312) * t20 + (-t310 - t313) * t19 + (-t374 * t102 + t381 * t103) * t321 + (-pkin(5) * t2 + qJ(6) * t1 - t17 * t20 + t18 * t347 - t28 * t65) * m(7) + (-Ifges(6,2) * t103 - t101 + t354) * t325 + (-t382 * t102 + t100 - t306 + t46) * t324; t103 * t66 - t123 * t72 + 0.2e1 * (t2 / 0.2e1 + t28 * t323 + t18 * t321) * m(7) + t23;];
tauc  = t27(:);
