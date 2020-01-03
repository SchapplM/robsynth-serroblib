% Calculate vector of inverse dynamics joint torques for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP6_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP6_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP6_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:56:46
% EndTime: 2019-12-31 19:57:15
% DurationCPUTime: 15.11s
% Computational Cost: add. (4677->527), mult. (10757->672), div. (0->0), fcn. (7376->10), ass. (0->249)
t374 = Ifges(5,4) + Ifges(6,4);
t169 = sin(qJ(2));
t242 = qJD(1) * qJD(2);
t222 = t169 * t242;
t172 = cos(qJ(2));
t241 = qJDD(1) * t172;
t144 = -t222 + t241;
t145 = qJDD(1) * t169 + t172 * t242;
t266 = sin(pkin(8));
t267 = cos(pkin(8));
t102 = t144 * t266 + t145 * t267;
t140 = t267 * t169 + t266 * t172;
t126 = t140 * qJD(1);
t168 = sin(qJ(4));
t171 = cos(qJ(4));
t108 = qJD(2) * t171 - t126 * t168;
t53 = qJD(4) * t108 + qJDD(2) * t168 + t102 * t171;
t317 = t53 / 0.2e1;
t109 = qJD(2) * t168 + t126 * t171;
t54 = -qJD(4) * t109 + qJDD(2) * t171 - t102 * t168;
t316 = t54 / 0.2e1;
t101 = t144 * t267 - t145 * t266;
t98 = qJDD(4) - t101;
t315 = t98 / 0.2e1;
t357 = Ifges(5,1) + Ifges(6,1);
t355 = Ifges(5,5) + Ifges(6,5);
t354 = Ifges(5,2) + Ifges(6,2);
t353 = Ifges(6,6) + Ifges(5,6);
t352 = Ifges(6,3) + Ifges(5,3);
t293 = t171 * pkin(4);
t159 = pkin(3) + t293;
t199 = -mrSges(6,1) * t171 + mrSges(6,2) * t168;
t201 = -mrSges(5,1) * t171 + mrSges(5,2) * t168;
t373 = -m(5) * pkin(3) - m(6) * t159 + t199 + t201;
t245 = qJD(4) * t168;
t139 = t169 * t266 - t172 * t267;
t125 = t139 * qJD(1);
t264 = t125 * t168;
t372 = t245 + t264;
t328 = m(4) + m(5) + m(6);
t371 = -mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t370 = t355 * t315 + t374 * t316 + t357 * t317;
t369 = t374 * t108;
t151 = -mrSges(3,1) * t172 + mrSges(3,2) * t169;
t165 = qJ(2) + pkin(8);
t163 = sin(t165);
t164 = cos(t165);
t368 = -t164 * mrSges(4,1) - t371 * t163 + t151;
t367 = t374 * t109;
t366 = t374 * t171;
t365 = t374 * t168;
t170 = sin(qJ(1));
t173 = cos(qJ(1));
t364 = g(1) * t173 + g(2) * t170;
t363 = qJD(4) + t125;
t318 = m(6) * pkin(4);
t362 = t353 * t98 + t354 * t54 + t374 * t53;
t361 = t144 / 0.2e1;
t290 = qJD(2) / 0.2e1;
t359 = -mrSges(6,1) - mrSges(5,1);
t358 = mrSges(5,2) + mrSges(6,2);
t351 = t353 * t108 + t355 * t109 + t352 * t363;
t350 = t354 * t108 + t353 * t363 + t367;
t349 = t357 * t109 + t355 * t363 + t369;
t348 = Ifges(4,4) * t125;
t347 = t125 * Ifges(4,2);
t300 = pkin(2) * t172;
t160 = pkin(1) + t300;
t146 = -qJD(1) * t160 + qJD(3);
t346 = t146 * mrSges(4,2);
t268 = qJDD(2) / 0.2e1;
t345 = mrSges(6,1) + t318;
t167 = -qJ(3) - pkin(6);
t150 = t167 * t169;
t152 = t167 * t172;
t105 = t150 * t266 - t152 * t267;
t103 = t171 * t105;
t91 = pkin(3) * t139 - pkin(7) * t140 - t160;
t47 = t168 * t91 + t103;
t230 = t266 * pkin(2);
t157 = t230 + pkin(7);
t251 = qJ(5) + t157;
t207 = qJD(4) * t251;
t263 = t125 * t171;
t143 = qJD(1) * t152;
t129 = t266 * t143;
t142 = qJD(1) * t150;
t100 = t142 * t267 + t129;
t250 = qJD(1) * t169;
t236 = pkin(2) * t250;
t74 = pkin(3) * t126 + pkin(7) * t125 + t236;
t33 = -t100 * t168 + t171 * t74;
t344 = -pkin(4) * t126 - qJ(5) * t263 - qJD(5) * t168 - t171 * t207 - t33;
t212 = t267 * t143;
t99 = t142 * t266 - t212;
t343 = t372 * pkin(4) - t99;
t34 = t171 * t100 + t168 * t74;
t342 = -qJ(5) * t264 + qJD(5) * t171 - t168 * t207 - t34;
t272 = t126 * mrSges(4,3);
t341 = -qJD(2) * mrSges(4,1) - mrSges(5,1) * t108 + mrSges(5,2) * t109 + t272;
t340 = Ifges(4,5) * qJD(2);
t339 = Ifges(4,6) * qJD(2);
t286 = mrSges(6,2) * t171;
t198 = mrSges(6,1) * t168 + t286;
t200 = mrSges(5,1) * t168 + mrSges(5,2) * t171;
t134 = qJD(2) * pkin(2) + t142;
t92 = t134 * t267 + t129;
t82 = -qJD(2) * pkin(3) - t92;
t55 = -t108 * pkin(4) + qJD(5) + t82;
t338 = t55 * t198 + t82 * t200;
t337 = -t353 * t168 + t355 * t171;
t336 = -t354 * t168 + t366;
t335 = t357 * t171 - t365;
t334 = t352 * t98 + t353 * t54 + t355 * t53;
t244 = qJD(4) * t171;
t333 = -t244 - t263;
t161 = pkin(6) * t241;
t137 = -pkin(6) * t222 + t161;
t138 = t145 * pkin(6);
t331 = t137 * t172 + t138 * t169;
t265 = qJDD(1) * pkin(1);
t116 = -pkin(2) * t144 + qJDD(3) - t265;
t31 = -pkin(3) * t101 - pkin(7) * t102 + t116;
t247 = qJD(3) * t169;
t88 = qJDD(2) * pkin(2) - qJ(3) * t145 - qJD(1) * t247 - t138;
t248 = qJD(2) * t169;
t233 = pkin(6) * t248;
t246 = qJD(3) * t172;
t96 = qJ(3) * t144 + t161 + (-t233 + t246) * qJD(1);
t39 = t266 * t88 + t267 * t96;
t37 = qJDD(2) * pkin(7) + t39;
t65 = pkin(3) * t125 - pkin(7) * t126 + t146;
t93 = t266 * t134 - t212;
t83 = qJD(2) * pkin(7) + t93;
t3 = t168 * t31 + t171 * t37 + t65 * t244 - t245 * t83;
t27 = t168 * t65 + t171 * t83;
t4 = -t27 * qJD(4) - t168 * t37 + t171 * t31;
t330 = -t4 * t168 + t3 * t171;
t327 = mrSges(5,1) + t345;
t326 = -m(3) * pkin(6) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t325 = 0.2e1 * t268;
t322 = m(3) * pkin(1) + mrSges(2,1) - t368;
t2 = qJ(5) * t54 + qJD(5) * t108 + t3;
t321 = t4 * mrSges(5,1) - t3 * mrSges(5,2) - t2 * mrSges(6,2);
t26 = -t168 * t83 + t171 * t65;
t18 = -qJ(5) * t109 + t26;
t16 = pkin(4) * t363 + t18;
t19 = qJ(5) * t108 + t27;
t320 = t146 * mrSges(4,1) + t26 * mrSges(5,1) + t16 * mrSges(6,1) - t27 * mrSges(5,2) - t19 * mrSges(6,2);
t314 = -t108 / 0.2e1;
t313 = t108 / 0.2e1;
t312 = -t109 / 0.2e1;
t311 = t109 / 0.2e1;
t310 = -t363 / 0.2e1;
t309 = t363 / 0.2e1;
t308 = t125 / 0.2e1;
t307 = t126 / 0.2e1;
t306 = -t126 / 0.2e1;
t298 = pkin(6) * t172;
t297 = pkin(7) * t163;
t285 = mrSges(5,3) * t108;
t284 = mrSges(5,3) * t109;
t283 = mrSges(6,3) * t108;
t282 = mrSges(6,3) * t109;
t281 = Ifges(3,4) * t169;
t280 = Ifges(3,4) * t172;
t273 = t125 * mrSges(4,3);
t271 = t126 * Ifges(4,4);
t128 = t139 * qJD(2);
t262 = t128 * t168;
t261 = t128 * t171;
t258 = t140 * t168;
t257 = t140 * t171;
t166 = -qJ(5) - pkin(7);
t256 = t163 * t166;
t255 = t168 * t173;
t254 = t170 * t168;
t253 = t170 * t171;
t252 = t171 * t173;
t249 = qJD(1) * t172;
t213 = qJD(2) * t167;
t123 = t169 * t213 + t246;
t124 = t172 * t213 - t247;
t73 = t123 * t267 + t124 * t266;
t127 = t140 * qJD(2);
t235 = pkin(2) * t248;
t75 = pkin(3) * t127 + pkin(7) * t128 + t235;
t238 = t168 * t75 + t171 * t73 + t91 * t244;
t237 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t250) * t298;
t231 = t267 * pkin(2);
t229 = t140 * t244;
t14 = -t54 * mrSges(6,1) + t53 * mrSges(6,2);
t217 = -t245 / 0.2e1;
t215 = -t101 * mrSges(4,1) + t102 * mrSges(4,2);
t214 = -t168 * t73 + t171 * t75;
t46 = -t105 * t168 + t171 * t91;
t158 = -t231 - pkin(3);
t204 = pkin(3) * t164 + t297;
t38 = -t266 * t96 + t267 * t88;
t203 = mrSges(3,1) * t169 + mrSges(3,2) * t172;
t195 = t172 * Ifges(3,2) + t281;
t192 = Ifges(3,5) * t172 - Ifges(3,6) * t169;
t72 = t123 * t266 - t267 * t124;
t104 = -t267 * t150 - t152 * t266;
t187 = t159 * t164 - t256;
t186 = qJ(5) * t128 - qJD(5) * t140;
t185 = pkin(1) * t203;
t121 = -t164 * t255 + t253;
t119 = t164 * t254 + t252;
t182 = t229 - t262;
t181 = t140 * t245 + t261;
t180 = t169 * (Ifges(3,1) * t172 - t281);
t36 = -qJDD(2) * pkin(3) - t38;
t162 = Ifges(3,4) * t249;
t149 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t249;
t147 = t158 - t293;
t136 = t251 * t171;
t135 = t251 * t168;
t133 = Ifges(3,1) * t250 + Ifges(3,5) * qJD(2) + t162;
t132 = Ifges(3,6) * qJD(2) + qJD(1) * t195;
t122 = t164 * t252 + t254;
t120 = -t164 * t253 + t255;
t112 = -qJD(2) * mrSges(4,2) - t273;
t89 = mrSges(4,1) * t125 + mrSges(4,2) * t126;
t85 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t102;
t84 = -qJDD(2) * mrSges(4,2) + mrSges(4,3) * t101;
t79 = t126 * Ifges(4,1) + t340 - t348;
t78 = t271 + t339 - t347;
t71 = pkin(4) * t258 + t104;
t64 = mrSges(5,1) * t363 - t284;
t63 = mrSges(6,1) * t363 - t282;
t62 = -mrSges(5,2) * t363 + t285;
t61 = -mrSges(6,2) * t363 + t283;
t56 = -mrSges(6,1) * t108 + mrSges(6,2) * t109;
t35 = pkin(4) * t182 + t72;
t30 = -qJ(5) * t258 + t47;
t24 = pkin(4) * t139 - qJ(5) * t257 + t46;
t23 = -mrSges(5,2) * t98 + mrSges(5,3) * t54;
t22 = -mrSges(6,2) * t98 + mrSges(6,3) * t54;
t21 = mrSges(5,1) * t98 - mrSges(5,3) * t53;
t20 = mrSges(6,1) * t98 - mrSges(6,3) * t53;
t15 = -mrSges(5,1) * t54 + mrSges(5,2) * t53;
t13 = -t54 * pkin(4) + qJDD(5) + t36;
t12 = -qJD(4) * t47 + t214;
t11 = -t105 * t245 + t238;
t6 = -qJ(5) * t229 + (-qJD(4) * t105 + t186) * t168 + t238;
t5 = pkin(4) * t127 + t186 * t171 + (-t103 + (qJ(5) * t140 - t91) * t168) * qJD(4) + t214;
t1 = pkin(4) * t98 - qJ(5) * t53 - qJD(5) * t109 + t4;
t7 = [(t334 / 0.2e1 + t116 * mrSges(4,1) + t1 * mrSges(6,1) - t39 * mrSges(4,3) - Ifges(4,4) * t102 - Ifges(4,2) * t101 - t325 * Ifges(4,6) + t352 * t315 + t353 * t316 + t355 * t317 + t321) * t139 + m(5) * (t11 * t27 + t12 * t26 + t3 * t47 + t4 * t46) + m(4) * (t105 * t39 - t116 * t160 + t146 * t235 + t73 * t93) + (-t229 / 0.2e1 + t262 / 0.2e1) * t350 + (t181 * t26 - t182 * t27 - t257 * t4 - t258 * t3) * mrSges(5,3) + t172 * t133 * t290 + t257 * t370 + t55 * (mrSges(6,1) * t182 - mrSges(6,2) * t181) + t82 * (mrSges(5,1) * t182 - mrSges(5,2) * t181) + (-t93 * mrSges(4,3) - Ifges(4,4) * t307 - Ifges(4,6) * t290 - t78 / 0.2e1 + t347 / 0.2e1 + t353 * t313 + t355 * t311 + t352 * t309 + t320 + t351 / 0.2e1) * t127 + (-t181 * t374 - t182 * t354) * t313 + (-t181 * t357 - t182 * t374) * t311 - t160 * t215 + Ifges(3,6) * t172 * t268 + t145 * t280 / 0.2e1 + (t180 + t172 * (-Ifges(3,2) * t169 + t280)) * t242 / 0.2e1 + (-m(4) * t38 + m(5) * t36 + t15 - t85) * t104 + (t92 * mrSges(4,3) - Ifges(4,1) * t307 - Ifges(4,5) * t290 - t79 / 0.2e1 + t348 / 0.2e1 - t346) * t128 + t24 * t20 + (t144 * t298 + t331) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t331) + t195 * t361 + t89 * t235 + (-t1 * t257 + t16 * t181 - t182 * t19 - t2 * t258) * mrSges(6,3) - t362 * t258 / 0.2e1 + (-t254 * t318 - t328 * (t173 * t160 - t170 * t167) + t359 * t122 - t358 * t121 + t326 * t170 + (-m(5) * t204 - m(6) * t187 - t322) * t173) * g(2) + (t359 * t120 - t358 * t119 + (t167 * t328 - t168 * t318 + t326) * t173 + (m(4) * t160 - m(6) * (-t160 - t187) - m(5) * (-t160 - t204) + t322) * t170) * g(1) + (Ifges(3,1) * t145 - pkin(6) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t145) + Ifges(3,4) * t361 + t325 * Ifges(3,5)) * t169 + (t192 * t290 - t237) * qJD(2) + (-t181 * t355 - t182 * t353) * t309 - t349 * t261 / 0.2e1 + (-m(4) * t92 + m(5) * t82 + t341) * t72 + m(6) * (t1 * t24 + t13 * t71 + t16 * t5 + t19 * t6 + t2 * t30 + t35 * t55) + (t116 * mrSges(4,2) - t38 * mrSges(4,3) + Ifges(4,1) * t102 + Ifges(4,4) * t101 + Ifges(4,5) * t325 + t13 * t198 + t36 * t200 + t217 * t349 + t315 * t337 + t316 * t336 + t317 * t335) * t140 - t132 * t248 / 0.2e1 - t185 * t242 - t151 * t265 + t30 * t22 + t46 * t21 + t47 * t23 + t35 * t56 + t6 * t61 + t11 * t62 + t5 * t63 + t12 * t64 + t71 * t14 + t105 * t84 + t73 * t112 - pkin(1) * (-mrSges(3,1) * t144 + mrSges(3,2) * t145) - qJDD(2) * mrSges(3,2) * t298 + t172 * (Ifges(3,4) * t145 + Ifges(3,2) * t144 + Ifges(3,6) * qJDD(2)) / 0.2e1 + Ifges(2,3) * qJDD(1) - t149 * t233; (t237 + (-t180 / 0.2e1 + t185) * qJD(1)) * qJD(1) + (-t271 + t351) * t306 + ((t266 * t39 + t267 * t38) * pkin(2) - t100 * t93 - t146 * t236 + t92 * t99) * m(4) + t168 * t370 + (t158 * t36 - t26 * t33 - t27 * t34 - t82 * t99) * m(5) + (-t1 * t168 + t16 * t333 + t171 * t2 - t19 * t372) * mrSges(6,3) + (t26 * t333 - t27 * t372 + t330) * mrSges(5,3) + (t171 * t354 + t365) * t316 + (t168 * t357 + t366) * t317 - (-Ifges(3,2) * t250 + t133 + t162) * t249 / 0.2e1 + t364 * (t203 + t328 * pkin(2) * t169 + (-m(5) * pkin(7) + m(6) * t166 - t371) * t164 + (mrSges(4,1) - t373) * t163) + (-m(4) * t300 - m(5) * (t297 + t300) - m(6) * (-t256 + t300) + t373 * t164 + t368) * g(3) + (t217 - t264 / 0.2e1) * t350 + (pkin(6) * t149 + t132 / 0.2e1) * t250 + t13 * t199 + t36 * t201 + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + (m(5) * ((-t27 * t168 - t26 * t171) * qJD(4) + t330) - t168 * t21 + t171 * t23 - t64 * t244 - t62 * t245) * t157 + t338 * qJD(4) + t93 * t272 + (-t348 + t79) * t308 + (t108 * t336 + t109 * t335 + t337 * t363) * qJD(4) / 0.2e1 + t362 * t171 / 0.2e1 + (t168 * t355 + t171 * t353) * t315 - (Ifges(4,2) * t308 - t339 / 0.2e1 - t353 * t314 - t355 * t312 - t352 * t310 + t320) * t126 - (Ifges(4,1) * t306 - t340 / 0.2e1 - t346 + t336 * t314 + t335 * t312 + t337 * t310 - t338) * t125 - t341 * t99 + t342 * t61 + t343 * t56 + t344 * t63 + (-t1 * t135 + t13 * t147 + t136 * t2 + t16 * t344 + t19 * t342 + t343 * t55) * m(6) + (t244 / 0.2e1 + t263 / 0.2e1) * t349 + t78 * t307 - t192 * t242 / 0.2e1 + t38 * mrSges(4,1) - t39 * mrSges(4,2) - t34 * t62 - t33 * t64 + Ifges(4,6) * t101 + Ifges(4,5) * t102 - t100 * t112 + t84 * t230 + t85 * t231 - t92 * t273 - t135 * t20 + t136 * t22 - t137 * mrSges(3,2) - t138 * mrSges(3,1) + Ifges(3,6) * t144 + Ifges(3,5) * t145 + t147 * t14 + t158 * t15 - t89 * t236; t125 * t112 - (t56 + t341) * t126 + (t20 + t21 + t363 * (t61 + t62)) * t171 + (t22 + t23 - t363 * (t63 + t64)) * t168 + t215 + (-g(1) * t170 + g(2) * t173) * t328 + (t1 * t171 - t126 * t55 + t168 * t2 + t363 * (-t16 * t168 + t171 * t19)) * m(6) + (-t126 * t82 + t168 * t3 + t171 * t4 + t363 * (-t168 * t26 + t171 * t27)) * m(5) + (t125 * t93 + t126 * t92 + t116) * m(4); t350 * t311 + t321 + (t285 - t62) * t26 + (-t109 * t354 + t349 + t369) * t314 + (t284 + t64) * t27 + (-t121 * t327 + t122 * t358) * g(1) + t345 * t1 + (t119 * t327 - t120 * t358) * g(2) + (-m(6) * (-t16 + t18) + t282 + t63) * t19 + t16 * t283 + (t108 * t355 - t109 * t353) * t310 + (t108 * t357 - t367) * t312 + (t168 * t345 + t200 + t286) * g(3) * t163 - t18 * t61 - t55 * (mrSges(6,1) * t109 + mrSges(6,2) * t108) - t82 * (mrSges(5,1) * t109 + mrSges(5,2) * t108) + t334 + ((-m(6) * t55 - t56) * t109 + t20) * pkin(4); -t108 * t61 + t109 * t63 + (g(3) * t164 - t19 * t108 + t16 * t109 - t364 * t163 + t13) * m(6) + t14;];
tau = t7;
