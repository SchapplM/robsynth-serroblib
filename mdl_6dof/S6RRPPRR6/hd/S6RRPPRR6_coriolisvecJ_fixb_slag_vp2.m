% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:13:06
% EndTime: 2019-03-09 09:13:25
% DurationCPUTime: 9.41s
% Computational Cost: add. (9202->534), mult. (22146->702), div. (0->0), fcn. (15501->8), ass. (0->238)
t354 = Ifges(6,1) / 0.2e1;
t208 = -qJD(2) + qJD(5);
t210 = sin(pkin(10));
t211 = cos(pkin(10));
t218 = cos(qJ(2));
t263 = qJD(1) * t218;
t215 = sin(qJ(2));
t264 = qJD(1) * t215;
t136 = -t210 * t264 - t211 * t263;
t137 = t210 * t263 - t211 * t264;
t214 = sin(qJ(5));
t217 = cos(qJ(5));
t322 = t217 * t136 + t214 * t137;
t267 = qJD(6) - t322;
t213 = sin(qJ(6));
t216 = cos(qJ(6));
t235 = t136 * t214 - t217 * t137;
t66 = t208 * t213 + t216 * t235;
t296 = Ifges(7,4) * t66;
t65 = t208 * t216 - t213 * t235;
t25 = Ifges(7,2) * t65 + Ifges(7,6) * t267 + t296;
t64 = Ifges(7,4) * t65;
t26 = t66 * Ifges(7,1) + Ifges(7,5) * t267 + t64;
t298 = t216 / 0.2e1;
t299 = -t213 / 0.2e1;
t228 = t25 * t299 + t26 * t298;
t240 = Ifges(7,5) * t216 - Ifges(7,6) * t213;
t281 = Ifges(7,4) * t216;
t241 = -Ifges(7,2) * t213 + t281;
t282 = Ifges(7,4) * t213;
t242 = Ifges(7,1) * t216 - t282;
t309 = -t267 / 0.2e1;
t311 = -t66 / 0.2e1;
t312 = -t65 / 0.2e1;
t243 = mrSges(7,1) * t213 + mrSges(7,2) * t216;
t294 = pkin(8) * t137;
t200 = pkin(7) * t264;
t164 = qJ(4) * t264 - t200;
t219 = -pkin(2) - pkin(3);
t252 = t219 * qJD(2);
t122 = qJD(3) + t252 - t164;
t201 = pkin(7) * t263;
t167 = -qJ(4) * t263 + t201;
t209 = qJD(2) * qJ(3);
t148 = t167 + t209;
t74 = t211 * t122 - t148 * t210;
t57 = -qJD(2) * pkin(4) + t294 + t74;
t295 = pkin(8) * t136;
t75 = t210 * t122 + t211 * t148;
t59 = t75 + t295;
t33 = -t214 * t59 + t217 * t57;
t31 = -pkin(5) * t208 - t33;
t351 = t31 * t243;
t352 = t235 * t354;
t153 = -qJD(1) * pkin(1) - pkin(2) * t263 - qJ(3) * t264;
t113 = pkin(3) * t263 + qJD(4) - t153;
t83 = -pkin(4) * t136 + t113;
t319 = t83 * mrSges(6,2) + Ifges(6,4) * t322 + t208 * Ifges(6,5) - t240 * t309 - t241 * t312 - t242 * t311 + t228 + t351 + t352;
t353 = -t33 * mrSges(6,3) + t319;
t310 = t66 / 0.2e1;
t349 = t267 / 0.2e1;
t346 = -mrSges(3,1) - mrSges(4,1);
t34 = t214 * t57 + t217 * t59;
t32 = pkin(9) * t208 + t34;
t35 = -pkin(5) * t322 - pkin(9) * t235 + t83;
t10 = -t213 * t32 + t216 * t35;
t11 = t213 * t35 + t216 * t32;
t239 = -t10 * t216 - t11 * t213;
t345 = t239 * mrSges(7,3);
t337 = -t208 * Ifges(6,6) / 0.2e1 - Ifges(6,4) * t235 / 0.2e1;
t340 = -t322 * Ifges(6,2) / 0.2e1;
t331 = t340 + t337;
t52 = pkin(5) * t235 - pkin(9) * t322;
t186 = Ifges(7,2) * t216 + t282;
t187 = Ifges(7,1) * t213 + t281;
t328 = -t216 / 0.2e1;
t300 = Ifges(7,5) * t299 + Ifges(7,6) * t328;
t233 = t210 * t218 - t211 * t215;
t142 = t233 * qJD(2);
t128 = qJD(1) * t142;
t232 = t210 * t215 + t211 * t218;
t143 = t232 * qJD(2);
t129 = qJD(1) * t143;
t49 = qJD(5) * t322 - t128 * t214 + t129 * t217;
t30 = -qJD(6) * t66 - t213 * t49;
t313 = t30 / 0.2e1;
t324 = qJD(6) * t65;
t29 = t216 * t49 + t324;
t314 = t29 / 0.2e1;
t50 = qJD(5) * t235 + t217 * t128 + t129 * t214;
t5 = t29 * Ifges(7,4) + t30 * Ifges(7,2) + t50 * Ifges(7,6);
t6 = Ifges(7,1) * t29 + Ifges(7,4) * t30 + Ifges(7,5) * t50;
t285 = pkin(7) - qJ(4);
t184 = t285 * t218;
t323 = t285 * t215;
t108 = t211 * t184 + t210 * t323;
t259 = qJD(2) * qJD(3);
t61 = -t210 * t259 + (qJD(2) * t108 + t233 * qJD(4)) * qJD(1);
t221 = -t129 * pkin(8) + t61;
t134 = qJD(2) * t184 - qJD(4) * t215;
t118 = t210 * t134;
t132 = -qJD(2) * t323 - qJD(4) * t218;
t62 = t211 * (qJD(1) * t132 + t259) + qJD(1) * t118;
t56 = -pkin(8) * t128 + t62;
t8 = qJD(5) * t33 + t214 * t221 + t217 * t56;
t344 = -(Ifges(6,6) + t300) * t50 + Ifges(6,5) * t49 + t213 * t6 / 0.2e1 + t5 * t298 + t187 * t314 + t186 * t313 + t241 * t324 / 0.2e1 - t8 * mrSges(6,2) + (t240 * t349 + t242 * t310 + t351) * qJD(6);
t248 = t215 * t252;
t196 = qJ(3) * t263;
t204 = t215 * qJD(3);
t266 = qJD(1) * t204 + qJD(2) * t196;
t104 = qJD(1) * t248 + t266;
t70 = pkin(4) * t128 + t104;
t16 = pkin(5) * t50 - pkin(9) * t49 + t70;
t1 = qJD(6) * t10 + t16 * t213 + t216 * t8;
t2 = -qJD(6) * t11 + t16 * t216 - t213 * t8;
t246 = t1 * t216 - t2 * t213;
t343 = t83 * mrSges(6,1) + t10 * mrSges(7,1) - t11 * mrSges(7,2) + t331;
t172 = -qJ(3) * t210 + t211 * t219;
t163 = -pkin(4) + t172;
t173 = t211 * qJ(3) + t210 * t219;
t101 = t214 * t163 + t217 * t173;
t159 = t210 * t217 + t211 * t214;
t98 = -t164 * t210 + t211 * t167;
t229 = t98 + t295;
t99 = t211 * t164 + t210 * t167;
t67 = t99 - t294;
t338 = -t159 * qJD(3) - qJD(5) * t101 + t214 * t67 - t217 * t229;
t179 = t201 + t209;
t181 = mrSges(4,2) * t263 + qJD(2) * mrSges(4,3);
t283 = Ifges(3,4) * t215;
t326 = qJD(2) / 0.2e1;
t327 = -qJD(2) / 0.2e1;
t330 = -Ifges(4,5) * t264 / 0.2e1;
t333 = Ifges(4,3) / 0.2e1;
t336 = -(m(4) * t179 - qJD(2) * mrSges(3,2) + mrSges(3,3) * t263 + t181) * pkin(7) - t179 * mrSges(4,2) - Ifges(4,6) * t327 - t263 * t333 - t330 - Ifges(3,6) * t326 - (t218 * Ifges(3,2) + t283) * qJD(1) / 0.2e1 + t153 * mrSges(4,1);
t279 = qJD(2) * pkin(2);
t175 = qJD(3) + t200 - t279;
t280 = Ifges(4,5) * t218;
t329 = -Ifges(3,4) * t263 / 0.2e1;
t334 = -Ifges(3,1) / 0.2e1;
t335 = (m(4) * t175 + (mrSges(4,2) + mrSges(3,3)) * t264 + t346 * qJD(2)) * pkin(7) - t153 * mrSges(4,3) + (t215 * Ifges(4,1) - t280) * qJD(1) / 0.2e1 - t264 * t334 - t329 + t175 * mrSges(4,2) - (Ifges(4,4) + Ifges(3,5)) * t327;
t332 = Ifges(7,3) * t349 + Ifges(7,5) * t310 + t65 * Ifges(7,6) / 0.2e1;
t325 = mrSges(6,3) * t235;
t284 = -mrSges(6,1) * t208 - mrSges(7,1) * t65 + mrSges(7,2) * t66 + t325;
t321 = t208 * t159;
t320 = -t10 * t213 + t11 * t216;
t318 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t29 + Ifges(7,6) * t30;
t307 = pkin(1) * mrSges(3,1);
t306 = pkin(1) * mrSges(3,2);
t107 = -t184 * t210 + t211 * t323;
t78 = pkin(8) * t233 + t107;
t79 = -pkin(8) * t232 + t108;
t43 = t214 * t79 - t217 * t78;
t9 = qJD(5) * t34 + t214 * t56 - t217 * t221;
t305 = t43 * t9;
t304 = -t137 / 0.2e1;
t302 = -t142 / 0.2e1;
t301 = t143 / 0.2e1;
t156 = t210 * t214 - t217 * t211;
t291 = t156 * t9;
t276 = t136 * Ifges(5,4);
t244 = mrSges(7,1) * t216 - mrSges(7,2) * t213;
t269 = mrSges(6,1) + t244;
t77 = t211 * t132 + t118;
t262 = qJD(2) * t218;
t265 = qJ(3) * t262 + t204;
t261 = qJD(6) * t213;
t260 = qJD(6) * t216;
t176 = -t218 * pkin(2) - t215 * qJ(3) - pkin(1);
t257 = t215 * t279;
t256 = 0.3e1 / 0.2e1 * Ifges(3,4) - 0.3e1 / 0.2e1 * Ifges(4,5);
t255 = Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t254 = Ifges(4,6) / 0.2e1 - Ifges(3,6) / 0.2e1;
t253 = m(4) * pkin(7) + mrSges(4,2);
t251 = t50 * mrSges(6,1) + t49 * mrSges(6,2);
t249 = t128 * mrSges(5,1) + t129 * mrSges(5,2);
t76 = -t132 * t210 + t211 * t134;
t155 = t218 * pkin(3) - t176;
t245 = -t1 * t213 - t2 * t216;
t12 = mrSges(7,1) * t50 - mrSges(7,3) * t29;
t13 = -mrSges(7,2) * t50 + mrSges(7,3) * t30;
t237 = -t213 * t12 + t216 * t13;
t236 = -t210 * t74 + t211 * t75;
t102 = pkin(4) * t232 + t155;
t234 = t214 * t233 - t217 * t232;
t95 = -t214 * t232 - t217 * t233;
t40 = -pkin(5) * t234 - pkin(9) * t95 + t102;
t44 = t214 * t78 + t217 * t79;
t21 = -t213 * t44 + t216 * t40;
t22 = t213 * t40 + t216 * t44;
t100 = t163 * t217 - t173 * t214;
t130 = t219 * t264 + t196;
t41 = -mrSges(7,2) * t267 + mrSges(7,3) * t65;
t42 = mrSges(7,1) * t267 - mrSges(7,3) * t66;
t71 = -mrSges(6,2) * t208 + mrSges(6,3) * t322;
t231 = -t213 * t42 + t216 * t41 + t71;
t230 = -pkin(8) * t143 + t76;
t115 = qJD(2) * mrSges(5,2) + t136 * mrSges(5,3);
t116 = -qJD(2) * mrSges(5,1) + mrSges(5,3) * t137;
t227 = t211 * t115 - t210 * t116 + t181;
t112 = t248 + t265;
t93 = pkin(4) * t137 + t130;
t82 = pkin(4) * t142 + t112;
t224 = qJD(6) * t239 + t246;
t223 = 0.2e1 * t332 + t337 + t343;
t174 = (qJD(3) - t200) * qJD(2);
t165 = (-t218 * mrSges(4,1) - mrSges(4,3) * t215) * qJD(1);
t140 = t156 * qJD(5);
t139 = t156 * qJD(2);
t133 = t257 - t265;
t131 = Ifges(5,4) * t137;
t117 = qJD(1) * t257 - t266;
t106 = -t139 * t216 + t213 * t264;
t105 = t139 * t213 + t216 * t264;
t97 = -pkin(9) + t101;
t96 = pkin(5) - t100;
t92 = -mrSges(5,1) * t136 - mrSges(5,2) * t137;
t85 = -t137 * Ifges(5,1) - Ifges(5,5) * qJD(2) + t276;
t84 = t136 * Ifges(5,2) - Ifges(5,6) * qJD(2) - t131;
t68 = -qJD(3) * t156 + qJD(5) * t100;
t60 = -pkin(8) * t142 + t77;
t55 = qJD(5) * t95 + t217 * t142 + t143 * t214;
t54 = qJD(5) * t234 - t142 * t214 + t143 * t217;
t51 = -mrSges(6,1) * t322 + mrSges(6,2) * t235;
t48 = Ifges(7,3) * t50;
t38 = t214 * t229 + t217 * t67;
t36 = -t52 + t93;
t23 = pkin(5) * t55 - pkin(9) * t54 + t82;
t20 = qJD(5) * t44 + t214 * t60 - t217 * t230;
t19 = -qJD(5) * t43 + t214 * t230 + t217 * t60;
t18 = t213 * t52 + t216 * t33;
t17 = -t213 * t33 + t216 * t52;
t15 = t213 * t36 + t216 * t38;
t14 = -t213 * t38 + t216 * t36;
t7 = -mrSges(7,1) * t30 + mrSges(7,2) * t29;
t4 = -qJD(6) * t22 - t19 * t213 + t216 * t23;
t3 = qJD(6) * t21 + t19 * t216 + t213 * t23;
t24 = [(t6 * t298 + t5 * t299 + Ifges(6,1) * t49 + t70 * mrSges(6,2) + t242 * t314 + t241 * t313 + (mrSges(6,3) + t243) * t9 + t245 * mrSges(7,3) + (-mrSges(7,3) * t320 + t186 * t312 + t187 * t311 + t31 * t244 + t25 * t328 + t26 * t299 + t267 * t300) * qJD(6) + (-Ifges(6,4) + t240 / 0.2e1) * t50) * t95 + (t352 + t345 + t319) * t54 - (-t8 * mrSges(6,3) + t48 / 0.2e1 - Ifges(6,4) * t49 + t70 * mrSges(6,1) + (Ifges(6,2) + Ifges(7,3) / 0.2e1) * t50 + t318) * t234 + (-t117 * mrSges(4,3) + (t254 * qJD(2) + (t176 * mrSges(4,1) - t215 * t256 - 0.2e1 * t307) * qJD(1) + t336) * qJD(2)) * t215 + t104 * (mrSges(5,1) * t232 - mrSges(5,2) * t233) + (-t129 * t233 + t143 * t304) * Ifges(5,1) + (-t33 * t54 - t34 * t55 + t43 * t49 - t44 * t50) * mrSges(6,3) + (Ifges(5,5) * t143 - Ifges(5,6) * t142) * t327 + (-t117 * mrSges(4,1) + t253 * t174 + (t255 * qJD(2) + (-0.2e1 * t306 - t176 * mrSges(4,3) + t256 * t218 + (-0.3e1 / 0.2e1 * Ifges(4,3) - 0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1) + 0.3e1 / 0.2e1 * Ifges(4,1) + t253 * pkin(7)) * t215) * qJD(1) + t335) * qJD(2)) * t218 + m(4) * (t117 * t176 + t133 * t153) + m(5) * (t104 * t155 + t107 * t61 + t108 * t62 + t112 * t113 + t74 * t76 + t75 * t77) + t85 * t301 + t84 * t302 + (t128 * t232 + t136 * t302) * Ifges(5,2) + (-t107 * t129 - t108 * t128 - t142 * t75 - t143 * t74 - t232 * t62 + t233 * t61) * mrSges(5,3) + (t128 * t233 - t129 * t232 + t136 * t301 - t142 * t304) * Ifges(5,4) + t155 * t249 + t102 * t251 + t113 * (mrSges(5,1) * t142 + mrSges(5,2) * t143) + t21 * t12 + t22 * t13 + t3 * t41 + t4 * t42 + t43 * t7 + t284 * t20 + t19 * t71 + t82 * t51 + m(6) * (t102 * t70 + t19 * t34 - t20 * t33 + t44 * t8 + t82 * t83 + t305) + m(7) * (t1 * t22 + t10 * t4 + t11 * t3 + t2 * t21 + t20 * t31 + t305) + t112 * t92 + t77 * t115 + t76 * t116 + t133 * t165 + (t340 + t223) * t55; -t344 + ((t329 + (t306 + t280 / 0.2e1) * qJD(1) + (-pkin(2) * mrSges(4,2) + (-m(4) * pkin(2) + t346) * pkin(7) + t255) * qJD(2) - t335) * t218 + (t330 + (t307 + t283 / 0.2e1) * qJD(1) + (Ifges(3,2) / 0.2e1 - Ifges(4,1) / 0.2e1 + t333 + t334) * t263 + (pkin(7) * mrSges(3,2) - qJ(3) * mrSges(4,2) + t254) * qJD(2) - t336) * t215) * qJD(1) - t246 * mrSges(7,3) + (t345 + t353) * t322 + t237 * t97 + t231 * t68 + t227 * qJD(3) + (-m(4) * t153 - t165) * (pkin(2) * t264 - t196) + (-Ifges(5,1) * t136 - t131 + t84) * t137 / 0.2e1 + (-t100 * t9 + t101 * t8 - t83 * t93 + (-t38 + t68) * t34 + t338 * t33) * m(6) - t338 * t284 + (-t10 * t14 - t11 * t15 + t224 * t97 - t31 * t338 + t320 * t68 + t9 * t96) * m(7) + (-t34 * mrSges(6,3) - Ifges(7,5) * t311 - Ifges(7,6) * t312 - Ifges(7,3) * t309 + t322 * t354 + t331 + t332 + t343) * t235 + ((-t26 / 0.2e1 + t10 * mrSges(7,3) - t97 * t42) * t216 + (t11 * mrSges(7,3) - t97 * t41 + t25 / 0.2e1) * t213) * qJD(6) + (-Ifges(5,5) * t136 - Ifges(5,6) * t137) * t326 + (-t100 * t49 - t101 * t50) * mrSges(6,3) + (qJD(3) * t236 - t113 * t130 + t172 * t61 + t173 * t62 - t74 * t98 - t75 * t99) * m(5) + (-t128 * t173 - t129 * t172 - t136 * t74 + t137 * t75) * mrSges(5,3) + m(4) * (qJ(3) * t174 + qJD(3) * t179) + t269 * t9 - t15 * t41 - t14 * t42 - t136 * (-Ifges(5,2) * t137 - t276) / 0.2e1 - t61 * mrSges(5,1) + t62 * mrSges(5,2) - t38 * t71 - t93 * t51 + t96 * t7 - t99 * t115 - t98 * t116 + Ifges(5,6) * t128 - Ifges(5,5) * t129 - t130 * t92 + t136 * t85 / 0.2e1 - t113 * (mrSges(5,1) * t137 - mrSges(5,2) * t136) + t174 * mrSges(4,3); -t105 * t42 - t106 * t41 + t139 * t71 + (t49 * mrSges(6,3) + t7) * t156 + (-t128 * t210 - t129 * t211) * mrSges(5,3) - t231 * t140 - t227 * qJD(2) + (-t50 * mrSges(6,3) + (-t213 * t41 - t216 * t42) * qJD(6) + t237) * t159 + (t253 * t262 + (t165 - t51 - t92) * t215) * qJD(1) - m(4) * (qJD(2) * t179 - t153 * t264) + t321 * t284 + (-t10 * t105 - t106 * t11 - t140 * t320 + t224 * t159 + t31 * t321 + t291) * m(7) + (t159 * t8 - t83 * t264 + t291 + (t139 - t140) * t34 - t321 * t33) * m(6) + (-qJD(2) * t236 - t113 * t264 + t210 * t62 + t211 * t61) * m(5); -t136 * t115 - t137 * t116 - t322 * t71 - t284 * t235 + (t267 * t41 + t12) * t216 + (-t267 * t42 + t13) * t213 + t249 + t251 + (-t31 * t235 + t267 * t320 - t245) * m(7) + (t235 * t33 - t322 * t34 + t70) * m(6) + (-t136 * t75 - t137 * t74 + t104) * m(5); -m(7) * (t10 * t17 + t11 * t18) + (-m(7) * pkin(5) - t269) * t9 + (m(7) * (-t10 * t260 - t11 * t261 + t246) - t42 * t260 - t41 * t261 + t237) * pkin(9) + ((-t10 * t267 + t1) * t216 + (-t11 * t267 - t2) * t213) * mrSges(7,3) + t228 * qJD(6) - t223 * t235 - pkin(5) * t7 - t18 * t41 - t17 * t42 - t33 * t71 + ((Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t235 - t353) * t322 + (-m(7) * t31 - t284 + t325) * t34 + t344; t48 - t31 * (mrSges(7,1) * t66 + mrSges(7,2) * t65) + (Ifges(7,1) * t65 - t296) * t311 + t25 * t310 + (Ifges(7,5) * t65 - Ifges(7,6) * t66) * t309 - t10 * t41 + t11 * t42 + (t10 * t65 + t11 * t66) * mrSges(7,3) + (-Ifges(7,2) * t66 + t26 + t64) * t312 + t318;];
tauc  = t24(:);
