% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:45:53
% EndTime: 2019-03-09 08:46:15
% DurationCPUTime: 11.02s
% Computational Cost: add. (9050->526), mult. (22925->684), div. (0->0), fcn. (16512->8), ass. (0->235)
t199 = sin(qJ(6));
t202 = cos(qJ(6));
t224 = mrSges(7,1) * t199 + mrSges(7,2) * t202;
t197 = sin(pkin(10));
t204 = cos(qJ(2));
t263 = cos(pkin(10));
t233 = t263 * t204;
t201 = sin(qJ(2));
t253 = qJD(1) * t201;
t148 = -qJD(1) * t233 + t197 * t253;
t164 = t197 * t204 + t201 * t263;
t150 = t164 * qJD(1);
t200 = sin(qJ(5));
t203 = cos(qJ(5));
t92 = t148 * t203 - t200 * t150;
t256 = qJD(6) - t92;
t194 = -qJD(2) + qJD(5);
t94 = t148 * t200 + t203 * t150;
t73 = t194 * t199 + t202 * t94;
t308 = Ifges(7,4) * t73;
t71 = t194 * t202 - t199 * t94;
t25 = Ifges(7,2) * t71 + Ifges(7,6) * t256 + t308;
t68 = Ifges(7,4) * t71;
t26 = t73 * Ifges(7,1) + Ifges(7,5) * t256 + t68;
t310 = t202 / 0.2e1;
t311 = -t199 / 0.2e1;
t289 = -qJ(3) - pkin(7);
t180 = t289 * t201;
t168 = qJD(1) * t180;
t160 = qJD(2) * pkin(2) + t168;
t181 = t289 * t204;
t169 = qJD(1) * t181;
t258 = t197 * t169;
t99 = t160 * t263 + t258;
t214 = qJD(4) - t99;
t302 = t150 * pkin(8);
t61 = -t302 + (-pkin(3) - pkin(4)) * qJD(2) + t214;
t304 = pkin(8) * t148;
t234 = t263 * t169;
t100 = t197 * t160 - t234;
t95 = qJD(2) * qJ(4) + t100;
t65 = t95 + t304;
t35 = -t200 * t65 + t203 * t61;
t33 = -pkin(5) * t194 - t35;
t363 = t33 * t224 + t25 * t311 + t26 * t310;
t221 = Ifges(7,5) * t202 - Ifges(7,6) * t199;
t281 = Ifges(7,4) * t202;
t222 = -Ifges(7,2) * t199 + t281;
t282 = Ifges(7,4) * t199;
t223 = Ifges(7,1) * t202 - t282;
t317 = -t256 / 0.2e1;
t319 = -t73 / 0.2e1;
t320 = -t71 / 0.2e1;
t350 = -Ifges(6,4) * t92 / 0.2e1;
t356 = -Ifges(6,1) / 0.2e1;
t331 = t94 * t356 + t350;
t252 = qJD(1) * t204;
t174 = -qJD(1) * pkin(1) - pkin(2) * t252 + qJD(3);
t77 = t148 * pkin(3) - t150 * qJ(4) + t174;
t59 = -pkin(4) * t148 - t77;
t205 = -t59 * mrSges(6,2) - t194 * Ifges(6,5) + t221 * t317 + t222 * t320 + t223 * t319 + t331 - t363;
t362 = t35 * mrSges(6,3) + t205;
t344 = Ifges(4,1) + Ifges(5,1);
t361 = Ifges(5,4) + Ifges(4,5);
t212 = -t197 * t201 + t233;
t151 = t212 * qJD(2);
t134 = qJD(1) * t151;
t235 = qJD(2) * t289;
t145 = t204 * qJD(3) + t201 * t235;
t210 = qJD(1) * t145;
t146 = -t201 * qJD(3) + t204 * t235;
t211 = qJD(1) * t146;
t70 = t197 * t210 - t263 * t211;
t208 = -t134 * pkin(8) + t70;
t149 = t164 * qJD(2);
t133 = qJD(1) * t149;
t72 = t197 * t211 + t263 * t210;
t66 = qJD(2) * qJD(4) + t72;
t56 = pkin(8) * t133 + t66;
t10 = qJD(5) * t35 + t200 * t208 + t203 * t56;
t182 = Ifges(7,5) * t199 + Ifges(7,6) * t202;
t183 = Ifges(7,2) * t202 + t282;
t184 = Ifges(7,1) * t199 + t281;
t318 = t73 / 0.2e1;
t49 = qJD(5) * t92 + t133 * t200 + t134 * t203;
t31 = -qJD(6) * t73 - t199 * t49;
t323 = t31 / 0.2e1;
t339 = qJD(6) * t71;
t30 = t202 * t49 + t339;
t324 = t30 / 0.2e1;
t359 = t221 / 0.2e1;
t50 = qJD(5) * t94 - t203 * t133 + t134 * t200;
t5 = t30 * Ifges(7,4) + t31 * Ifges(7,2) + t50 * Ifges(7,6);
t6 = Ifges(7,1) * t30 + Ifges(7,4) * t31 + Ifges(7,5) * t50;
t360 = -(Ifges(6,6) - t182 / 0.2e1) * t50 + Ifges(6,5) * t49 - t10 * mrSges(6,2) + t199 * t6 / 0.2e1 + t5 * t310 + t184 * t324 + t183 * t323 + t222 * t339 / 0.2e1 + (t223 * t318 + t256 * t359 + t363) * qJD(6);
t358 = -mrSges(5,1) - mrSges(4,1);
t294 = t92 * Ifges(6,2);
t346 = t294 / 0.2e1;
t27 = -pkin(5) * t92 - pkin(9) * t94 + t59;
t36 = t200 * t61 + t203 * t65;
t34 = pkin(9) * t194 + t36;
t8 = -t199 * t34 + t202 * t27;
t9 = t199 * t27 + t202 * t34;
t207 = -t59 * mrSges(6,1) - t8 * mrSges(7,1) + t9 * mrSges(7,2) + Ifges(6,4) * t94 + 0.2e1 * Ifges(7,5) * t319 + t194 * Ifges(6,6) + 0.2e1 * Ifges(7,6) * t320 + 0.2e1 * Ifges(7,3) * t317 + t346;
t357 = t36 * mrSges(6,3) + t207;
t290 = -Ifges(4,4) + Ifges(5,5);
t141 = Ifges(4,4) * t148;
t277 = t148 * Ifges(5,5);
t355 = t361 * qJD(2) + t344 * t150 - t141 + t277;
t239 = t263 * pkin(2);
t189 = -t239 - pkin(3);
t186 = -pkin(4) + t189;
t307 = pkin(2) * t197;
t187 = qJ(4) + t307;
t122 = t200 * t186 + t203 * t187;
t103 = t168 * t197 - t234;
t216 = t103 + t304;
t104 = t263 * t168 + t258;
t74 = t104 + t302;
t340 = -qJD(5) * t122 - t203 * t216 + (-qJD(4) + t74) * t200;
t15 = mrSges(7,1) * t50 - mrSges(7,3) * t30;
t16 = -mrSges(7,2) * t50 + mrSges(7,3) * t31;
t220 = -t199 * t15 + t202 * t16;
t248 = qJD(6) * t202;
t249 = qJD(6) * t199;
t41 = -mrSges(7,2) * t256 + mrSges(7,3) * t71;
t42 = mrSges(7,1) * t256 - mrSges(7,3) * t73;
t354 = -t42 * t248 - t41 * t249 + t220;
t52 = pkin(5) * t94 - pkin(9) * t92;
t342 = -Ifges(4,6) + Ifges(5,6);
t288 = -mrSges(6,1) * t194 - mrSges(7,1) * t71 + mrSges(7,2) * t73 + mrSges(6,3) * t94;
t225 = mrSges(7,1) * t202 - mrSges(7,2) * t199;
t341 = mrSges(6,1) + t225;
t287 = mrSges(5,2) * t148;
t120 = qJD(2) * mrSges(5,3) - t287;
t278 = t148 * mrSges(4,3);
t255 = -qJD(2) * mrSges(4,2) + t120 - t278;
t276 = t150 * mrSges(4,3);
t286 = mrSges(5,2) * t150;
t254 = -qJD(2) * t358 - t276 - t286;
t336 = t8 * t248 + t9 * t249;
t261 = Ifges(3,6) * qJD(2);
t283 = Ifges(3,4) * t201;
t334 = t261 / 0.2e1 + (t204 * Ifges(3,2) + t283) * qJD(1) / 0.2e1 + pkin(7) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t252);
t247 = qJD(1) * qJD(2);
t238 = t201 * t247;
t232 = pkin(2) * t238;
t58 = t133 * pkin(3) - t134 * qJ(4) - t150 * qJD(4) + t232;
t48 = pkin(4) * t133 + t58;
t12 = pkin(5) * t50 - pkin(9) * t49 - t48;
t1 = qJD(6) * t8 + t10 * t202 + t12 * t199;
t2 = -qJD(6) * t9 - t10 * t199 + t12 * t202;
t333 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t30 + Ifges(7,6) * t31;
t75 = -mrSges(6,2) * t194 + mrSges(6,3) * t92;
t332 = -t199 * t42 + t202 * t41 + t75;
t327 = -0.2e1 * pkin(1);
t316 = -t148 / 0.2e1;
t315 = t148 / 0.2e1;
t313 = t150 / 0.2e1;
t306 = pkin(7) * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t253);
t11 = qJD(5) * t36 + t200 * t56 - t203 * t208;
t109 = -t263 * t180 - t181 * t197;
t82 = -pkin(8) * t164 + t109;
t110 = t197 * t180 - t263 * t181;
t83 = -pkin(8) * t212 + t110;
t43 = t200 * t83 - t203 * t82;
t303 = t11 * t43;
t299 = t49 * mrSges(6,3);
t298 = t50 * mrSges(6,3);
t291 = mrSges(5,2) + mrSges(4,3);
t190 = -t204 * pkin(2) - pkin(1);
t285 = mrSges(7,3) * t199;
t284 = mrSges(7,3) * t202;
t280 = t109 * t70;
t279 = t134 * mrSges(5,2);
t275 = t150 * Ifges(4,4);
t274 = t164 * t70;
t262 = Ifges(3,5) * qJD(2);
t81 = t263 * t145 + t197 * t146;
t251 = qJD(2) * t201;
t250 = qJD(2) * t203;
t246 = pkin(2) * t253;
t237 = t204 * t247;
t80 = t145 * t197 - t263 * t146;
t98 = -pkin(3) * t212 - t164 * qJ(4) + t190;
t230 = t50 * mrSges(6,1) + t49 * mrSges(6,2);
t229 = t1 * t202 - t199 * t2;
t228 = -t1 * t199 - t2 * t202;
t227 = -t199 * t9 - t202 * t8;
t226 = t199 * t8 - t202 * t9;
t102 = t164 * t203 - t200 * t212;
t218 = -t164 * t200 - t203 * t212;
t69 = pkin(4) * t212 - t98;
t37 = -pkin(5) * t218 - pkin(9) * t102 + t69;
t44 = t200 * t82 + t203 * t83;
t18 = -t199 * t44 + t202 * t37;
t19 = t199 * t37 + t202 * t44;
t121 = t186 * t203 - t187 * t200;
t79 = t150 * pkin(3) + t148 * qJ(4) + t246;
t217 = -pkin(8) * t151 + t80;
t64 = pkin(2) * t251 + t149 * pkin(3) - t151 * qJ(4) - t164 * qJD(4);
t62 = -pkin(4) * t150 - t79;
t55 = -pkin(4) * t149 - t64;
t209 = qJD(6) * t227 + t229;
t191 = Ifges(3,4) * t252;
t158 = Ifges(3,1) * t253 + t191 + t262;
t140 = Ifges(5,5) * t150;
t128 = t134 * mrSges(4,2);
t127 = t133 * mrSges(5,1);
t115 = pkin(5) - t121;
t108 = t150 * t199 + t202 * t250;
t107 = t150 * t202 - t199 * t250;
t97 = mrSges(4,1) * t148 + mrSges(4,2) * t150;
t96 = mrSges(5,1) * t148 - mrSges(5,3) * t150;
t88 = -t148 * Ifges(4,2) + Ifges(4,6) * qJD(2) + t275;
t87 = Ifges(5,6) * qJD(2) + t148 * Ifges(5,3) + t140;
t86 = -qJD(2) * pkin(3) + t214;
t63 = pkin(8) * t149 + t81;
t54 = qJD(5) * t102 - t203 * t149 + t151 * t200;
t53 = qJD(5) * t218 + t149 * t200 + t151 * t203;
t51 = -mrSges(6,1) * t92 + mrSges(6,2) * t94;
t47 = Ifges(7,3) * t50;
t39 = t200 * t216 + t203 * t74;
t32 = -t52 + t62;
t23 = t199 * t52 + t202 * t35;
t22 = -t199 * t35 + t202 * t52;
t21 = qJD(5) * t44 + t200 * t63 - t203 * t217;
t20 = -qJD(5) * t43 + t200 * t217 + t203 * t63;
t17 = pkin(5) * t54 - pkin(9) * t53 + t55;
t14 = t199 * t32 + t202 * t39;
t13 = -t199 * t39 + t202 * t32;
t7 = -mrSges(7,1) * t31 + mrSges(7,2) * t30;
t4 = -qJD(6) * t19 + t17 * t202 - t199 * t20;
t3 = qJD(6) * t18 + t17 * t199 + t20 * t202;
t24 = [(Ifges(5,3) * t315 - Ifges(4,2) * t316 - t100 * mrSges(4,3) - t95 * mrSges(5,2) + (Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1) * qJD(2) + t87 / 0.2e1 - t88 / 0.2e1 + t77 * mrSges(5,1) + t174 * mrSges(4,1) + t290 * t313) * t149 + (t227 * mrSges(7,3) - t205 - t331) * t53 + (t355 / 0.2e1 + Ifges(5,5) * t315 + Ifges(4,4) * t316 - t99 * mrSges(4,3) + t86 * mrSges(5,2) + (Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * qJD(2) - t77 * mrSges(5,3) + t174 * mrSges(4,2) + t344 * t313) * t151 + (t5 * t311 + t6 * t310 + Ifges(6,1) * t49 - t48 * mrSges(6,2) + t223 * t324 + t222 * t323 + (mrSges(6,3) + t224) * t11 + t228 * mrSges(7,3) + (t26 * t311 - t202 * t25 / 0.2e1 + t33 * t225 + t183 * t320 + t184 * t319 + t182 * t317 + t226 * mrSges(7,3)) * qJD(6) + (-Ifges(6,4) + t359) * t50) * t102 + (-t35 * t53 - t36 * t54 + t43 * t49 - t44 * t50) * mrSges(6,3) - (-t10 * mrSges(6,3) + t47 / 0.2e1 - Ifges(6,4) * t49 - t48 * mrSges(6,1) + (Ifges(6,2) + Ifges(7,3) / 0.2e1) * t50 + t333) * t218 + (-mrSges(5,3) * t98 + t291 * t109 + t164 * t344 - t212 * t290) * t134 + (-t261 / 0.2e1 + (mrSges(3,1) * t327 - 0.3e1 / 0.2e1 * t283 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t204) * qJD(1) + (qJD(1) * (-mrSges(4,1) * t212 + mrSges(4,2) * t164) + m(4) * (qJD(1) * t190 + t174) + t97) * pkin(2) - t334) * t251 + (t212 * t66 + t274) * mrSges(5,2) + (t212 * t72 + t274) * mrSges(4,3) + (mrSges(4,1) * t190 + t290 * t164 - (Ifges(4,2) + Ifges(5,3)) * t212 - t291 * t110) * t133 + t58 * (-mrSges(5,1) * t212 - mrSges(5,3) * t164) + t98 * t127 + t190 * t128 + t69 * t230 + (t158 / 0.2e1 - t306 + t262 / 0.2e1 + (mrSges(3,2) * t327 + 0.3e1 / 0.2e1 * Ifges(3,4) * t204) * qJD(1)) * t204 * qJD(2) - t254 * t80 + t255 * t81 + m(4) * (t100 * t81 + t110 * t72 - t80 * t99 + t280) + m(5) * (t110 * t66 + t58 * t98 + t64 * t77 + t80 * t86 + t81 * t95 + t280) + t288 * t21 + (-t294 / 0.2e1 - t207) * t54 + m(7) * (t1 * t19 + t18 * t2 + t21 * t33 + t3 * t9 + t4 * t8 + t303) + m(6) * (t10 * t44 + t20 * t36 - t21 * t35 - t48 * t69 + t55 * t59 + t303) + t18 * t15 + t19 * t16 + t3 * t41 + t4 * t42 + t43 * t7 + t55 * t51 + t20 * t75 + t64 * t96; -(t8 * t284 + t9 * t285 + t331 + t362) * t92 + (-Ifges(4,2) * t150 - t141 + t355) * t315 - t340 * t288 + (-t103 * t86 + t187 * t66 + t189 * t70 - t77 * t79 + (qJD(4) - t104) * t95) * m(5) + (m(7) * t209 + t354) * (-pkin(9) + t122) + Ifges(3,5) * t237 + t358 * t70 + (-mrSges(5,2) * t187 - mrSges(4,3) * t307 + t342) * t133 - (-t148 * t344 + t140 - t275 + t87) * t150 / 0.2e1 + (-mrSges(4,3) * t239 + t361) * t134 - (-t148 * t361 + t150 * t342) * qJD(2) / 0.2e1 + (-mrSges(3,1) * t237 + mrSges(3,2) * t238) * pkin(7) + ((t197 * t72 - t263 * t70) * pkin(2) - t100 * t104 + t103 * t99 - t174 * t246) * m(4) + t254 * t103 - t255 * t104 + t336 * mrSges(7,3) + (m(6) * t36 - m(7) * t226 + t332) * (qJD(4) * t203 + qJD(5) * t121) + t334 * t253 + (t10 * t122 - t11 * t121 + t340 * t35 - t36 * t39 - t59 * t62) * m(6) + (t11 * t115 - t13 * t8 - t14 * t9 - t33 * t340) * m(7) + t341 * t11 - Ifges(3,6) * t238 - t97 * t246 - (t346 + t357) * t94 - t360 + (Ifges(5,3) * t150 - t277) * t316 + t88 * t313 - (Ifges(3,5) * t204 - Ifges(3,6) * t201) * t247 / 0.2e1 - t99 * t278 - t1 * t284 + t100 * t276 + t189 * t279 + t2 * t285 + t95 * t286 + t86 * t287 + t252 * t306 - t122 * t298 - t121 * t299 - t14 * t41 - t13 * t42 - t62 * t51 + t66 * mrSges(5,3) - t72 * mrSges(4,2) - t39 * t75 - (-Ifges(3,2) * t253 + t158 + t191) * t252 / 0.2e1 + (pkin(1) * (mrSges(3,1) * t201 + mrSges(3,2) * t204) - t201 * (Ifges(3,1) * t204 - t283) / 0.2e1) * qJD(1) ^ 2 - t79 * t96 + t115 * t7 + qJD(4) * t120 - t77 * (mrSges(5,1) * t150 + mrSges(5,3) * t148) - t174 * (mrSges(4,1) * t150 - mrSges(4,2) * t148); t133 * mrSges(4,1) - t134 * mrSges(5,3) + t92 * t75 + t127 + t128 + t288 * t94 + t254 * t150 + t255 * t148 + (-t256 * t41 - t15) * t202 + (t256 * t42 - t16) * t199 - t230 + (t226 * t256 + t33 * t94 + t228) * m(7) + (-t35 * t94 + t36 * t92 + t48) * m(6) + (t148 * t95 - t150 * t86 + t58) * m(5) + (t100 * t148 + t150 * t99 + t232) * m(4); t279 - qJD(2) * t120 - t107 * t42 - t108 * t41 + (-t51 + t96) * t150 + (-qJD(2) * t75 + t332 * qJD(5) - t299 - t7) * t203 + (-t298 + (-t199 * t41 - t202 * t42) * qJD(6) + t220 + t194 * t288) * t200 + (-t107 * t8 - t108 * t9 + (-qJD(5) * t226 - t11) * t203 + (t194 * t33 + t209) * t200) * m(7) + (t10 * t200 - t11 * t203 - t150 * t59 + t194 * (-t200 * t35 + t203 * t36)) * m(6) + (-qJD(2) * t95 + t77 * t150 + t70) * m(5); ((-t256 * t8 + t1) * t202 + (-t256 * t9 - t2) * t199) * mrSges(7,3) + (t350 + (Ifges(6,2) / 0.2e1 + t356) * t94 + t362) * t92 - m(7) * (t22 * t8 + t23 * t9) + (-m(7) * pkin(5) - t341) * t11 + (m(7) * (t229 - t336) + t354) * pkin(9) + t357 * t94 - pkin(5) * t7 - t23 * t41 - t22 * t42 - t35 * t75 + (-m(7) * t33 - t288) * t36 + t360; t47 - t33 * (mrSges(7,1) * t73 + mrSges(7,2) * t71) + (Ifges(7,1) * t71 - t308) * t319 + t25 * t318 + (Ifges(7,5) * t71 - Ifges(7,6) * t73) * t317 - t8 * t41 + t9 * t42 + (t71 * t8 + t73 * t9) * mrSges(7,3) + (-Ifges(7,2) * t73 + t26 + t68) * t320 + t333;];
tauc  = t24(:);
