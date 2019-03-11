% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:10:20
% EndTime: 2019-03-08 22:10:36
% DurationCPUTime: 7.44s
% Computational Cost: add. (5554->500), mult. (13267->669), div. (0->0), fcn. (8977->10), ass. (0->227)
t181 = sin(qJ(5));
t182 = sin(qJ(3));
t185 = cos(qJ(5));
t186 = cos(qJ(3));
t118 = t181 * t182 + t185 * t186;
t243 = qJD(3) * t182;
t284 = pkin(8) - pkin(9);
t130 = t284 * t243;
t187 = cos(qJ(2));
t177 = sin(pkin(6));
t248 = qJD(1) * t177;
t151 = t187 * t248;
t145 = t284 * t182;
t146 = t284 * t186;
t204 = t185 * t145 - t146 * t181;
t220 = qJD(3) * t146;
t323 = qJD(5) * t204 - t118 * t151 - t185 * t130 + t181 * t220;
t183 = sin(qJ(2));
t226 = t183 * t248;
t296 = qJD(3) - qJD(5);
t72 = t296 * t118;
t257 = t181 * t186;
t241 = qJD(3) * t186;
t242 = qJD(3) * t185;
t256 = t182 * t185;
t299 = qJD(5) * t256 + t181 * t241 - t182 * t242;
t73 = -qJD(5) * t257 + t299;
t188 = -pkin(3) - pkin(4);
t227 = t188 * qJD(3);
t219 = t182 * t227;
t169 = t182 * qJD(4);
t249 = qJ(4) * t241 + t169;
t92 = t219 + t249;
t326 = pkin(5) * t73 - pkin(10) * t72 + t226 + t92;
t180 = sin(qJ(6));
t184 = cos(qJ(6));
t136 = -t186 * pkin(3) - t182 * qJ(4) - pkin(2);
t116 = t186 * pkin(4) - t136;
t119 = t256 - t257;
t51 = pkin(5) * t118 - pkin(10) * t119 + t116;
t86 = t145 * t181 + t146 * t185;
t26 = -t180 * t86 + t184 * t51;
t325 = qJD(6) * t26 + t326 * t180 + t323 * t184;
t27 = t180 * t51 + t184 * t86;
t324 = -qJD(6) * t27 - t323 * t180 + t326 * t184;
t322 = qJD(5) * t86 + t119 * t151 - t130 * t181 - t185 * t220;
t246 = qJD(2) * t182;
t132 = qJD(2) * pkin(8) + t226;
t178 = cos(pkin(6));
t258 = t178 * t186;
t152 = qJD(1) * t258;
t87 = -t182 * t132 + t152;
t76 = pkin(9) * t246 + t87;
t321 = qJD(4) - t76;
t211 = Ifges(7,5) * t184 - Ifges(7,6) * t180;
t320 = t211 / 0.2e1;
t319 = -mrSges(5,1) - mrSges(4,1);
t310 = mrSges(5,2) + mrSges(4,3);
t244 = qJD(2) * t186;
t223 = t181 * t244;
t110 = t185 * t246 - t223;
t272 = Ifges(7,4) * t180;
t143 = Ifges(7,2) * t184 + t272;
t271 = Ifges(7,4) * t184;
t144 = Ifges(7,1) * t180 + t271;
t108 = -t181 * t246 - t185 * t244;
t100 = qJD(6) - t108;
t63 = t227 + t321;
t176 = qJD(3) * qJ(4);
t247 = qJD(1) * t182;
t150 = t178 * t247;
t88 = t186 * t132 + t150;
t202 = -pkin(9) * t244 + t88;
t69 = t176 + t202;
t30 = t181 * t63 + t185 * t69;
t25 = -pkin(10) * t296 + t30;
t133 = -qJD(2) * pkin(2) - t151;
t91 = -pkin(3) * t244 - qJ(4) * t246 + t133;
t77 = pkin(4) * t244 - t91;
t40 = -pkin(5) * t108 - pkin(10) * t110 + t77;
t11 = -t180 * t25 + t184 * t40;
t12 = t180 * t40 + t184 * t25;
t312 = -t108 * Ifges(6,2) / 0.2e1;
t83 = -t110 * t180 - t184 * t296;
t84 = t110 * t184 - t180 * t296;
t191 = t77 * mrSges(6,1) + t11 * mrSges(7,1) - t12 * mrSges(7,2) - Ifges(6,4) * t110 + t84 * Ifges(7,5) + Ifges(6,6) * t296 + t83 * Ifges(7,6) + t100 * Ifges(7,3) + t312;
t278 = -t180 / 0.2e1;
t315 = -t184 / 0.2e1;
t279 = Ifges(7,5) * t278 + Ifges(7,6) * t315;
t64 = t72 * qJD(2);
t37 = -qJD(6) * t84 - t180 * t64;
t290 = t37 / 0.2e1;
t36 = qJD(6) * t83 + t184 * t64;
t291 = t36 / 0.2e1;
t214 = mrSges(7,1) * t180 + mrSges(7,2) * t184;
t29 = -t181 * t69 + t185 * t63;
t24 = pkin(5) * t296 - t29;
t306 = t214 * t24;
t212 = -Ifges(7,2) * t180 + t271;
t213 = Ifges(7,1) * t184 - t272;
t286 = t84 / 0.2e1;
t317 = t213 * t286 + t212 * t83 / 0.2e1 + t100 * t320;
t65 = qJD(2) * t299 - qJD(5) * t223;
t218 = qJD(2) * t151;
t221 = pkin(9) * qJD(2) - t132;
t300 = (-t186 * t221 + t150) * qJD(3) + t182 * t218;
t175 = qJD(3) * qJD(4);
t253 = qJD(3) * t152 + t186 * t218;
t47 = t221 * t243 + t175 + t253;
t7 = qJD(5) * t29 + t181 * t300 + t185 * t47;
t318 = -(Ifges(6,6) + t279) * t65 - t7 * mrSges(6,2) + Ifges(6,5) * t64 + t143 * t290 + t144 * t291 + (t306 + t317) * qJD(6) - t191 * t110;
t316 = -Ifges(4,6) / 0.2e1;
t97 = Ifges(6,4) * t108;
t228 = t97 / 0.2e1;
t314 = t246 / 0.2e1;
t313 = t110 * Ifges(6,1) / 0.2e1;
t311 = qJD(3) / 0.2e1;
t39 = t181 * t202 + t185 * t76;
t134 = -qJ(4) * t181 + t185 * t188;
t94 = qJD(4) * t185 + qJD(5) * t134;
t309 = -t39 + t94;
t307 = mrSges(6,3) * t110;
t273 = mrSges(6,1) * t296 - mrSges(7,1) * t83 + mrSges(7,2) * t84 + t307;
t135 = t185 * qJ(4) + t181 * t188;
t308 = qJD(5) * t135 + t181 * t321 + t185 * t202;
t301 = (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t110;
t254 = qJD(4) - t87;
t210 = -t11 * t184 - t12 * t180;
t298 = t210 * mrSges(7,3) + t228;
t297 = -t11 * t180 + t12 * t184;
t162 = qJ(4) * t244;
t250 = -qJD(2) * t169 - qJD(3) * t162;
t62 = (t219 - t226) * qJD(2) - t250;
t18 = pkin(5) * t65 - pkin(10) * t64 + t62;
t1 = qJD(6) * t11 + t18 * t180 + t184 * t7;
t2 = -qJD(6) * t12 + t18 * t184 - t180 * t7;
t295 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t36 + Ifges(7,6) * t37;
t277 = t184 / 0.2e1;
t276 = Ifges(7,4) * t84;
t32 = Ifges(7,2) * t83 + Ifges(7,6) * t100 + t276;
t78 = Ifges(7,4) * t83;
t33 = Ifges(7,1) * t84 + Ifges(7,5) * t100 + t78;
t201 = t277 * t33 + t278 * t32;
t293 = -t77 * mrSges(6,2) + Ifges(6,5) * t296 - t201 - t228 - t306 - t313;
t289 = -t83 / 0.2e1;
t287 = -t84 / 0.2e1;
t260 = t177 * t183;
t101 = t182 * t260 - t258;
t102 = t178 * t182 + t186 * t260;
t206 = t101 * t185 - t102 * t181;
t8 = t30 * qJD(5) + t181 * t47 - t185 * t300;
t283 = t206 * t8;
t282 = t8 * t204;
t281 = -t100 / 0.2e1;
t259 = t177 * t187;
t224 = qJD(2) * t259;
t60 = t132 * t241 + (qJD(3) * t178 + t224) * t247;
t269 = t101 * t60;
t215 = mrSges(7,1) * t184 - mrSges(7,2) * t180;
t262 = mrSges(6,1) + t215;
t121 = (-mrSges(5,1) * t186 - mrSges(5,3) * t182) * qJD(2);
t66 = -mrSges(6,1) * t108 + mrSges(6,2) * t110;
t261 = t121 - t66;
t255 = t186 * t187;
t252 = -qJD(3) * t319 - t246 * t310;
t139 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t244;
t140 = mrSges(5,2) * t244 + qJD(3) * mrSges(5,3);
t251 = t139 + t140;
t245 = qJD(2) * t183;
t240 = qJD(6) * t180;
t239 = qJD(6) * t184;
t237 = qJD(2) * qJD(3);
t234 = pkin(3) * t243;
t233 = Ifges(4,5) / 0.2e1 + Ifges(5,4) / 0.2e1;
t232 = 0.3e1 / 0.2e1 * Ifges(5,5) - 0.3e1 / 0.2e1 * Ifges(4,4);
t230 = t316 + Ifges(5,6) / 0.2e1;
t229 = (-mrSges(4,1) * t186 + mrSges(4,2) * t182) * qJD(2) + t261;
t225 = t177 * t245;
t67 = pkin(5) * t110 - pkin(10) * t108;
t216 = t1 * t184 - t2 * t180;
t19 = mrSges(7,1) * t65 - mrSges(7,3) * t36;
t20 = -mrSges(7,2) * t65 + mrSges(7,3) * t37;
t208 = -t180 * t19 + t184 * t20;
t58 = t101 * t181 + t102 * t185;
t96 = t188 * t246 + t162;
t59 = -t132 * t243 + t253;
t44 = -t180 * t58 + t184 * t259;
t45 = t180 * t259 + t184 * t58;
t197 = qJD(6) * t210 + t216;
t164 = Ifges(5,5) * t246;
t82 = t176 + t88;
t196 = t133 * mrSges(4,1) + t91 * mrSges(5,1) + Ifges(5,6) * t311 - Ifges(5,3) * t244 / 0.2e1 + t164 / 0.2e1 + qJD(3) * t316 - (Ifges(4,4) * t182 + t186 * Ifges(4,2)) * qJD(2) / 0.2e1 - t82 * mrSges(5,2) - t88 * mrSges(4,3);
t165 = Ifges(4,4) * t244;
t79 = -qJD(3) * pkin(3) + t254;
t195 = t133 * mrSges(4,2) + t79 * mrSges(5,2) + (t182 * Ifges(5,1) - Ifges(5,5) * t186) * qJD(2) / 0.2e1 + Ifges(4,1) * t314 + t165 / 0.2e1 - t87 * mrSges(4,3) - t91 * mrSges(5,3) + (Ifges(5,4) + Ifges(4,5)) * t311;
t190 = t211 * t281 + t212 * t289 + t213 * t287 + t293;
t189 = qJD(2) ^ 2;
t129 = -pkin(10) + t135;
t128 = pkin(5) - t134;
t123 = pkin(3) * t246 - t162;
t113 = (mrSges(4,1) * t182 + mrSges(4,2) * t186) * t237;
t112 = (mrSges(5,1) * t182 - mrSges(5,3) * t186) * t237;
t109 = t180 * t246 + t184 * t242;
t107 = -t180 * t242 + t184 * t246;
t99 = t234 - t249;
t89 = mrSges(6,2) * t296 + mrSges(6,3) * t108;
t75 = -qJD(3) * t101 + t186 * t224;
t74 = qJD(3) * t102 + t182 * t224;
t71 = (t226 + t234) * qJD(2) + t250;
t61 = Ifges(7,3) * t65;
t52 = t175 + t59;
t50 = mrSges(7,1) * t100 - mrSges(7,3) * t84;
t49 = -mrSges(7,2) * t100 + mrSges(7,3) * t83;
t46 = -t67 + t96;
t28 = mrSges(6,1) * t65 + mrSges(6,2) * t64;
t22 = qJD(5) * t206 + t181 * t74 + t185 * t75;
t21 = qJD(5) * t58 + t181 * t75 - t74 * t185;
t17 = t180 * t67 + t184 * t29;
t16 = -t180 * t29 + t184 * t67;
t15 = t180 * t46 + t184 * t39;
t14 = -t180 * t39 + t184 * t46;
t13 = -mrSges(7,1) * t37 + mrSges(7,2) * t36;
t10 = qJD(6) * t44 - t180 * t225 + t184 * t22;
t9 = -qJD(6) * t45 - t180 * t22 - t184 * t225;
t6 = t36 * Ifges(7,1) + t37 * Ifges(7,4) + t65 * Ifges(7,5);
t5 = t36 * Ifges(7,4) + t37 * Ifges(7,2) + t65 * Ifges(7,6);
t3 = [t10 * t49 - t206 * t13 + t44 * t19 + t45 * t20 + t22 * t89 + t9 * t50 + t251 * t75 - t252 * t74 + t273 * t21 + (-t206 * t64 - t58 * t65) * mrSges(6,3) + ((-mrSges(3,2) * t189 - t112 - t113 + t28) * t187 + (-mrSges(3,1) * t189 + qJD(2) * t229) * t183) * t177 + m(7) * (t1 * t45 + t10 * t12 + t11 * t9 + t2 * t44 + t21 * t24 - t283) + m(6) * (-t21 * t29 + t22 * t30 - t283 + t58 * t7 + (t187 * t62 - t245 * t77) * t177) + m(5) * (t269 + t102 * t52 + t74 * t79 + t75 * t82 + (-t187 * t71 + t245 * t91) * t177) + m(4) * (t269 + t102 * t59 - t74 * t87 + t75 * t88 + (t133 - t151) * t225) + t310 * t237 * (t101 * t186 - t102 * t182); (t213 * t291 + t212 * t290 + t5 * t278 + t6 * t277 + Ifges(6,1) * t64 + t62 * mrSges(6,2) + (mrSges(6,3) + t214) * t8 + (-t1 * t180 - t2 * t184) * mrSges(7,3) + (-mrSges(7,3) * t297 + t100 * t279 + t143 * t289 + t144 * t287 + t24 * t215 + t33 * t278 + t315 * t32) * qJD(6) + (t320 - Ifges(6,4)) * t65) * t119 + (-m(5) * (t183 * t91 + t255 * t82) + m(6) * t183 * t77 - (pkin(2) * t245 + t133 * t183 + t255 * t88) * m(4) + (mrSges(4,2) * t245 + t187 * t252) * t182 + (-mrSges(4,1) * t245 - t187 * t251) * t186) * t248 + (-t204 * t64 - t29 * t72 - t30 * t73 - t65 * t86) * mrSges(6,3) - t204 * t13 + (t116 * t62 - t29 * t322 + t30 * t323 + t7 * t86 + t77 * t92 - t282) * m(6) + t323 * t89 + t324 * t50 + (t1 * t27 + t324 * t11 + t12 * t325 + t2 * t26 + t322 * t24 - t282) * m(7) + t325 * t49 + (-t71 * mrSges(5,3) + (t230 * qJD(3) + (-m(4) * t88 - m(5) * t82 - t251) * pkin(8) + (t232 * t182 + (-0.3e1 / 0.2e1 * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t186) * qJD(2) + t196) * qJD(3) + (m(4) * t87 - m(5) * t79) * t151 + (t310 + (m(5) + m(4)) * pkin(8)) * t60) * t182 + (-t7 * mrSges(6,3) + t61 / 0.2e1 - Ifges(6,4) * t64 + t62 * mrSges(6,1) + (Ifges(7,3) / 0.2e1 + Ifges(6,2)) * t65 + t295) * t118 + 0.2e1 * (m(5) * (qJD(3) * t79 + t52) / 0.2e1 + m(4) * (-qJD(3) * t87 + t59) / 0.2e1) * t186 * pkin(8) + (-t71 * mrSges(5,1) + t52 * mrSges(5,2) + t59 * mrSges(4,3) + (-pkin(8) * t252 + t233 * qJD(3) - t232 * t244 + t195) * qJD(3)) * t186 + m(5) * (t136 * t71 + t91 * t99) + t322 * t273 + (t312 + t191) * t73 + (t313 - t190 + t298) * t72 - t229 * t226 + t26 * t19 + t27 * t20 + t92 * t66 - pkin(2) * t113 + t116 * t28 + t99 * t121 + t136 * t112; -t318 + t319 * t60 + (-t293 + t298 - t301 + t317) * t108 + (-t108 * t29 - t110 * t30 - t134 * t64 - t135 * t65) * mrSges(6,3) + (-t5 / 0.2e1 - t1 * mrSges(7,3) + t94 * t49 + t129 * t20) * t184 + (-t6 / 0.2e1 + t2 * mrSges(7,3) - t94 * t50 - t129 * t19) * t180 + (-pkin(3) * t60 + qJ(4) * t52 - t123 * t91 + t254 * t82 - t79 * t88) * m(5) + t308 * t273 + (-t11 * t14 - t12 * t15 + t128 * t8 + t197 * t129 + t24 * t308 + t297 * t94) * m(7) + t309 * t89 + (-t134 * t8 + t135 * t7 - t29 * t308 + t30 * t309 - t77 * t96) * m(6) + ((-t165 / 0.2e1 + Ifges(5,5) * t244 / 0.2e1 + (-pkin(3) * mrSges(5,2) + t233) * qJD(3) - t195) * t186 + (-t164 / 0.2e1 + Ifges(4,4) * t314 + (-qJ(4) * mrSges(5,2) + t230) * qJD(3) + (Ifges(5,3) / 0.2e1 - Ifges(4,1) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t244 - t196) * t182) * qJD(2) + ((t11 * mrSges(7,3) - t129 * t50 - t33 / 0.2e1) * t184 + (t12 * mrSges(7,3) - t129 * t49 + t32 / 0.2e1) * t180) * qJD(6) + t252 * t88 + t254 * t140 - t15 * t49 - t14 * t50 + t52 * mrSges(5,3) + t262 * t8 - t59 * mrSges(4,2) - t96 * t66 - t123 * t121 + t128 * t13 - t87 * t139; -qJD(3) * t140 - t107 * t50 - t109 * t49 + (mrSges(5,2) * t241 + t182 * t261) * qJD(2) + (-t64 * mrSges(6,3) - qJD(3) * t89 - t13 + (-t180 * t50 + t184 * t49 + t89) * qJD(5)) * t185 + (-t65 * mrSges(6,3) + (-t180 * t49 - t184 * t50) * qJD(6) + t208 - t296 * t273) * t181 + (-t107 * t11 - t109 * t12 + (qJD(5) * t297 - t8) * t185 + (-t24 * t296 + t197) * t181) * m(7) + (t181 * t7 - t185 * t8 - t246 * t77 - t296 * (-t181 * t29 + t185 * t30)) * m(6) + (-t82 * qJD(3) + t246 * t91 + t60) * m(5); -m(7) * (t11 * t16 + t12 * t17) + (-m(7) * pkin(5) - t262) * t8 + (-t49 * t240 - t50 * t239 + m(7) * (-t11 * t239 - t12 * t240 + t216) + t208) * pkin(10) + t201 * qJD(6) + ((-t100 * t11 + t1) * t184 + (-t100 * t12 - t2) * t180) * mrSges(7,3) + t5 * t277 + t180 * t6 / 0.2e1 + (t29 * mrSges(6,3) - t97 / 0.2e1 + t301 + t190) * t108 - pkin(5) * t13 - t17 * t49 - t16 * t50 - t29 * t89 + (-m(7) * t24 - t273 + t307) * t30 + t318; t61 - t24 * (mrSges(7,1) * t84 + mrSges(7,2) * t83) + (Ifges(7,1) * t83 - t276) * t287 + t32 * t286 + (Ifges(7,5) * t83 - Ifges(7,6) * t84) * t281 - t11 * t49 + t12 * t50 + (t11 * t83 + t12 * t84) * mrSges(7,3) + (-Ifges(7,2) * t84 + t33 + t78) * t289 + t295;];
tauc  = t3(:);
