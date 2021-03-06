% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:48:04
% EndTime: 2019-03-09 03:48:19
% DurationCPUTime: 8.17s
% Computational Cost: add. (8581->478), mult. (22234->605), div. (0->0), fcn. (16542->8), ass. (0->209)
t183 = cos(qJ(6));
t261 = t183 / 0.2e1;
t177 = sin(pkin(10));
t182 = sin(qJ(3));
t178 = cos(pkin(10));
t260 = cos(qJ(3));
t215 = t260 * t178;
t194 = -t182 * t177 + t215;
t134 = t194 * qJD(1);
t146 = t177 * t260 + t182 * t178;
t135 = t146 * qJD(1);
t181 = sin(qJ(5));
t184 = cos(qJ(5));
t92 = -t134 * t181 + t184 * t135;
t315 = t92 * Ifges(6,1) / 0.2e1;
t291 = Ifges(4,1) + Ifges(5,1);
t290 = Ifges(4,5) + Ifges(5,4);
t180 = sin(qJ(6));
t206 = mrSges(7,1) * t180 + mrSges(7,2) * t183;
t174 = -qJD(3) + qJD(5);
t185 = -pkin(3) - pkin(4);
t245 = pkin(7) + qJ(2);
t158 = t245 * t177;
t147 = qJD(1) * t158;
t159 = t245 * t178;
t148 = qJD(1) * t159;
t100 = -t260 * t147 - t182 * t148;
t72 = pkin(8) * t135 + t100;
t306 = qJD(4) - t72;
t62 = qJD(3) * t185 + t306;
t176 = qJD(3) * qJ(4);
t101 = -t182 * t147 + t260 * t148;
t211 = -pkin(8) * t134 + t101;
t66 = t176 + t211;
t36 = -t181 * t66 + t184 * t62;
t33 = -pkin(5) * t174 - t36;
t314 = t33 * t206;
t246 = -t178 * pkin(2) - pkin(1);
t155 = qJD(1) * t246 + qJD(2);
t77 = -t134 * pkin(3) - t135 * qJ(4) + t155;
t59 = pkin(4) * t134 - t77;
t90 = -t134 * t184 - t181 * t135;
t27 = -pkin(5) * t90 - pkin(9) * t92 + t59;
t37 = t181 * t62 + t184 * t66;
t34 = pkin(9) * t174 + t37;
t10 = -t180 * t34 + t183 * t27;
t11 = t180 * t27 + t183 * t34;
t301 = -t90 * Ifges(6,2) / 0.2e1;
t313 = t59 * mrSges(6,1) + t10 * mrSges(7,1) - t11 * mrSges(7,2) - Ifges(6,4) * t92 - t174 * Ifges(6,6) + t301;
t224 = qJD(6) - t90;
t74 = t174 * t180 + t183 * t92;
t255 = t74 * Ifges(7,4);
t73 = t174 * t183 - t180 * t92;
t25 = t73 * Ifges(7,2) + Ifges(7,6) * t224 + t255;
t69 = Ifges(7,4) * t73;
t26 = Ifges(7,1) * t74 + Ifges(7,5) * t224 + t69;
t262 = -t180 / 0.2e1;
t195 = t25 * t262 + t26 * t261;
t312 = t59 * mrSges(6,2) + Ifges(6,4) * t90 + t174 * Ifges(6,5) + t195 + t314 + t315;
t311 = t73 / 0.2e1;
t273 = t74 / 0.2e1;
t310 = t224 / 0.2e1;
t308 = t36 * mrSges(6,3);
t202 = -t10 * t183 - t11 * t180;
t307 = t202 * mrSges(7,3);
t305 = Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1;
t240 = Ifges(7,4) * t180;
t162 = Ifges(7,2) * t183 + t240;
t239 = Ifges(7,4) * t183;
t163 = Ifges(7,1) * t180 + t239;
t294 = -t183 / 0.2e1;
t263 = Ifges(7,5) * t262 + Ifges(7,6) * t294;
t136 = t194 * qJD(3);
t116 = qJD(1) * t136;
t137 = t146 * qJD(3);
t117 = qJD(1) * t137;
t48 = qJD(5) * t90 + t116 * t184 + t117 * t181;
t31 = -qJD(6) * t74 - t180 * t48;
t277 = t31 / 0.2e1;
t30 = qJD(6) * t73 + t183 * t48;
t278 = t30 / 0.2e1;
t203 = Ifges(7,5) * t183 - Ifges(7,6) * t180;
t204 = -Ifges(7,2) * t180 + t239;
t205 = Ifges(7,1) * t183 - t240;
t298 = t203 * t310 + t204 * t311 + t205 * t273;
t49 = qJD(5) * t92 + t116 * t181 - t184 * t117;
t5 = t30 * Ifges(7,4) + t31 * Ifges(7,2) + t49 * Ifges(7,6);
t6 = t30 * Ifges(7,1) + t31 * Ifges(7,4) + t49 * Ifges(7,5);
t193 = t146 * qJD(2);
t65 = qJD(1) * t193 + qJD(3) * t101;
t190 = -t116 * pkin(8) + t65;
t165 = qJD(2) * t215;
t214 = qJD(3) * t260;
t217 = qJD(1) * qJD(2);
t64 = -t147 * t214 + qJD(1) * t165 + (-qJD(3) * t148 - t177 * t217) * t182;
t63 = qJD(3) * qJD(4) + t64;
t53 = pkin(8) * t117 + t63;
t8 = qJD(5) * t36 + t181 * t190 + t184 * t53;
t304 = -(Ifges(6,6) + t263) * t49 + Ifges(6,5) * t48 + t180 * t6 / 0.2e1 + t5 * t261 + t163 * t278 + t162 * t277 - t8 * mrSges(6,2) + (t314 + t298) * qJD(6);
t197 = qJ(4) * t116 + qJD(4) * t135;
t52 = t117 * t185 + t197;
t12 = pkin(5) * t49 - pkin(9) * t48 + t52;
t1 = qJD(6) * t10 + t12 * t180 + t183 * t8;
t2 = -qJD(6) * t11 + t12 * t183 - t180 * t8;
t209 = t1 * t183 - t2 * t180;
t300 = -mrSges(5,1) - mrSges(4,1);
t247 = Ifges(5,5) - Ifges(4,4);
t289 = Ifges(5,6) - Ifges(4,6);
t124 = Ifges(4,4) * t134;
t238 = Ifges(5,5) * t134;
t299 = t290 * qJD(3) + t291 * t135 + t124 - t238;
t51 = pkin(5) * t92 - pkin(9) * t90;
t296 = Ifges(7,5) * t273 + Ifges(7,6) * t311 + Ifges(7,3) * t310;
t292 = mrSges(6,3) * t92;
t244 = -mrSges(6,1) * t174 - mrSges(7,1) * t73 + mrSges(7,2) * t74 + t292;
t157 = t184 * qJ(4) + t181 * t185;
t288 = qJD(5) * t157 + t181 * t306 + t184 * t211;
t103 = -t182 * t158 + t260 * t159;
t284 = -t10 * t180 + t11 * t183;
t283 = qJD(4) - t100;
t282 = t298 + t312;
t281 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t30 + Ifges(7,6) * t31;
t280 = (m(3) * qJ(2) + mrSges(3,3)) * (t177 ^ 2 + t178 ^ 2);
t276 = -t73 / 0.2e1;
t274 = -t74 / 0.2e1;
t272 = -t224 / 0.2e1;
t102 = t260 * t158 + t159 * t182;
t81 = -pkin(8) * t146 + t102;
t82 = -pkin(8) * t194 + t103;
t43 = t181 * t82 - t184 * t81;
t9 = qJD(5) * t37 + t181 * t53 - t184 * t190;
t269 = t43 * t9;
t268 = t134 / 0.2e1;
t267 = -t134 / 0.2e1;
t265 = t135 / 0.2e1;
t249 = qJD(3) / 0.2e1;
t248 = mrSges(5,2) + mrSges(4,3);
t243 = mrSges(4,3) * t134;
t242 = mrSges(4,3) * t135;
t241 = Ifges(4,4) * t135;
t236 = t102 * t65;
t234 = t146 * t65;
t207 = mrSges(7,1) * t183 - mrSges(7,2) * t180;
t229 = mrSges(6,1) + t207;
t110 = mrSges(5,2) * t134 + qJD(3) * mrSges(5,3);
t223 = -qJD(3) * mrSges(4,2) + t110 + t243;
t222 = -mrSges(5,2) * t135 - qJD(3) * t300 - t242;
t95 = t135 * pkin(3) - t134 * qJ(4);
t220 = qJD(3) * t184;
t219 = qJD(6) * t180;
t218 = qJD(6) * t183;
t70 = t137 * pkin(3) - t136 * qJ(4) - t146 * qJD(4);
t97 = -pkin(3) * t194 - t146 * qJ(4) + t246;
t68 = -pkin(4) * t135 - t95;
t210 = t49 * mrSges(6,1) + t48 * mrSges(6,2);
t208 = -t1 * t180 - t2 * t183;
t15 = mrSges(7,1) * t49 - mrSges(7,3) * t30;
t16 = -mrSges(7,2) * t49 + mrSges(7,3) * t31;
t200 = -t180 * t15 + t183 * t16;
t198 = -t146 * t181 - t184 * t194;
t71 = pkin(4) * t194 - t97;
t99 = t146 * t184 - t181 * t194;
t35 = -pkin(5) * t198 - pkin(9) * t99 + t71;
t44 = t181 * t81 + t184 * t82;
t20 = -t180 * t44 + t183 * t35;
t21 = t180 * t35 + t183 * t44;
t156 = -qJ(4) * t181 + t184 * t185;
t57 = -pkin(4) * t137 - t70;
t41 = -mrSges(7,2) * t224 + mrSges(7,3) * t73;
t42 = mrSges(7,1) * t224 - mrSges(7,3) * t74;
t75 = -mrSges(6,2) * t174 + mrSges(6,3) * t90;
t196 = -t180 * t42 + t183 * t41 + t75;
t78 = -t158 * t214 + t165 + (-qJD(2) * t177 - qJD(3) * t159) * t182;
t192 = qJD(6) * t202 + t209;
t79 = qJD(3) * t103 + t193;
t191 = -t136 * pkin(8) + t79;
t189 = 0.2e1 * t296 + t313;
t154 = -pkin(9) + t157;
t153 = pkin(5) - t156;
t123 = Ifges(5,5) * t135;
t118 = qJD(4) * t184 + qJD(5) * t156;
t112 = t117 * mrSges(5,1);
t111 = t116 * mrSges(4,2);
t105 = t135 * t180 + t183 * t220;
t104 = t135 * t183 - t180 * t220;
t96 = -mrSges(5,1) * t134 - mrSges(5,3) * t135;
t94 = t176 + t101;
t93 = -qJD(3) * pkin(3) + t283;
t86 = t134 * Ifges(4,2) + Ifges(4,6) * qJD(3) + t241;
t85 = Ifges(5,6) * qJD(3) - t134 * Ifges(5,3) + t123;
t61 = pkin(3) * t117 - t197;
t60 = pkin(8) * t137 + t78;
t56 = qJD(5) * t99 + t136 * t181 - t184 * t137;
t55 = qJD(5) * t198 + t136 * t184 + t137 * t181;
t50 = -mrSges(6,1) * t90 + mrSges(6,2) * t92;
t47 = Ifges(7,3) * t49;
t39 = t181 * t211 + t184 * t72;
t32 = -t51 + t68;
t23 = t180 * t51 + t183 * t36;
t22 = -t180 * t36 + t183 * t51;
t19 = pkin(5) * t56 - pkin(9) * t55 + t57;
t18 = qJD(5) * t44 + t181 * t60 - t184 * t191;
t17 = -qJD(5) * t43 + t181 * t191 + t184 * t60;
t14 = t180 * t32 + t183 * t39;
t13 = -t180 * t39 + t183 * t32;
t7 = -mrSges(7,1) * t31 + mrSges(7,2) * t30;
t4 = -qJD(6) * t21 - t17 * t180 + t183 * t19;
t3 = qJD(6) * t20 + t17 * t183 + t180 * t19;
t24 = [(t301 + t189) * t56 + (-t36 * t55 - t37 * t56 + t43 * t48 - t44 * t49) * mrSges(6,3) + m(6) * (t17 * t37 - t18 * t36 + t44 * t8 + t52 * t71 + t57 * t59 + t269) + m(7) * (t1 * t21 + t10 * t4 + t11 * t3 + t18 * t33 + t2 * t20 + t269) + t71 * t210 + (Ifges(5,3) * t267 - Ifges(4,2) * t268 - t101 * mrSges(4,3) - t94 * mrSges(5,2) + t85 / 0.2e1 + t77 * mrSges(5,1) - t86 / 0.2e1 + t155 * mrSges(4,1) + t247 * t265 + t289 * t249) * t137 + (t315 + t307 + t282) * t55 + t244 * t18 + t246 * t111 + m(5) * (t103 * t63 + t61 * t97 + t70 * t77 + t78 * t94 + t79 * t93 + t236) + m(4) * (-t100 * t79 + t101 * t78 + t103 * t64 + t236) - t222 * t79 + t223 * t78 + t97 * t112 + 0.2e1 * t280 * t217 - (-t8 * mrSges(6,3) + t47 / 0.2e1 - Ifges(6,4) * t48 + t52 * mrSges(6,1) + (Ifges(6,2) + Ifges(7,3) / 0.2e1) * t49 + t281) * t198 + (t194 * t63 + t234) * mrSges(5,2) + (t194 * t64 + t234) * mrSges(4,3) + (-mrSges(5,3) * t97 + t248 * t102 + t146 * t291 - t194 * t247) * t116 + (t246 * mrSges(4,1) + t247 * t146 - (Ifges(5,3) + Ifges(4,2)) * t194 - t248 * t103) * t117 + t61 * (-mrSges(5,1) * t194 - mrSges(5,3) * t146) + t20 * t15 + t21 * t16 + t3 * t41 + t4 * t42 + t43 * t7 + t57 * t50 + t17 * t75 + t70 * t96 + (t299 / 0.2e1 + t155 * mrSges(4,2) + t93 * mrSges(5,2) - t100 * mrSges(4,3) - t77 * mrSges(5,3) + Ifges(4,4) * t268 + Ifges(5,5) * t267 + t249 * t290 + t265 * t291) * t136 + (t5 * t262 + t6 * t261 + Ifges(6,1) * t48 + t52 * mrSges(6,2) + t205 * t278 + t204 * t277 + (mrSges(6,3) + t206) * t9 + t208 * mrSges(7,3) + (-mrSges(7,3) * t284 + t162 * t276 + t163 * t274 + t33 * t207 + t224 * t263 + t25 * t294 + t26 * t262) * qJD(6) + (-Ifges(6,4) + t203 / 0.2e1) * t49) * t99; t117 * mrSges(4,1) - t116 * mrSges(5,3) + t90 * t75 + t111 + t112 + t244 * t92 + t222 * t135 - t223 * t134 + (-t224 * t41 - t15) * t183 + (t224 * t42 - t16) * t180 - m(4) * (-t100 * t135 + t101 * t134) - t210 - t280 * qJD(1) ^ 2 + (-t224 * t284 + t33 * t92 + t208) * m(7) + (-t36 * t92 + t37 * t90 - t52) * m(6) + (-t134 * t94 - t135 * t93 + t61) * m(5); ((t10 * mrSges(7,3) - t154 * t42 - t26 / 0.2e1) * t183 + (t11 * mrSges(7,3) - t154 * t41 + t25 / 0.2e1) * t180) * qJD(6) + t200 * t154 + t196 * t118 + (-Ifges(4,2) * t135 + t124 + t299) * t267 + (-t203 * t272 - t204 * t276 - t205 * t274 + t307 - t308 + t312) * t90 + (-t37 * mrSges(6,3) - Ifges(7,5) * t274 - Ifges(7,6) * t276 - Ifges(7,3) * t272 - t305 * t90 + t296 + t313) * t92 - t209 * mrSges(7,3) + t300 * t65 + t289 * t117 + t290 * t116 + t86 * t265 + (t222 + t242) * t101 + t229 * t9 + t288 * t244 + (-t10 * t13 - t11 * t14 + t118 * t284 + t153 * t9 + t154 * t192 + t288 * t33) * m(7) + (-t156 * t9 + t157 * t8 - t59 * t68 + (t118 - t39) * t37 - t288 * t36) * m(6) + (-pkin(3) * t65 + qJ(4) * t63 - t101 * t93 + t283 * t94 - t77 * t95) * m(5) + (-t156 * t48 - t157 * t49) * mrSges(6,3) - t77 * (mrSges(5,1) * t135 - mrSges(5,3) * t134) - t155 * (mrSges(4,1) * t135 + mrSges(4,2) * t134) + (Ifges(5,3) * t135 + t238) * t268 + (-t223 + t243) * t100 - (t291 * t134 + t123 - t241 + t85) * t135 / 0.2e1 - (t134 * t290 + t135 * t289) * qJD(3) / 0.2e1 + (-pkin(3) * t116 - qJ(4) * t117 - t134 * t93 + t135 * t94) * mrSges(5,2) - t304 - t14 * t41 - t13 * t42 + t63 * mrSges(5,3) - t64 * mrSges(4,2) - t68 * t50 - t39 * t75 - t95 * t96 + qJD(4) * t110 + t153 * t7; t116 * mrSges(5,2) - qJD(3) * t110 - t104 * t42 - t105 * t41 + (-t50 + t96) * t135 + (-t48 * mrSges(6,3) - qJD(3) * t75 + t196 * qJD(5) - t7) * t184 + (-t49 * mrSges(6,3) + (-t180 * t41 - t183 * t42) * qJD(6) + t200 + t174 * t244) * t181 + (-t10 * t104 - t105 * t11 + (qJD(5) * t284 - t9) * t184 + (t174 * t33 + t192) * t181) * m(7) + (-t135 * t59 + t181 * t8 - t184 * t9 + t174 * (-t181 * t36 + t184 * t37)) * m(6) + (-qJD(3) * t94 + t135 * t77 + t65) * m(5); -m(7) * (t10 * t22 + t11 * t23) + (-m(7) * pkin(5) - t229) * t9 + (-t42 * t218 - t41 * t219 + m(7) * (-t10 * t218 - t11 * t219 + t209) + t200) * pkin(9) + (t305 * t92 - t282 + t308) * t90 - t189 * t92 + t195 * qJD(6) + ((-t10 * t224 + t1) * t183 + (-t11 * t224 - t2) * t180) * mrSges(7,3) - pkin(5) * t7 - t23 * t41 - t22 * t42 - t36 * t75 + (-m(7) * t33 - t244 + t292) * t37 + t304; t47 - t33 * (mrSges(7,1) * t74 + mrSges(7,2) * t73) + (Ifges(7,1) * t73 - t255) * t274 + t25 * t273 + (Ifges(7,5) * t73 - Ifges(7,6) * t74) * t272 - t10 * t41 + t11 * t42 + (t10 * t73 + t11 * t74) * mrSges(7,3) + (-Ifges(7,2) * t74 + t26 + t69) * t276 + t281;];
tauc  = t24(:);
