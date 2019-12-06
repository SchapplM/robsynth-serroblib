% Calculate vector of inverse dynamics joint torques for
% S5PPRRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:18:38
% EndTime: 2019-12-05 15:18:55
% DurationCPUTime: 5.82s
% Computational Cost: add. (3189->412), mult. (8363->618), div. (0->0), fcn. (7388->14), ass. (0->206)
t274 = m(6) + m(5);
t133 = sin(qJ(4));
t280 = pkin(8) * t133;
t136 = cos(qJ(4));
t174 = pkin(4) * t133 - pkin(9) * t136;
t227 = cos(pkin(5));
t122 = qJD(1) * t227 + qJD(2);
t129 = sin(pkin(5));
t134 = sin(qJ(3));
t127 = sin(pkin(11));
t243 = cos(qJ(3));
t194 = t243 * t127;
t130 = cos(pkin(11));
t131 = cos(pkin(6));
t220 = t130 * t131;
t147 = (t134 * t220 + t194) * t129;
t128 = sin(pkin(6));
t222 = t128 * t134;
t54 = qJD(1) * t147 + t122 * t222;
t279 = -pkin(8) * qJD(5) * t136 + t174 * qJD(4) - t54;
t132 = sin(qJ(5));
t135 = cos(qJ(5));
t210 = qJD(4) * t135;
t213 = qJD(3) * t133;
t103 = -t132 * t213 + t210;
t278 = -t103 / 0.2e1;
t104 = qJD(4) * t132 + t135 * t213;
t277 = -t104 / 0.2e1;
t204 = t136 * qJD(3);
t123 = qJD(5) - t204;
t276 = -t123 / 0.2e1;
t199 = mrSges(5,3) * t213;
t266 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t103 + mrSges(6,2) * t104 + t199;
t170 = -t135 * mrSges(6,1) + t132 * mrSges(6,2);
t153 = m(6) * pkin(4) - t170;
t171 = -mrSges(5,1) * t136 + mrSges(5,2) * t133;
t200 = m(6) * pkin(9) + mrSges(6,3);
t275 = pkin(3) * t274 + t133 * t200 + t136 * t153 + mrSges(4,1) - t171;
t203 = qJD(3) * qJD(4);
t112 = qJDD(3) * t136 - t133 * t203;
t273 = t112 / 0.2e1;
t113 = qJDD(3) * t133 + t136 * t203;
t272 = t113 / 0.2e1;
t117 = -pkin(4) * t136 - pkin(9) * t133 - pkin(3);
t179 = t243 * t220;
t159 = t129 * t179;
t155 = qJD(1) * t159;
t221 = t129 * t134;
t193 = qJD(1) * t221;
t195 = t128 * t243;
t139 = -t122 * t195 + t127 * t193 - t155;
t206 = qJD(5) * t135;
t216 = t135 * t136;
t271 = t117 * t206 + t279 * t132 + t216 * t139 - t210 * t280;
t208 = qJD(5) * t132;
t211 = qJD(4) * t133;
t218 = t132 * t136;
t270 = t132 * t211 * pkin(8) - t117 * t208 + t279 * t135 - t218 * t139;
t62 = qJD(5) * t103 + qJDD(4) * t132 + t113 * t135;
t63 = -qJD(5) * t104 + qJDD(4) * t135 - t113 * t132;
t31 = -mrSges(6,1) * t63 + mrSges(6,2) * t62;
t269 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t113 + t31;
t268 = mrSges(5,1) + t153;
t267 = mrSges(5,2) - t200;
t74 = -mrSges(5,1) * t112 + mrSges(5,2) * t113;
t265 = mrSges(4,1) * qJDD(3) - t74;
t125 = Ifges(5,4) * t204;
t99 = Ifges(6,4) * t103;
t57 = t104 * Ifges(6,1) + t123 * Ifges(6,5) + t99;
t264 = Ifges(5,1) * t213 + Ifges(5,5) * qJD(4) + t135 * t57 + t125;
t105 = t171 * qJD(3);
t263 = qJD(3) * mrSges(4,1) - t105;
t225 = sin(pkin(10));
t172 = t227 * t225;
t226 = cos(pkin(10));
t143 = t127 * t226 + t130 * t172;
t184 = t129 * t225;
t262 = -t128 * t184 + t143 * t131;
t173 = t227 * t226;
t142 = t127 * t225 - t130 * t173;
t185 = t129 * t226;
t261 = t128 * t185 + t131 * t142;
t209 = qJD(4) * t136;
t121 = qJDD(1) * t227 + qJDD(2);
t177 = qJD(3) * t193;
t178 = qJD(3) * t195;
t180 = t129 * t194;
t190 = qJDD(1) * t221;
t29 = qJD(3) * t155 + qJDD(1) * t180 + t121 * t222 + t122 * t178 - t127 * t177 + t190 * t220;
t27 = qJDD(3) * pkin(8) + t29;
t52 = qJD(3) * pkin(8) + t54;
t197 = t129 * t130 * t128;
t81 = -qJDD(1) * t197 + t121 * t131;
t82 = -qJD(1) * t197 + t122 * t131;
t7 = t133 * t81 + t136 * t27 + t82 * t209 - t211 * t52;
t37 = t133 * t82 + t136 * t52;
t8 = -t37 * qJD(4) - t133 * t27 + t136 * t81;
t259 = -t133 * t8 + t136 * t7;
t212 = qJD(3) * t134;
t192 = t128 * t212;
t30 = -qJD(3) * qJD(1) * t180 + qJDD(1) * t159 + t121 * t195 - t122 * t192 - t127 * t190 - t177 * t220;
t28 = -qJDD(3) * pkin(3) - t30;
t11 = -t112 * pkin(4) - t113 * pkin(9) + t28;
t5 = qJDD(4) * pkin(9) + t7;
t35 = qJD(4) * pkin(9) + t37;
t48 = qJD(3) * t117 + t139;
t9 = -t132 * t35 + t135 * t48;
t1 = qJD(5) * t9 + t11 * t132 + t135 * t5;
t10 = t132 * t48 + t135 * t35;
t2 = -qJD(5) * t10 + t11 * t135 - t132 * t5;
t258 = t1 * t135 - t132 * t2;
t234 = Ifges(5,4) * t133;
t166 = t136 * Ifges(5,2) + t234;
t256 = Ifges(5,6) * qJD(4) / 0.2e1 + qJD(3) * t166 / 0.2e1 + Ifges(6,5) * t277 + Ifges(6,6) * t278 + Ifges(6,3) * t276;
t254 = -g(1) * t184 + g(2) * t185 - g(3) * t227;
t169 = t132 * mrSges(6,1) + t135 * mrSges(6,2);
t253 = -t274 * pkin(8) + mrSges(4,2) - t169;
t36 = -t133 * t52 + t136 * t82;
t34 = -qJD(4) * pkin(4) - t36;
t252 = -m(6) * t34 - t266;
t251 = -t2 * mrSges(6,1) + t1 * mrSges(6,2);
t137 = qJD(3) ^ 2;
t232 = Ifges(6,4) * t104;
t56 = t103 * Ifges(6,2) + t123 * Ifges(6,6) + t232;
t250 = -t56 / 0.2e1;
t249 = t62 / 0.2e1;
t248 = t63 / 0.2e1;
t102 = qJDD(5) - t112;
t247 = t102 / 0.2e1;
t245 = t104 / 0.2e1;
t6 = -qJDD(4) * pkin(4) - t8;
t239 = t133 * t6;
t233 = Ifges(5,4) * t136;
t231 = Ifges(6,4) * t132;
t230 = Ifges(6,4) * t135;
t229 = t133 * t139;
t223 = qJD(3) * mrSges(4,2);
t219 = t132 * t133;
t217 = t133 * t135;
t214 = mrSges(4,2) * qJDD(3);
t207 = qJD(5) * t133;
t202 = Ifges(6,5) * t62 + Ifges(6,6) * t63 + Ifges(6,3) * t102;
t198 = mrSges(5,3) * t204;
t191 = t132 * t209;
t186 = t128 * t227;
t183 = t203 / 0.2e1;
t168 = Ifges(6,1) * t135 - t231;
t167 = Ifges(6,1) * t132 + t230;
t165 = -Ifges(6,2) * t132 + t230;
t164 = Ifges(6,2) * t135 + t231;
t163 = Ifges(5,5) * t136 - Ifges(5,6) * t133;
t162 = Ifges(6,5) * t135 - Ifges(6,6) * t132;
t161 = Ifges(6,5) * t132 + Ifges(6,6) * t135;
t148 = t131 * t227 - t197;
t66 = t134 * t186 + t147;
t45 = t133 * t148 + t66 * t136;
t158 = t243 * t186;
t65 = t127 * t221 - t158 - t159;
t19 = t132 * t65 + t135 * t45;
t18 = -t132 * t45 + t135 * t65;
t91 = t131 * t133 + t136 * t222;
t90 = -t136 * t131 + t133 * t222;
t51 = -qJD(3) * pkin(3) + t139;
t156 = t51 * (mrSges(5,1) * t133 + mrSges(5,2) * t136);
t154 = t133 * (Ifges(5,1) * t136 - t234);
t72 = -t132 * t91 - t135 * t195;
t151 = t132 * t195 - t135 * t91;
t150 = -t132 * t207 + t135 * t209;
t149 = t133 * t206 + t191;
t146 = Ifges(6,5) * t133 + t136 * t168;
t145 = Ifges(6,6) * t133 + t136 * t165;
t144 = Ifges(6,3) * t133 + t136 * t162;
t44 = t66 * t133 - t136 * t148;
t119 = -qJD(4) * mrSges(5,2) + t198;
t110 = t174 * qJD(3);
t89 = -t127 * t172 + t130 * t226;
t88 = t127 * t173 + t130 * t225;
t86 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t112;
t84 = pkin(8) * t216 + t117 * t132;
t83 = -pkin(8) * t218 + t117 * t135;
t78 = mrSges(6,1) * t123 - mrSges(6,3) * t104;
t77 = -mrSges(6,2) * t123 + mrSges(6,3) * t103;
t71 = qJD(4) * t91 + t133 * t178;
t70 = -qJD(4) * t90 + t136 * t178;
t68 = t128 * t143 + t131 * t184;
t67 = t128 * t142 - t131 * t185;
t61 = t66 * qJD(3);
t60 = (t158 + (-t127 * t134 + t179) * t129) * qJD(3);
t47 = -mrSges(6,2) * t102 + mrSges(6,3) * t63;
t46 = mrSges(6,1) * t102 - mrSges(6,3) * t62;
t43 = -t262 * t134 + t89 * t243;
t42 = t134 * t89 + t262 * t243;
t41 = -t261 * t134 + t243 * t88;
t40 = t88 * t134 + t261 * t243;
t33 = qJD(5) * t151 - t132 * t70 + t135 * t192;
t32 = qJD(5) * t72 + t132 * t192 + t135 * t70;
t23 = t62 * Ifges(6,1) + t63 * Ifges(6,4) + t102 * Ifges(6,5);
t22 = t62 * Ifges(6,4) + t63 * Ifges(6,2) + t102 * Ifges(6,6);
t21 = t110 * t132 + t135 * t36;
t20 = t110 * t135 - t132 * t36;
t17 = t133 * t68 + t136 * t43;
t15 = t133 * t67 + t136 * t41;
t13 = -qJD(4) * t44 + t60 * t136;
t4 = qJD(5) * t18 + t13 * t135 + t132 * t61;
t3 = -qJD(5) * t19 - t13 * t132 + t135 * t61;
t12 = [t13 * t119 + t45 * t86 + t4 * t77 + t3 * t78 + t18 * t46 + t19 * t47 - t60 * t223 + m(3) * (t121 * t227 + (t127 ^ 2 + t130 ^ 2) * t129 ^ 2 * qJDD(1)) - t66 * t214 + m(4) * (t148 * t81 + t29 * t66 + t54 * t60) + m(5) * (t13 * t37 + t45 * t7) + m(6) * (t1 * t19 + t10 * t4 + t18 * t2 + t3 * t9) + m(2) * qJDD(1) + (-m(4) * t30 + m(5) * t28 - t265) * t65 + (m(4) * t139 + m(5) * t51 - t263) * t61 + (-m(5) * t8 + m(6) * t6 + t269) * t44 + (-m(5) * t36 - t252) * (qJD(4) * t45 + t60 * t133) + (-m(2) - m(3) - m(4) - t274) * g(3); t105 * t192 + t70 * t119 + t32 * t77 + t33 * t78 + t72 * t46 - t151 * t47 + t91 * t86 + t269 * t90 + t266 * t71 + (-mrSges(4,1) * t137 - t214) * t222 + (-mrSges(4,2) * t137 + t265) * t195 + (-t1 * t151 + t10 * t32 + t2 * t72 + t33 * t9 + t34 * t71 + t6 * t90 + t254) * m(6) + (-t36 * t71 + t37 * t70 + t7 * t91 - t8 * t90 + (t212 * t51 - t243 * t28) * t128 + t254) * m(5) + (t81 * t131 + (t243 * t30 + t134 * t29 + (t134 * t139 + t243 * t54) * qJD(3)) * t128 + t254) * m(4) + (t121 + t254) * m(3); t271 * t77 + (Ifges(5,1) * t113 + Ifges(5,4) * t273 + t162 * t247 + t165 * t248 + t168 * t249) * t133 + t270 * t78 + (t1 * t84 + t2 * t83 + (t209 * t34 + t239) * pkin(8) + t229 * t34 + t270 * t9 + t271 * t10) * m(6) + t264 * t209 / 0.2e1 + t263 * t54 + (-g(1) * t43 - g(2) * t41 - g(3) * t66 - t209 * t36 - t211 * t37 + t259) * mrSges(5,3) + (t9 * mrSges(6,1) - t10 * mrSges(6,2) - pkin(8) * t119 - t256) * t211 - (t132 * t57 + t135 * t56) * t207 / 0.2e1 + (t139 * t119 + pkin(8) * t86 - Ifges(6,3) * t247 - Ifges(6,6) * t248 - Ifges(6,5) * t249 - t202 / 0.2e1 + Ifges(5,4) * t272 + Ifges(5,2) * t273 + (-Ifges(5,2) * t133 + t233) * t183 + t251) * t136 + (-pkin(3) * t28 + ((-t37 * t133 - t36 * t136) * qJD(4) + t259) * pkin(8) - t51 * t54 + (-t133 * t36 + t136 * t37) * t139) * m(5) - t139 * t223 + t233 * t272 + t166 * t273 + qJDD(4) * (Ifges(5,5) * t133 + Ifges(5,6) * t136) + (qJD(4) * t146 - t167 * t207) * t245 + t191 * t250 + t83 * t46 + t84 * t47 - pkin(3) * t74 - t29 * mrSges(4,2) + t30 * mrSges(4,1) + t123 * (qJD(4) * t144 - t161 * t207) / 0.2e1 - t22 * t219 / 0.2e1 + t23 * t217 / 0.2e1 + t103 * (qJD(4) * t145 - t164 * t207) / 0.2e1 + t169 * t239 + (t253 * t43 + t275 * t42) * g(1) + (t253 * t66 + t275 * t65) * g(3) + (t253 * t41 + t275 * t40) * g(2) + t266 * (pkin(8) * t209 + t229) + t34 * (mrSges(6,1) * t149 + mrSges(6,2) * t150) + t28 * t171 + qJD(4) ^ 2 * t163 / 0.2e1 + (-t1 * t219 - t10 * t149 - t9 * t150 - t2 * t217) * mrSges(6,3) + qJD(4) * t156 + Ifges(4,3) * qJDD(3) + t154 * t183 + t269 * t280; (t267 * t17 - t268 * (-t133 * t43 + t136 * t68)) * g(1) - (-Ifges(5,2) * t213 + t125 + t264) * t204 / 0.2e1 + (t267 * t45 + t268 * t44) * g(3) + (t267 * t15 - t268 * (-t133 * t41 + t136 * t67)) * g(2) + (t135 * t47 + m(6) * ((-t10 * t132 - t9 * t135) * qJD(5) + t258) - t132 * t46 - t78 * t206 - t77 * t208) * pkin(9) + (-t10 * t208 - t206 * t9 + t258) * mrSges(6,3) + t256 * t213 + (t199 + t252) * t37 + (t103 * t165 + t104 * t168 + t123 * t162) * qJD(5) / 0.2e1 - (t103 * t145 + t104 * t146 + t123 * t144) * qJD(3) / 0.2e1 + (t204 * t56 + t23) * t132 / 0.2e1 + t123 * t34 * t169 + (-t119 + t198) * t36 + t135 * t22 / 0.2e1 + t161 * t247 + t164 * t248 + t167 * t249 + t208 * t250 - t163 * t203 / 0.2e1 + Ifges(5,6) * t112 + Ifges(5,5) * t113 - t21 * t77 - t20 * t78 - pkin(4) * t31 - t7 * mrSges(5,2) + t8 * mrSges(5,1) + t57 * t206 / 0.2e1 + (-t156 - t9 * (mrSges(6,1) * t133 - mrSges(6,3) * t216) - t10 * (-mrSges(6,2) * t133 - mrSges(6,3) * t218)) * qJD(3) + (-pkin(4) * t6 - t10 * t21 - t20 * t9) * m(6) - t137 * t154 / 0.2e1 + t6 * t170 + Ifges(5,3) * qJDD(4); -t34 * (mrSges(6,1) * t104 + mrSges(6,2) * t103) + (Ifges(6,1) * t103 - t232) * t277 + t56 * t245 + (Ifges(6,5) * t103 - Ifges(6,6) * t104) * t276 - t9 * t77 + t10 * t78 - g(1) * ((-t132 * t17 + t135 * t42) * mrSges(6,1) + (-t132 * t42 - t135 * t17) * mrSges(6,2)) - g(2) * ((-t132 * t15 + t135 * t40) * mrSges(6,1) + (-t132 * t40 - t135 * t15) * mrSges(6,2)) - g(3) * (mrSges(6,1) * t18 - mrSges(6,2) * t19) + (t10 * t104 + t103 * t9) * mrSges(6,3) + t202 + (-Ifges(6,2) * t104 + t57 + t99) * t278 - t251;];
tau = t12;
