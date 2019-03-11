% Calculate vector of inverse dynamics joint torques for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPPRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:31:38
% EndTime: 2019-03-09 01:31:49
% DurationCPUTime: 6.24s
% Computational Cost: add. (4302->425), mult. (8267->565), div. (0->0), fcn. (5407->14), ass. (0->188)
t138 = cos(pkin(9));
t115 = -pkin(1) * t138 - pkin(2);
t107 = -qJ(4) + t115;
t267 = qJD(1) * qJD(4) - qJDD(1) * t107;
t140 = sin(qJ(6));
t143 = cos(qJ(6));
t135 = sin(pkin(10));
t141 = sin(qJ(5));
t137 = cos(pkin(10));
t235 = cos(qJ(5));
t183 = t235 * t137;
t152 = -t141 * t135 + t183;
t153 = -t135 * t235 - t141 * t137;
t85 = t153 * qJD(1);
t53 = qJD(5) * t85 + qJDD(1) * t152;
t198 = qJD(1) * t135;
t86 = qJD(1) * t183 - t141 * t198;
t63 = qJD(5) * t143 - t140 * t86;
t26 = qJD(6) * t63 + qJDD(5) * t140 + t143 * t53;
t247 = t26 / 0.2e1;
t64 = qJD(5) * t140 + t143 * t86;
t27 = -qJD(6) * t64 + qJDD(5) * t143 - t140 * t53;
t246 = t27 / 0.2e1;
t54 = -t152 * qJD(1) * qJD(5) + qJDD(1) * t153;
t51 = qJDD(6) - t54;
t245 = t51 / 0.2e1;
t136 = sin(pkin(9));
t110 = pkin(1) * t136 + qJ(3);
t266 = t135 ^ 2 + t137 ^ 2;
t131 = qJ(1) + pkin(9);
t123 = sin(t131);
t125 = cos(t131);
t258 = -g(1) * t123 + g(2) * t125;
t265 = Ifges(7,1) * t247 + Ifges(7,4) * t246 + Ifges(7,5) * t245;
t264 = -m(7) - m(6);
t9 = -mrSges(7,1) * t27 + mrSges(7,2) * t26;
t236 = -qJDD(5) * mrSges(6,1) + mrSges(6,3) * t53 + t9;
t90 = qJD(1) * t107 + qJD(3);
t65 = -qJD(2) * t135 + t137 * t90;
t58 = -pkin(7) * qJD(1) * t137 + t65;
t66 = t137 * qJD(2) + t135 * t90;
t59 = -pkin(7) * t198 + t66;
t31 = t141 * t58 + t235 * t59;
t29 = qJD(5) * pkin(8) + t31;
t98 = t110 * qJD(1);
t96 = qJD(4) + t98;
t82 = pkin(4) * t198 + t96;
t35 = -pkin(5) * t85 - pkin(8) * t86 + t82;
t12 = -t140 * t29 + t143 * t35;
t263 = t12 * mrSges(7,1);
t13 = t140 * t35 + t143 * t29;
t262 = t13 * mrSges(7,2);
t221 = t86 * mrSges(6,3);
t218 = -qJD(5) * mrSges(6,1) - mrSges(7,1) * t63 + mrSges(7,2) * t64 + t221;
t261 = qJD(6) - t85;
t77 = qJDD(3) - t267;
t60 = -qJDD(2) * t135 + t137 * t77;
t61 = t137 * qJDD(2) + t135 * t77;
t164 = t61 * t135 + t60 * t137;
t14 = mrSges(7,1) * t51 - mrSges(7,3) * t26;
t15 = -mrSges(7,2) * t51 + mrSges(7,3) * t27;
t259 = -t140 * t14 + t143 * t15;
t191 = qJDD(1) * t135;
t206 = pkin(1) * qJDD(1);
t117 = t136 * t206;
t133 = qJD(1) * qJD(3);
t91 = -qJDD(1) * qJ(3) - t117 - t133;
t89 = qJDD(4) - t91;
t75 = pkin(4) * t191 + t89;
t18 = -pkin(5) * t54 - pkin(8) * t53 + t75;
t179 = qJD(5) * t235;
t197 = qJD(5) * t141;
t190 = qJDD(1) * t137;
t56 = -pkin(7) * t190 + t60;
t57 = -pkin(7) * t191 + t61;
t10 = t141 * t56 + t58 * t179 - t197 * t59 + t235 * t57;
t7 = qJDD(5) * pkin(8) + t10;
t1 = qJD(6) * t12 + t140 * t18 + t143 * t7;
t2 = -qJD(6) * t13 - t140 * t7 + t143 * t18;
t175 = -t1 * t143 + t140 * t2;
t256 = mrSges(3,1) + mrSges(5,3) + mrSges(6,3) - mrSges(4,2);
t11 = -qJD(5) * t31 - t141 * t57 + t235 * t56;
t36 = -mrSges(7,2) * t261 + mrSges(7,3) * t63;
t37 = mrSges(7,1) * t261 - mrSges(7,3) * t64;
t162 = t140 * t36 + t143 * t37;
t40 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t54;
t255 = qJD(6) * t162 - t259 - t40;
t166 = t12 * t143 + t13 * t140;
t254 = qJD(6) * t166 + t175;
t30 = -t141 * t59 + t235 * t58;
t28 = -qJD(5) * pkin(5) - t30;
t253 = -m(7) * t28 - t218;
t87 = -t135 * t179 - t137 * t197;
t88 = -t135 * t197 + t137 * t179;
t252 = t10 * t153 - t11 * t152 - t30 * t87 - t31 * t88;
t251 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t130 = pkin(10) + qJ(5);
t122 = sin(t130);
t124 = cos(t130);
t172 = mrSges(6,1) * t122 + mrSges(6,2) * t124;
t174 = mrSges(5,1) * t135 + mrSges(5,2) * t137;
t250 = mrSges(3,2) - t174 - mrSges(4,3) - t172 - m(7) * (pkin(5) * t122 - pkin(8) * t124) + t124 * mrSges(7,3);
t248 = m(4) * t98 + m(5) * t96 + m(6) * t82 - mrSges(6,1) * t85 + mrSges(6,2) * t86 + t174 * qJD(1);
t244 = -t63 / 0.2e1;
t243 = -t64 / 0.2e1;
t242 = t64 / 0.2e1;
t241 = -t261 / 0.2e1;
t238 = t86 / 0.2e1;
t8 = -qJDD(5) * pkin(5) - t11;
t237 = t8 * t152;
t234 = mrSges(6,3) * t85;
t233 = Ifges(6,4) * t86;
t232 = Ifges(7,4) * t64;
t142 = sin(qJ(1));
t230 = pkin(1) * t142;
t126 = t135 * pkin(4);
t144 = cos(qJ(1));
t127 = t144 * pkin(1);
t219 = -pkin(7) + t107;
t217 = Ifges(7,4) * t140;
t216 = Ifges(7,4) * t143;
t214 = t140 * t85;
t213 = t140 * t152;
t211 = t143 * t85;
t210 = t143 * t87;
t209 = t143 * t152;
t205 = t123 * t140;
t204 = t123 * t143;
t203 = t125 * t140;
t202 = t125 * t143;
t199 = mrSges(5,1) * t191 + mrSges(5,2) * t190;
t196 = qJD(6) * t140;
t195 = qJD(6) * t143;
t194 = -m(5) + t264;
t189 = Ifges(7,5) * t26 + Ifges(7,6) * t27 + Ifges(7,3) * t51;
t188 = m(4) - t194;
t187 = m(7) * pkin(8) + mrSges(7,3);
t22 = Ifges(7,2) * t63 + Ifges(7,6) * t261 + t232;
t186 = -t140 * t22 / 0.2e1;
t184 = t125 * pkin(2) + t123 * qJ(3) + t127;
t180 = -t54 * mrSges(6,1) + t53 * mrSges(6,2);
t178 = t195 / 0.2e1;
t177 = t125 * qJ(3) - t230;
t176 = t266 * mrSges(5,3);
t97 = t110 + t126;
t171 = -mrSges(7,1) * t143 + mrSges(7,2) * t140;
t170 = mrSges(7,1) * t140 + mrSges(7,2) * t143;
t169 = Ifges(7,1) * t143 - t217;
t168 = -Ifges(7,2) * t140 + t216;
t167 = Ifges(7,5) * t143 - Ifges(7,6) * t140;
t165 = -t12 * t140 + t13 * t143;
t163 = -t135 * t66 - t137 * t65;
t38 = -pkin(5) * t153 - pkin(8) * t152 + t97;
t83 = t219 * t135;
t84 = t219 * t137;
t42 = t141 * t84 + t235 * t83;
t19 = -t140 * t42 + t143 * t38;
t20 = t140 * t38 + t143 * t42;
t161 = -pkin(2) * t123 + t177;
t73 = -qJD(5) * mrSges(6,2) + t234;
t159 = -t140 * t37 + t143 * t36 + t73;
t157 = -t141 * t83 + t235 * t84;
t156 = -t140 * t87 - t152 * t195;
t155 = -t152 * t196 + t210;
t154 = t28 * t170;
t150 = m(7) * pkin(5) - t171;
t145 = qJD(1) ^ 2;
t139 = -pkin(7) - qJ(4);
t95 = qJDD(1) * t115 + qJDD(3);
t80 = Ifges(6,4) * t85;
t72 = t122 * t202 - t205;
t71 = t122 * t203 + t204;
t70 = t122 * t204 + t203;
t69 = -t122 * t205 + t202;
t62 = Ifges(7,4) * t63;
t52 = pkin(5) * t86 - pkin(8) * t85;
t48 = pkin(5) * t88 - pkin(8) * t87 + qJD(3);
t44 = t86 * Ifges(6,1) + Ifges(6,5) * qJD(5) + t80;
t43 = t85 * Ifges(6,2) + Ifges(6,6) * qJD(5) + t233;
t32 = qJD(4) * t153 + qJD(5) * t157;
t23 = Ifges(7,1) * t64 + Ifges(7,5) * t261 + t62;
t21 = t64 * Ifges(7,5) + t63 * Ifges(7,6) + Ifges(7,3) * t261;
t17 = t140 * t52 + t143 * t30;
t16 = -t140 * t30 + t143 * t52;
t5 = t26 * Ifges(7,4) + t27 * Ifges(7,2) + t51 * Ifges(7,6);
t4 = -qJD(6) * t20 - t140 * t32 + t143 * t48;
t3 = qJD(6) * t19 + t140 * t48 + t143 * t32;
t6 = [t170 * t237 + (-m(6) * t30 - t253) * (qJD(4) * t152 + qJD(5) * t42) + t261 * (Ifges(7,5) * t155 + Ifges(7,6) * t156 + Ifges(7,3) * t88) / 0.2e1 - (-t75 * mrSges(6,2) - Ifges(6,1) * t53 - Ifges(6,4) * t54 - Ifges(6,5) * qJDD(5) - t167 * t245 - t168 * t246 - t169 * t247 + t178 * t22) * t152 - (t189 / 0.2e1 - Ifges(6,2) * t54 - Ifges(6,6) * qJDD(5) - Ifges(6,4) * t53 + t75 * mrSges(6,1) + Ifges(7,5) * t247 + Ifges(7,3) * t245 + Ifges(7,6) * t246 + t251) * t153 + 0.2e1 * t138 * mrSges(3,1) * t206 + t23 * t210 / 0.2e1 + t110 * t199 - (Ifges(5,4) * t137 - Ifges(5,2) * t135) * t191 - 0.2e1 * mrSges(3,2) * t117 + t252 * mrSges(6,3) + t248 * qJD(3) + t97 * t180 + t89 * t174 + t63 * (Ifges(7,4) * t155 + Ifges(7,2) * t156 + Ifges(7,6) * t88) / 0.2e1 + t28 * (-mrSges(7,1) * t156 + mrSges(7,2) * t155) + (-t91 + t133) * mrSges(4,3) + (-t1 * t213 - t12 * t155 + t13 * t156 - t2 * t209) * mrSges(7,3) + (t266 * t267 - t164) * mrSges(5,3) - (qJD(6) * t23 + t5) * t213 / 0.2e1 - (-m(6) * t11 + m(7) * t8 + t236) * t157 + (m(3) * t230 - m(4) * t161 - m(5) * t177 + mrSges(2,1) * t142 - t72 * mrSges(7,1) + mrSges(2,2) * t144 + t71 * mrSges(7,2) + t264 * (t123 * t139 + t125 * t126 + t161) + (-m(5) * (-pkin(2) - qJ(4)) + t256) * t123 + t250 * t125) * g(1) + (-m(3) * t127 - mrSges(2,1) * t144 - t70 * mrSges(7,1) + mrSges(2,2) * t142 - t69 * mrSges(7,2) + (-m(5) - m(4)) * t184 + t264 * (t123 * t126 - t125 * t139 + t184) + (-m(5) * qJ(4) - t256) * t125 + t250 * t123) * g(2) - t88 * t262 + t88 * t263 + (t110 * mrSges(4,3) + t115 * mrSges(4,2) + m(3) * (t136 ^ 2 + t138 ^ 2) * pkin(1) ^ 2 + Ifges(4,1) + Ifges(3,3) + Ifges(2,3)) * qJDD(1) + m(6) * (t10 * t42 + t31 * t32 + t75 * t97) + m(7) * (t1 * t20 + t12 * t4 + t13 * t3 + t19 * t2) + m(4) * (-t110 * t91 + t115 * t95) + m(5) * (qJD(4) * t163 + t107 * t164 + t110 * t89) + t95 * mrSges(4,2) + t87 * t44 / 0.2e1 + t88 * t21 / 0.2e1 + t85 * (Ifges(6,4) * t87 - Ifges(6,2) * t88) / 0.2e1 + t82 * (mrSges(6,1) * t88 + mrSges(6,2) * t87) + qJD(5) * (Ifges(6,5) * t87 - Ifges(6,6) * t88) / 0.2e1 - t88 * t43 / 0.2e1 + t32 * t73 + t42 * t40 + t3 * t36 + t4 * t37 + t19 * t14 + t20 * t15 + t87 * t186 + t209 * t265 + (Ifges(5,1) * t137 - Ifges(5,4) * t135) * t190 + (Ifges(6,1) * t87 - Ifges(6,4) * t88) * t238 + (Ifges(7,1) * t155 + Ifges(7,4) * t156 + Ifges(7,5) * t88) * t242; -t236 * t153 + t218 * t88 + (m(4) + m(3)) * qJDD(2) + t159 * t87 - t255 * t152 + m(5) * (-t135 * t60 + t137 * t61) + m(6) * (t10 * t152 + t11 * t153 - t30 * t88 + t31 * t87) + m(7) * (-t152 * t254 - t153 * t8 + t165 * t87 + t28 * t88) + (-m(3) - t188) * g(3); -t145 * mrSges(4,3) - t236 * t152 - t218 * t87 + t159 * t88 + (mrSges(4,2) - t176) * qJDD(1) + t255 * t153 + m(7) * (t153 * t254 + t165 * t88 - t28 * t87 - t237) + m(5) * t164 - m(6) * t252 + m(4) * t95 + (-m(7) * t166 - t162 - t248) * qJD(1) + t258 * t188; -t85 * t73 - t218 * t86 - t145 * t176 + (t261 * t36 + t14) * t143 + (-t261 * t37 + t15) * t140 + t180 + t199 + (g(1) * t125 + g(2) * t123) * t194 + (t1 * t140 + t143 * t2 + t165 * t261 - t28 * t86) * m(7) + (t30 * t86 - t31 * t85 + t75) * m(6) + (-qJD(1) * t163 + t89) * m(5); ((-t196 + t214) * t13 + (-t195 + t211) * t12 - t175) * mrSges(7,3) + (t167 * t261 + t168 * t63 + t169 * t64) * qJD(6) / 0.2e1 + t22 * t214 / 0.2e1 + (t221 + t253) * t31 + t8 * t171 - t85 * t154 + (-t73 + t234) * t30 + ((mrSges(6,1) + t150) * t124 + (-mrSges(6,2) + t187) * t122) * t258 + (-t211 / 0.2e1 + t178) * t23 - (Ifges(6,1) * t85 + t21 - t233) * t86 / 0.2e1 - (-Ifges(6,2) * t86 + t44 + t80) * t85 / 0.2e1 - t86 * t263 + (t186 + t154) * qJD(6) + (t122 * t150 - t124 * t187 + t172) * g(3) + t86 * t262 - t82 * (mrSges(6,1) * t86 + mrSges(6,2) * t85) - qJD(5) * (Ifges(6,5) * t85 - Ifges(6,6) * t86) / 0.2e1 + Ifges(6,5) * t53 + Ifges(6,6) * t54 - t17 * t36 - t16 * t37 - t10 * mrSges(6,2) + t11 * mrSges(6,1) - pkin(5) * t9 + (-pkin(5) * t8 - t12 * t16 - t13 * t17) * m(7) + (Ifges(7,1) * t140 + t216) * t247 + (-m(7) * t254 - t195 * t37 - t196 * t36 + t259) * pkin(8) + Ifges(6,3) * qJDD(5) + t140 * t265 + t143 * t5 / 0.2e1 + t43 * t238 + (Ifges(7,3) * t86 + t167 * t85) * t241 + (Ifges(7,5) * t86 + t169 * t85) * t243 + (Ifges(7,6) * t86 + t168 * t85) * t244 + (Ifges(7,5) * t140 + Ifges(7,6) * t143) * t245 + (Ifges(7,2) * t143 + t217) * t246; -t28 * (mrSges(7,1) * t64 + mrSges(7,2) * t63) + (Ifges(7,1) * t63 - t232) * t243 + t22 * t242 + (Ifges(7,5) * t63 - Ifges(7,6) * t64) * t241 - t12 * t36 + t13 * t37 - g(1) * (mrSges(7,1) * t69 - mrSges(7,2) * t70) - g(2) * (mrSges(7,1) * t71 + mrSges(7,2) * t72) + g(3) * t170 * t124 + (t12 * t63 + t13 * t64) * mrSges(7,3) + t189 + (-Ifges(7,2) * t64 + t23 + t62) * t244 + t251;];
tau  = t6;
