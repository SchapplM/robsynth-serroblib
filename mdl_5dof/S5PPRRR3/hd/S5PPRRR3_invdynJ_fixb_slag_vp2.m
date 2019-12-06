% Calculate vector of inverse dynamics joint torques for
% S5PPRRR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:23
% EndTime: 2019-12-05 15:16:32
% DurationCPUTime: 3.57s
% Computational Cost: add. (1931->323), mult. (4478->474), div. (0->0), fcn. (3250->12), ass. (0->158)
t120 = sin(qJ(4));
t123 = cos(qJ(4));
t117 = cos(pkin(9));
t159 = qJDD(1) * t117;
t121 = sin(qJ(3));
t124 = cos(qJ(3));
t115 = sin(pkin(9));
t168 = qJD(3) * t115;
t215 = -qJD(1) * t168 + qJDD(2);
t216 = qJD(2) * qJD(3) + qJDD(1) * t115;
t54 = t121 * t215 + t124 * t216;
t48 = qJDD(3) * pkin(6) + t54;
t169 = qJD(1) * t117;
t149 = t123 * t169;
t170 = qJD(1) * t115;
t84 = qJD(2) * t121 + t124 * t170;
t76 = qJD(3) * pkin(6) + t84;
t56 = -t120 * t76 - t149;
t14 = qJD(4) * t56 - t120 * t159 + t123 * t48;
t100 = t120 * t169;
t15 = -t120 * t48 + qJD(4) * t100 + (-qJD(4) * t76 - t159) * t123;
t139 = -t120 * t15 + t123 * t14;
t163 = qJD(4) * t123;
t164 = qJD(4) * t120;
t57 = t123 * t76 - t100;
t217 = -t56 * t163 - t57 * t164 + t139;
t161 = qJD(3) * qJD(4);
t90 = qJDD(3) * t123 - t120 * t161;
t214 = t90 / 0.2e1;
t119 = sin(qJ(5));
t122 = cos(qJ(5));
t125 = -pkin(7) - pkin(6);
t94 = t125 * t120;
t95 = t125 * t123;
t60 = t119 * t94 - t122 * t95;
t83 = qJD(2) * t124 - t121 * t170;
t86 = t119 * t123 + t120 * t122;
t150 = qJD(4) * t125;
t88 = t120 * t150;
t89 = t123 * t150;
t213 = -qJD(5) * t60 - t119 * t88 + t122 * t89 + t86 * t83;
t138 = t119 * t120 - t122 * t123;
t59 = t119 * t95 + t122 * t94;
t212 = qJD(5) * t59 + t119 * t89 + t122 * t88 + t138 * t83;
t202 = m(6) * pkin(4);
t211 = -mrSges(5,1) - t202;
t210 = -qJD(3) * t83 + t54;
t55 = -t121 * t216 + t124 * t215;
t209 = qJD(3) * t84 + t55;
t208 = pkin(4) * t164 - t84;
t67 = t138 * t121;
t167 = qJD(3) * t120;
t154 = mrSges(5,3) * t167;
t92 = qJD(4) * mrSges(5,1) - t154;
t166 = qJD(3) * t123;
t153 = mrSges(5,3) * t166;
t93 = -qJD(4) * mrSges(5,2) + t153;
t206 = -t120 * t92 + t123 * t93;
t68 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t90;
t91 = qJDD(3) * t120 + t123 * t161;
t69 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t91;
t205 = -t120 * t69 + t123 * t68;
t113 = qJD(4) + qJD(5);
t105 = pkin(4) * t123 + pkin(3);
t114 = qJ(4) + qJ(5);
t110 = sin(t114);
t111 = cos(t114);
t142 = -mrSges(5,1) * t123 + mrSges(5,2) * t120;
t204 = m(5) * pkin(3) + m(6) * t105 + mrSges(6,1) * t111 - mrSges(6,2) * t110 + mrSges(4,1) - t142;
t203 = -m(5) * pkin(6) + m(6) * t125 + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t126 = qJD(3) ^ 2;
t82 = t86 * qJD(3);
t200 = t82 / 0.2e1;
t144 = pkin(7) * qJD(3) + t76;
t47 = t123 * t144 - t100;
t189 = t119 * t47;
t46 = -t120 * t144 - t149;
t38 = qJD(4) * pkin(4) + t46;
t10 = t122 * t38 - t189;
t198 = t10 * mrSges(6,3);
t197 = t82 * mrSges(6,3);
t196 = t82 * Ifges(6,4);
t116 = sin(pkin(8));
t179 = t115 * t116;
t118 = cos(pkin(8));
t172 = t118 * t121;
t173 = t116 * t124;
t72 = t117 * t173 - t172;
t195 = (-t110 * t72 + t111 * t179) * mrSges(6,1) + (-t110 * t179 - t111 * t72) * mrSges(6,2);
t178 = t115 * t118;
t171 = t118 * t124;
t174 = t116 * t121;
t74 = t117 * t171 + t174;
t194 = (-t110 * t74 + t111 * t178) * mrSges(6,1) + (-t110 * t178 - t111 * t74) * mrSges(6,2);
t175 = t115 * t124;
t193 = (-t110 * t175 - t111 * t117) * mrSges(6,1) + (t110 * t117 - t111 * t175) * mrSges(6,2);
t192 = Ifges(5,4) * t120;
t191 = Ifges(5,4) * t123;
t185 = t122 * t47;
t177 = t115 * t120;
t176 = t115 * t123;
t165 = qJD(3) * t124;
t158 = pkin(4) * t167;
t156 = m(3) + m(4) + m(5) + m(6);
t148 = t121 * t168;
t141 = t123 * Ifges(5,2) + t192;
t140 = Ifges(5,5) * t123 - Ifges(5,6) * t120;
t11 = t119 * t38 + t185;
t77 = -t117 * t123 - t120 * t175;
t78 = -t117 * t120 + t123 * t175;
t34 = -t119 * t78 + t122 * t77;
t35 = t119 * t77 + t122 * t78;
t75 = -qJD(3) * pkin(3) - t83;
t135 = t75 * (mrSges(5,1) * t120 + mrSges(5,2) * t123);
t134 = t120 * (Ifges(5,1) * t123 - t192);
t58 = -mrSges(5,1) * t90 + mrSges(5,2) * t91;
t131 = t138 * qJD(5);
t28 = -qJD(3) * t131 + t119 * t90 + t122 * t91;
t132 = t86 * qJD(5);
t29 = -qJD(3) * t132 - t119 * t91 + t122 * t90;
t6 = -mrSges(6,1) * t29 + mrSges(6,2) * t28;
t133 = qJDD(3) * mrSges(4,1) - mrSges(4,2) * t126 - t58 - t6;
t81 = t138 * qJD(3);
t49 = -qJDD(3) * pkin(3) - t55;
t45 = mrSges(6,1) * t81 + mrSges(6,2) * t82;
t87 = t142 * qJD(3);
t128 = -mrSges(4,1) * t126 - qJDD(3) * mrSges(4,2) + (t45 + t87) * qJD(3);
t112 = qJDD(4) + qJDD(5);
t8 = qJDD(4) * pkin(4) - pkin(7) * t91 + t15;
t9 = pkin(7) * t90 + t14;
t2 = qJD(5) * t10 + t119 * t8 + t122 * t9;
t3 = -qJD(5) * t11 - t119 * t9 + t122 * t8;
t32 = -t81 * Ifges(6,2) + t113 * Ifges(6,6) + t196;
t70 = Ifges(6,4) * t81;
t33 = t82 * Ifges(6,1) + t113 * Ifges(6,5) - t70;
t63 = -qJD(3) * t105 - t83;
t127 = t3 * mrSges(6,1) - t2 * mrSges(6,2) - t81 * t198 + t32 * t200 - t63 * (mrSges(6,1) * t82 - mrSges(6,2) * t81) + Ifges(6,3) * t112 - t82 * (-Ifges(6,1) * t81 - t196) / 0.2e1 + Ifges(6,6) * t29 + Ifges(6,5) * t28 - t113 * (-Ifges(6,5) * t81 - Ifges(6,6) * t82) / 0.2e1 + (-Ifges(6,2) * t82 + t33 - t70) * t81 / 0.2e1;
t51 = -qJD(4) * t86 - t132;
t108 = t117 ^ 2 * qJDD(1);
t106 = Ifges(5,4) * t166;
t80 = Ifges(5,1) * t167 + Ifges(5,5) * qJD(4) + t106;
t79 = Ifges(5,6) * qJD(4) + qJD(3) * t141;
t66 = t86 * t121;
t62 = mrSges(6,1) * t113 - t197;
t61 = -mrSges(6,2) * t113 - mrSges(6,3) * t81;
t53 = qJD(4) * t77 - t123 * t148;
t52 = -qJD(4) * t78 + t120 * t148;
t50 = -qJD(4) * t138 - t131;
t30 = -pkin(4) * t90 + t49;
t23 = -mrSges(6,2) * t112 + mrSges(6,3) * t29;
t22 = mrSges(6,1) * t112 - mrSges(6,3) * t28;
t17 = t113 * t67 - t86 * t165;
t16 = t121 * t51 - t124 * t81;
t13 = t122 * t46 - t189;
t12 = -t119 * t46 - t185;
t5 = -qJD(5) * t35 - t119 * t53 + t122 * t52;
t4 = qJD(5) * t34 + t119 * t52 + t122 * t53;
t1 = [m(2) * qJDD(1) + t34 * t22 + t35 * t23 + t4 * t61 + t5 * t62 + t52 * t92 + t53 * t93 + t78 * t68 + t77 * t69 + (-m(2) - t156) * g(3) + (-t121 * t133 + t124 * t128) * t115 + m(4) * (t108 + (-t121 * t55 + t124 * t54 + (-t121 * t84 - t124 * t83) * qJD(3)) * t115) + m(5) * (t14 * t78 + t15 * t77 + t52 * t56 + t53 * t57 + (t121 * t49 + t165 * t75) * t115) + m(6) * (t10 * t5 + t11 * t4 + t2 * t35 + t3 * t34 + (t121 * t30 + t165 * t63) * t115) + m(3) * (qJDD(1) * t115 ^ 2 + t108); -t66 * t22 - t67 * t23 + t16 * t61 + t17 * t62 + m(6) * (t10 * t17 + t11 * t16 - t2 * t67 - t3 * t66) + m(3) * qJDD(2) + (t206 * qJD(3) + m(5) * (t166 * t57 - t167 * t56 - t49) + m(4) * t209 - m(6) * t30 + t133) * t124 + ((-t120 * t93 - t123 * t92) * qJD(4) + m(5) * (qJD(3) * t75 + t217) + m(4) * t210 + m(6) * qJD(3) * t63 + t128 + t205) * t121 + (-g(1) * t116 + g(2) * t118) * t156; (t121 * t204 + t124 * t203) * g(3) * t115 + t123 * (Ifges(5,4) * t91 + Ifges(5,2) * t90) / 0.2e1 + (Ifges(5,1) * t91 + Ifges(5,4) * t214) * t120 + (-pkin(3) * t49 - t75 * t84 - (-t120 * t56 + t123 * t57) * t83) * m(5) + t212 * t61 + (t10 * t213 - t105 * t30 + t11 * t212 + t2 * t60 + t208 * t63 + t3 * t59) * m(6) + t213 * t62 + t208 * t45 + t209 * mrSges(4,1) - t210 * mrSges(4,2) + (mrSges(6,2) * t30 - mrSges(6,3) * t3 + Ifges(6,1) * t28 + Ifges(6,4) * t29 + Ifges(6,5) * t112) * t86 + (t135 + t140 * qJD(4) / 0.2e1) * qJD(4) + (m(5) * ((-t120 * t57 - t123 * t56) * qJD(4) + t139) - t92 * t163 - t93 * t164 + t205) * pkin(6) - t206 * t83 - t50 * t198 + qJDD(4) * (Ifges(5,5) * t120 + Ifges(5,6) * t123) + t113 * (Ifges(6,5) * t50 + Ifges(6,6) * t51) / 0.2e1 - t105 * t6 - t84 * t87 - t81 * (Ifges(6,4) * t50 + Ifges(6,2) * t51) / 0.2e1 - pkin(3) * t58 + t59 * t22 + t60 * t23 + t63 * (-mrSges(6,1) * t51 + mrSges(6,2) * t50) + t50 * t33 / 0.2e1 + t51 * t32 / 0.2e1 + (t203 * t72 - t204 * (-t117 * t174 - t171)) * g(2) + (t203 * t74 - t204 * (-t117 * t172 + t173)) * g(1) + t11 * t51 * mrSges(6,3) - (-mrSges(6,1) * t30 + mrSges(6,3) * t2 + Ifges(6,4) * t28 + Ifges(6,2) * t29 + Ifges(6,6) * t112) * t138 + (t123 * (-Ifges(5,2) * t120 + t191) + t134) * t161 / 0.2e1 + t217 * mrSges(5,3) + t80 * t163 / 0.2e1 - t79 * t164 / 0.2e1 + t141 * t214 + t91 * t191 / 0.2e1 + t49 * t142 + (Ifges(6,1) * t50 + Ifges(6,4) * t51) * t200 + Ifges(4,3) * qJDD(3); -t126 * t134 / 0.2e1 + t127 + Ifges(5,6) * t90 + Ifges(5,5) * t91 - t13 * t61 - t12 * t62 - t14 * mrSges(5,2) + t15 * mrSges(5,1) + t79 * t167 / 0.2e1 - t140 * t161 / 0.2e1 - t45 * t158 - m(6) * (t10 * t12 + t11 * t13 + t158 * t63) - qJD(3) * t135 + t11 * t197 + (t119 * t2 + t122 * t3 + (-t10 * t119 + t11 * t122) * qJD(5)) * t202 + Ifges(5,3) * qJDD(4) + (t92 + t154) * t57 + (-t93 + t153) * t56 - (-Ifges(5,2) * t167 + t106 + t80) * t166 / 0.2e1 + (mrSges(5,2) * t78 + t211 * t77 - t193) * g(3) + (-t195 - (-t116 * t177 - t123 * t72) * mrSges(5,2) + t211 * (t116 * t176 - t120 * t72)) * g(2) + (-t194 - (-t118 * t177 - t123 * t74) * mrSges(5,2) + t211 * (t118 * t176 - t120 * t74)) * g(1) + ((-t119 * t62 + t122 * t61) * qJD(5) + t119 * t23 + t122 * t22) * pkin(4); t127 + (t62 + t197) * t11 - t10 * t61 - g(3) * t193 - g(2) * t195 - g(1) * t194;];
tau = t1;
