% Calculate vector of inverse dynamics joint torques for
% S5PRRPP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:12
% EndTime: 2019-12-05 16:06:28
% DurationCPUTime: 6.06s
% Computational Cost: add. (1526->335), mult. (3293->425), div. (0->0), fcn. (2002->8), ass. (0->143)
t186 = m(5) + m(6);
t97 = sin(qJ(3));
t158 = Ifges(4,4) * t97;
t190 = mrSges(5,1) + mrSges(6,1);
t189 = mrSges(5,2) - mrSges(6,3);
t185 = Ifges(5,1) + Ifges(6,1);
t183 = Ifges(6,4) + Ifges(5,5);
t184 = -Ifges(5,4) + Ifges(6,5);
t98 = cos(qJ(3));
t76 = -mrSges(4,1) * t98 + t97 * mrSges(4,2);
t94 = qJ(3) + pkin(8);
t87 = sin(t94);
t89 = cos(t94);
t188 = t189 * t87 - t190 * t89 + t76;
t93 = pkin(7) + qJ(2);
t86 = sin(t93);
t159 = g(2) * t86;
t88 = cos(t93);
t187 = g(1) * t88 + t159;
t142 = cos(pkin(8));
t116 = t142 * t97;
t95 = sin(pkin(8));
t66 = t95 * t98 + t116;
t165 = t66 / 0.2e1;
t182 = Ifges(6,6) - Ifges(5,6);
t115 = t142 * t98;
t141 = qJD(2) * t97;
t58 = -qJD(2) * t115 + t141 * t95;
t155 = Ifges(6,5) * t58;
t52 = Ifges(5,4) * t58;
t60 = t66 * qJD(2);
t181 = t183 * qJD(3) + t185 * t60 + t155 - t52;
t132 = qJD(2) * qJD(3);
t121 = t97 * t132;
t134 = qJDD(2) * t98;
t68 = -t121 + t134;
t122 = t98 * t132;
t69 = t97 * qJDD(2) + t122;
t35 = -t142 * t68 + t69 * t95;
t25 = -qJDD(3) * mrSges(5,2) - mrSges(5,3) * t35;
t28 = -mrSges(6,2) * t35 + qJDD(3) * mrSges(6,3);
t180 = t25 + t28;
t36 = t142 * t69 + t95 * t68;
t26 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t36;
t27 = -qJDD(3) * mrSges(6,1) + t36 * mrSges(6,2);
t179 = -t26 + t27;
t150 = t58 * mrSges(5,3);
t151 = t58 * mrSges(6,2);
t45 = qJD(3) * mrSges(6,3) - t151;
t145 = -qJD(3) * mrSges(5,2) - t150 + t45;
t149 = t60 * mrSges(6,2);
t144 = -mrSges(5,3) * t60 + qJD(3) * t190 - t149;
t136 = qJDD(2) * pkin(2);
t133 = qJD(1) * qJD(3);
t131 = pkin(6) * t134 + t97 * qJDD(1) + t98 * t133;
t39 = -pkin(6) * t121 + t131;
t91 = t98 * qJDD(1);
t40 = -pkin(6) * t69 - t133 * t97 + t91;
t177 = t39 * t98 - t40 * t97;
t175 = -m(4) * pkin(6) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t173 = m(4) * pkin(2) + mrSges(3,1) - t188;
t170 = -t58 / 0.2e1;
t169 = t58 / 0.2e1;
t167 = t60 / 0.2e1;
t164 = m(2) + m(3);
t163 = pkin(3) * t95;
t161 = pkin(3) * t98;
t11 = -pkin(6) * t122 + qJDD(3) * pkin(3) - t69 * qJ(4) + t91 + (-pkin(6) * qJDD(2) - qJD(2) * qJD(4) - t133) * t97;
t139 = qJD(3) * t97;
t128 = pkin(6) * t139;
t137 = qJD(4) * t98;
t15 = t68 * qJ(4) + (-t128 + t137) * qJD(2) + t131;
t5 = t95 * t11 + t142 * t15;
t157 = Ifges(4,4) * t98;
t156 = Ifges(5,4) * t60;
t135 = t97 * qJD(1);
t96 = -qJ(4) - pkin(6);
t77 = t96 * t98;
t56 = -qJD(2) * t77 + t135;
t46 = t142 * t56;
t124 = t96 * t97;
t92 = t98 * qJD(1);
t55 = qJD(2) * t124 + t92;
t49 = qJD(3) * pkin(3) + t55;
t14 = t95 * t49 + t46;
t154 = t14 * mrSges(5,3);
t148 = t95 * t56;
t146 = qJD(3) / 0.2e1;
t140 = qJD(2) * t98;
t138 = qJD(3) * t98;
t130 = pkin(3) * t141;
t129 = pkin(3) * t139;
t127 = mrSges(4,3) * t141;
t126 = mrSges(4,3) * t140;
t83 = pkin(2) + t161;
t123 = t142 * pkin(3);
t119 = t35 * mrSges(5,1) + t36 * mrSges(5,2);
t118 = t35 * mrSges(6,1) - t36 * mrSges(6,3);
t114 = qJD(3) * t96;
t112 = -g(1) * t86 + g(2) * t88;
t111 = mrSges(4,1) * t97 + mrSges(4,2) * t98;
t108 = t98 * Ifges(4,2) + t158;
t107 = Ifges(4,5) * t98 - Ifges(4,6) * t97;
t106 = t89 * pkin(4) + t87 * qJ(5);
t105 = pkin(2) * t111;
t104 = t97 * (Ifges(4,1) * t98 - t158);
t50 = -pkin(3) * t68 + qJDD(4) - t136;
t4 = t11 * t142 - t95 * t15;
t13 = t142 * t49 - t148;
t103 = -t95 * t97 + t115;
t73 = -qJD(2) * t83 + qJD(4);
t101 = -t97 * qJD(4) + t114 * t98;
t85 = Ifges(4,4) * t140;
t81 = -t123 - pkin(4);
t79 = qJ(5) + t163;
t75 = -qJD(3) * mrSges(4,2) + t126;
t74 = qJD(3) * mrSges(4,1) - t127;
t72 = pkin(6) * t140 + t135;
t71 = -pkin(6) * t141 + t92;
t63 = Ifges(4,1) * t141 + Ifges(4,5) * qJD(3) + t85;
t62 = Ifges(4,6) * qJD(3) + qJD(2) * t108;
t61 = t103 * qJD(3);
t59 = t66 * qJD(3);
t57 = t114 * t97 + t137;
t54 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t69;
t53 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t68;
t51 = Ifges(6,5) * t60;
t31 = -pkin(4) * t103 - t66 * qJ(5) - t83;
t30 = mrSges(5,1) * t58 + mrSges(5,2) * t60;
t29 = mrSges(6,1) * t58 - mrSges(6,3) * t60;
t22 = -Ifges(5,2) * t58 + Ifges(5,6) * qJD(3) + t156;
t21 = Ifges(6,6) * qJD(3) + Ifges(6,3) * t58 + t51;
t19 = t142 * t55 - t148;
t17 = t55 * t95 + t46;
t16 = pkin(4) * t60 + qJ(5) * t58 + t130;
t12 = t58 * pkin(4) - t60 * qJ(5) + t73;
t9 = qJD(3) * qJ(5) + t14;
t8 = -qJD(3) * pkin(4) + qJD(5) - t13;
t6 = pkin(4) * t59 - qJ(5) * t61 - qJD(5) * t66 + t129;
t3 = pkin(4) * t35 - qJ(5) * t36 - qJD(5) * t60 + t50;
t2 = -qJDD(3) * pkin(4) + qJDD(5) - t4;
t1 = qJDD(3) * qJ(5) + qJD(3) * qJD(5) + t5;
t7 = [t97 * t53 + t98 * t54 + t180 * t66 - t179 * t103 + t145 * t61 - t144 * t59 + t164 * qJDD(1) + (-t74 * t97 + t75 * t98) * qJD(3) + m(4) * (t39 * t97 + t40 * t98 + (-t71 * t97 + t72 * t98) * qJD(3)) + m(5) * (t103 * t4 - t13 * t59 + t14 * t61 + t5 * t66) + m(6) * (t1 * t66 - t103 * t2 + t59 * t8 + t61 * t9) + (-m(4) - t164 - t186) * g(3); t69 * Ifges(4,1) * t97 + (-(Ifges(5,2) + Ifges(6,3)) * t103 + 0.2e1 * t184 * t165) * t35 - t75 * t128 + t103 * (Ifges(5,4) * t36 + Ifges(5,6) * qJDD(3)) / 0.2e1 + t50 * (-mrSges(5,1) * t103 + mrSges(5,2) * t66) + t3 * (-mrSges(6,1) * t103 - mrSges(6,3) * t66) - t103 * (Ifges(6,5) * t36 + Ifges(6,6) * qJDD(3)) / 0.2e1 + (-t103 * t184 + t185 * t66) * t36 / 0.2e1 + t69 * t157 / 0.2e1 + (-t103 * t182 + t183 * t66) * qJDD(3) / 0.2e1 + (t158 + t108) * t68 / 0.2e1 + (-t13 * mrSges(5,3) + t181 / 0.2e1 + mrSges(5,2) * t73 - mrSges(6,3) * t12 + Ifges(5,4) * t170 + Ifges(6,5) * t169 + t183 * t146 + t185 * t167 + t8 * mrSges(6,2)) * t61 + (t98 * (-Ifges(4,2) * t97 + t157) + t104) * t132 / 0.2e1 + (t103 * t5 - t4 * t66) * mrSges(5,3) + t31 * t118 - t83 * t119 + t97 * Ifges(4,5) * qJDD(3) - t62 * t139 / 0.2e1 - t105 * t132 - t76 * t136 + t63 * t138 / 0.2e1 + t30 * t129 + qJD(3) ^ 2 * t107 / 0.2e1 + t6 * t29 + (t1 * t103 + t2 * t66 - t59 * t9 - t159) * mrSges(6,2) + m(6) * (t12 * t6 + t3 * t31) + m(5) * (t129 * t73 - t50 * t83) + Ifges(3,3) * qJDD(2) + (-t74 * t138 + m(4) * ((-t71 * t98 - t72 * t97) * qJD(3) + t177) + t98 * t53 - t97 * t54) * pkin(6) + (-t138 * t71 - t139 * t72 + t177) * mrSges(4,3) + (m(4) * t136 + mrSges(4,1) * t68 - mrSges(4,2) * t69) * pkin(2) + (-m(5) * t13 + m(6) * t8 - t144) * (-t101 * t142 + t57 * t95) + (m(5) * t14 + m(6) * t9 + t145) * (t101 * t95 + t142 * t57) + (-m(5) * t4 + m(6) * t2 + t179) * (-t116 * t96 - t77 * t95) + (m(5) * t5 + m(6) * t1 + t180) * (t124 * t95 - t142 * t77) + (-t154 + Ifges(6,3) * t169 - Ifges(5,2) * t170 + t73 * mrSges(5,1) + t12 * mrSges(6,1) + t21 / 0.2e1 - t22 / 0.2e1 + t184 * t167 + t182 * t146) * t59 + (t183 * qJDD(3) + t185 * t36) * t165 + ((t186 * t96 - mrSges(6,2) + t175) * t88 + (m(5) * t83 - m(6) * (-t106 - t83) + t173) * t86) * g(1) + (-t186 * (t88 * t83 - t86 * t96) + t175 * t86 + (-m(6) * t106 - t173) * t88) * g(2) + t98 * qJDD(3) * Ifges(4,6) + t98 * (Ifges(4,4) * t69 + Ifges(4,2) * t68) / 0.2e1; (t13 * t17 - t73 * t130 - t14 * t19 + (t142 * t4 + t5 * t95) * pkin(3)) * m(5) - t30 * t130 + t187 * (t111 + t186 * pkin(3) * t97 + (-m(6) * qJ(5) + t189) * t89 + (m(6) * pkin(4) + t190) * t87) - (-Ifges(4,2) * t141 + t63 + t85) * t140 / 0.2e1 + (-t104 / 0.2e1 + t105) * qJD(2) ^ 2 + (t74 + t127) * t72 + (t1 * t79 - t12 * t16 - t17 * t8 + t2 * t81 + (-t19 + qJD(5)) * t9) * m(6) + (Ifges(4,3) + Ifges(6,2) + Ifges(5,3)) * qJDD(3) + (-t75 + t126) * t71 - t13 * t150 + t25 * t163 + t22 * t167 + (Ifges(6,3) * t60 - t155) * t170 + t9 * t149 + t8 * t151 + t60 * t154 + t26 * t123 + t62 * t141 / 0.2e1 - t107 * t132 / 0.2e1 + t79 * t28 + t81 * t27 + Ifges(4,6) * t68 + Ifges(4,5) * t69 - t73 * (mrSges(5,1) * t60 - mrSges(5,2) * t58) - t12 * (mrSges(6,1) * t60 + mrSges(6,3) * t58) + qJD(5) * t45 - t39 * mrSges(4,2) + t40 * mrSges(4,1) - t16 * t29 + t4 * mrSges(5,1) - t5 * mrSges(5,2) + t1 * mrSges(6,3) - t2 * mrSges(6,1) + (-m(5) * t161 - m(6) * (t106 + t161) + t188) * g(3) + t144 * t17 - t145 * t19 + (-Ifges(5,2) * t60 + t181 - t52) * t169 + t182 * t35 + t183 * t36 - (t182 * t60 - t183 * t58) * qJD(3) / 0.2e1 - (-t185 * t58 - t156 + t21 + t51) * t60 / 0.2e1; t144 * t60 + t145 * t58 + t118 + t119 + (t58 * t9 - t60 * t8 + t112 + t3) * m(6) + (t13 * t60 + t14 * t58 + t112 + t50) * m(5); -qJD(3) * t45 + t60 * t29 + (g(3) * t89 - t9 * qJD(3) + t12 * t60 - t187 * t87 + t2) * m(6) + t27;];
tau = t7;
