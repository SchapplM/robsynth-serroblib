% Calculate vector of inverse dynamics joint torques for
% S5RPPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:41:38
% EndTime: 2019-12-05 17:41:49
% DurationCPUTime: 4.56s
% Computational Cost: add. (3080->315), mult. (6900->425), div. (0->0), fcn. (4932->16), ass. (0->139)
t127 = sin(pkin(8));
t102 = pkin(1) * t127 + qJ(3);
t90 = qJD(1) * qJD(3) + qJDD(1) * t102;
t119 = qJDD(4) + qJDD(5);
t124 = qJD(4) + qJD(5);
t131 = sin(qJ(5));
t134 = cos(qJ(5));
t126 = sin(pkin(9));
t132 = sin(qJ(4));
t128 = cos(pkin(9));
t135 = cos(qJ(4));
t161 = t128 * t135;
t92 = -t126 * t132 + t161;
t84 = t92 * qJD(1);
t93 = t126 * t135 + t128 * t132;
t85 = t93 * qJD(1);
t151 = -t131 * t85 + t134 * t84;
t50 = t131 * t84 + t134 * t85;
t178 = Ifges(6,4) * t50;
t86 = t92 * qJD(4);
t55 = qJD(1) * t86 + qJDD(1) * t93;
t87 = t93 * qJD(4);
t56 = -qJD(1) * t87 + qJDD(1) * t92;
t19 = qJD(5) * t151 + t131 * t56 + t134 * t55;
t112 = t128 * qJD(2);
t170 = pkin(6) * qJD(1);
t99 = t102 * qJD(1);
t67 = t112 + (-t99 - t170) * t126;
t75 = t126 * qJD(2) + t128 * t99;
t68 = t128 * t170 + t75;
t34 = t132 * t67 + t135 * t68;
t110 = t128 * qJDD(2);
t63 = t110 + (-pkin(6) * qJDD(1) - t90) * t126;
t156 = qJDD(1) * t128;
t70 = t126 * qJDD(2) + t128 * t90;
t64 = pkin(6) * t156 + t70;
t13 = -qJD(4) * t34 - t132 * t64 + t135 * t63;
t6 = qJDD(4) * pkin(4) - pkin(7) * t55 + t13;
t160 = qJD(4) * t135;
t166 = t132 * t68;
t12 = -qJD(4) * t166 + t132 * t63 + t135 * t64 + t67 * t160;
t7 = pkin(7) * t56 + t12;
t28 = pkin(7) * t84 + t34;
t168 = t131 * t28;
t33 = t135 * t67 - t166;
t27 = -pkin(7) * t85 + t33;
t26 = qJD(4) * pkin(4) + t27;
t8 = t134 * t26 - t168;
t2 = qJD(5) * t8 + t131 * t6 + t134 * t7;
t20 = -qJD(5) * t50 - t131 * t55 + t134 * t56;
t42 = Ifges(6,4) * t151;
t24 = Ifges(6,1) * t50 + Ifges(6,5) * t124 + t42;
t165 = t134 * t28;
t9 = t131 * t26 + t165;
t3 = -qJD(5) * t9 - t131 * t7 + t134 * t6;
t129 = cos(pkin(8));
t107 = -pkin(1) * t129 - pkin(2);
t118 = t128 * pkin(3);
t193 = t107 - t118;
t83 = qJD(1) * t193 + qJD(3);
t54 = -pkin(4) * t84 + t83;
t205 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t19 + Ifges(6,6) * t20 + Ifges(6,3) * t119 - (Ifges(6,5) * t151 - Ifges(6,6) * t50) * t124 / 0.2e1 + (t151 * t8 + t50 * t9) * mrSges(6,3) - (-Ifges(6,2) * t50 + t24 + t42) * t151 / 0.2e1 - t54 * (mrSges(6,1) * t50 + mrSges(6,2) * t151) - (Ifges(6,1) * t151 - t178) * t50 / 0.2e1;
t159 = -m(4) - m(5) - m(6);
t197 = -m(3) + t159;
t204 = -pkin(1) * t197 + mrSges(2,1);
t203 = t126 ^ 2 + t128 ^ 2;
t123 = pkin(9) + qJ(4);
t113 = sin(t123);
t115 = cos(t123);
t117 = qJ(5) + t123;
t104 = sin(t117);
t105 = cos(t117);
t150 = t105 * mrSges(6,1) - t104 * mrSges(6,2);
t202 = -t115 * mrSges(5,1) + t113 * mrSges(5,2) - t150;
t201 = mrSges(6,1) * t104 + mrSges(6,2) * t105;
t23 = Ifges(6,2) * t151 + Ifges(6,6) * t124 + t178;
t199 = t23 / 0.2e1;
t69 = -t126 * t90 + t110;
t191 = -t69 * t126 + t70 * t128;
t106 = t118 + pkin(2);
t148 = -mrSges(4,1) * t128 + mrSges(4,2) * t126;
t189 = m(5) * t106 + mrSges(3,1) + m(6) * (pkin(4) * t115 + t106) + m(4) * pkin(2) - t148 - t202;
t130 = -pkin(6) - qJ(3);
t187 = m(4) * qJ(3) - m(5) * t130 - m(6) * (-pkin(7) + t130) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t184 = t50 / 0.2e1;
t182 = t85 / 0.2e1;
t181 = pkin(4) * t87;
t179 = Ifges(5,4) * t85;
t125 = qJ(1) + pkin(8);
t116 = cos(t125);
t175 = g(3) * t116;
t174 = pkin(6) + t102;
t88 = t174 * t126;
t89 = t174 * t128;
t52 = -t132 * t88 + t135 * t89;
t114 = sin(t125);
t173 = t201 * t114;
t157 = qJDD(1) * t126;
t153 = -t56 * mrSges(5,1) + t55 * mrSges(5,2);
t152 = -t20 * mrSges(6,1) + t19 * mrSges(6,2);
t51 = -t132 * t89 - t135 * t88;
t149 = -mrSges(4,1) * t156 + mrSges(4,2) * t157;
t146 = -mrSges(5,1) * t113 - mrSges(5,2) * t115;
t144 = -t126 * (-t126 * t99 + t112) + t128 * t75;
t38 = -pkin(7) * t93 + t51;
t39 = pkin(7) * t92 + t52;
t21 = -t131 * t39 + t134 * t38;
t22 = t131 * t38 + t134 * t39;
t57 = -t131 * t93 + t134 * t92;
t58 = t131 * t92 + t134 * t93;
t81 = qJDD(1) * t193 + qJDD(3);
t35 = qJD(3) * t161 - t88 * t160 + (-qJD(3) * t126 - qJD(4) * t89) * t132;
t36 = -qJD(3) * t93 - qJD(4) * t52;
t136 = cos(qJ(1));
t133 = sin(qJ(1));
t97 = qJDD(1) * t107 + qJDD(3);
t80 = Ifges(5,4) * t84;
t72 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t85;
t71 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t84;
t66 = -pkin(4) * t92 + t193;
t46 = t85 * Ifges(5,1) + Ifges(5,5) * qJD(4) + t80;
t45 = t84 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t179;
t44 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t56;
t43 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t55;
t41 = mrSges(6,1) * t124 - mrSges(6,3) * t50;
t40 = -mrSges(6,2) * t124 + mrSges(6,3) * t151;
t37 = -pkin(4) * t56 + t81;
t32 = -pkin(7) * t86 + t36;
t31 = -pkin(7) * t87 + t35;
t30 = -qJD(5) * t58 - t131 * t86 - t134 * t87;
t29 = qJD(5) * t57 - t131 * t87 + t134 * t86;
t25 = -mrSges(6,1) * t151 + mrSges(6,2) * t50;
t15 = -mrSges(6,2) * t119 + mrSges(6,3) * t20;
t14 = mrSges(6,1) * t119 - mrSges(6,3) * t19;
t11 = t134 * t27 - t168;
t10 = -t131 * t27 - t165;
t5 = -qJD(5) * t22 - t131 * t31 + t134 * t32;
t4 = qJD(5) * t21 + t131 * t32 + t134 * t31;
t1 = [(-mrSges(5,1) * t81 + mrSges(5,3) * t12 + Ifges(5,4) * t55 + Ifges(5,2) * t56 + Ifges(5,6) * qJDD(4)) * t92 + (mrSges(5,2) * t81 - mrSges(5,3) * t13 + Ifges(5,1) * t55 + Ifges(5,4) * t56 + Ifges(5,5) * qJDD(4)) * t93 + t124 * (Ifges(6,5) * t29 + Ifges(6,6) * t30) / 0.2e1 + t86 * t46 / 0.2e1 - t87 * t45 / 0.2e1 + t35 * t71 + t36 * t72 + t51 * t43 + t52 * t44 + t54 * (-mrSges(6,1) * t30 + mrSges(6,2) * t29) + t4 * t40 + t5 * t41 + t29 * t24 / 0.2e1 + t21 * t14 + t22 * t15 + (-t29 * t8 + t30 * t9) * mrSges(6,3) + t25 * t181 + (-mrSges(2,2) * t133 + t114 * t187 + t116 * t189 + t136 * t204) * g(2) + (mrSges(2,2) * t136 + t114 * t189 - t116 * t187 + t133 * t204) * g(3) + m(4) * (t144 * qJD(3) + t102 * t191 + t107 * t97) + (-mrSges(6,1) * t37 + mrSges(6,3) * t2 + Ifges(6,4) * t19 + Ifges(6,2) * t20 + Ifges(6,6) * t119) * t57 + t30 * t199 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t129 - 0.2e1 * mrSges(3,2) * t127 + m(3) * (t127 ^ 2 + t129 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + m(5) * (t12 * t52 + t13 * t51 + t193 * t81 + t33 * t36 + t34 * t35) + t193 * t153 + (t203 * t90 + t191) * mrSges(4,3) + t97 * t148 + t107 * t149 + t66 * t152 + (Ifges(4,4) * t126 + Ifges(4,2) * t128) * t156 + (Ifges(4,1) * t126 + Ifges(4,4) * t128) * t157 + (Ifges(6,1) * t29 + Ifges(6,4) * t30) * t184 + t151 * (Ifges(6,4) * t29 + Ifges(6,2) * t30) / 0.2e1 + qJD(4) * (Ifges(5,5) * t86 - Ifges(5,6) * t87) / 0.2e1 + t83 * (mrSges(5,1) * t87 + mrSges(5,2) * t86) + t84 * (Ifges(5,4) * t86 - Ifges(5,2) * t87) / 0.2e1 + (-t33 * t86 - t34 * t87) * mrSges(5,3) + (Ifges(5,1) * t86 - Ifges(5,4) * t87) * t182 + m(6) * (t181 * t54 + t2 * t22 + t21 * t3 + t37 * t66 + t4 * t9 + t5 * t8) + (mrSges(6,2) * t37 - mrSges(6,3) * t3 + Ifges(6,1) * t19 + Ifges(6,4) * t20 + Ifges(6,5) * t119) * t58; m(3) * qJDD(2) + t57 * t14 + t58 * t15 + t29 * t40 + t30 * t41 + t92 * t43 + t93 * t44 + t86 * t71 - t87 * t72 + m(4) * (t126 * t70 + t128 * t69) + m(5) * (t12 * t93 + t13 * t92 - t33 * t87 + t34 * t86) + m(6) * (t2 * t58 + t29 * t9 + t3 * t57 + t30 * t8) + t197 * g(1); -t151 * t40 + t50 * t41 - t84 * t71 + t85 * t72 - t203 * qJD(1) ^ 2 * mrSges(4,3) + t149 + t152 + t153 + (g(2) * t116 + g(3) * t114) * t159 + (-t151 * t9 + t50 * t8 + t37) * m(6) + (t33 * t85 - t34 * t84 + t81) * m(5) + (-qJD(1) * t144 + t97) * m(4); (t201 - t146) * t175 - t83 * (mrSges(5,1) * t85 + mrSges(5,2) * t84) - qJD(4) * (Ifges(5,5) * t84 - Ifges(5,6) * t85) / 0.2e1 - t33 * t71 + t34 * t72 + Ifges(5,5) * t55 + Ifges(5,6) * t56 - t11 * t40 - t10 * t41 + t202 * g(1) - t85 * (Ifges(5,1) * t84 - t179) / 0.2e1 - t12 * mrSges(5,2) + t13 * mrSges(5,1) + (t131 * t15 + t134 * t14 - t85 * t25 + (-g(1) * t115 + t131 * t2 + t134 * t3 - t54 * t85 + (-g(2) * t114 + t175) * t113) * m(6) + (-t131 * t41 + t134 * t40 + (-t131 * t8 + t134 * t9) * m(6)) * qJD(5)) * pkin(4) - m(6) * (t10 * t8 + t11 * t9) + (t114 * t146 - t173) * g(2) - (-Ifges(5,2) * t85 + t46 + t80) * t84 / 0.2e1 + t50 * t199 + t205 + t45 * t182 + Ifges(5,3) * qJDD(4) + (t33 * t84 + t34 * t85) * mrSges(5,3); -g(1) * t150 - g(2) * t173 + t175 * t201 + t23 * t184 - t8 * t40 + t9 * t41 + t205;];
tau = t1;
