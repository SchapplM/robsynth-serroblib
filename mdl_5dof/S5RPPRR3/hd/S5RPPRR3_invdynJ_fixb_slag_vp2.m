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
% m [6x1]
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
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:14:09
% EndTime: 2022-01-23 09:14:15
% DurationCPUTime: 4.35s
% Computational Cost: add. (3080->312), mult. (6900->421), div. (0->0), fcn. (4932->16), ass. (0->137)
t126 = sin(pkin(8));
t100 = pkin(1) * t126 + qJ(3);
t90 = qJD(1) * qJD(3) + qJDD(1) * t100;
t118 = qJDD(4) + qJDD(5);
t123 = qJD(4) + qJD(5);
t130 = sin(qJ(5));
t133 = cos(qJ(5));
t125 = sin(pkin(9));
t131 = sin(qJ(4));
t127 = cos(pkin(9));
t134 = cos(qJ(4));
t160 = t127 * t134;
t92 = -t125 * t131 + t160;
t84 = t92 * qJD(1);
t93 = t125 * t134 + t127 * t131;
t85 = t93 * qJD(1);
t149 = -t130 * t85 + t133 * t84;
t50 = t130 * t84 + t133 * t85;
t174 = Ifges(6,4) * t50;
t86 = t92 * qJD(4);
t55 = qJD(1) * t86 + qJDD(1) * t93;
t87 = t93 * qJD(4);
t56 = -qJD(1) * t87 + qJDD(1) * t92;
t19 = qJD(5) * t149 + t130 * t56 + t133 * t55;
t110 = t127 * qJD(2);
t169 = pkin(6) * qJD(1);
t97 = t100 * qJD(1);
t67 = t110 + (-t97 - t169) * t125;
t75 = t125 * qJD(2) + t127 * t97;
t68 = t127 * t169 + t75;
t34 = t131 * t67 + t134 * t68;
t108 = t127 * qJDD(2);
t63 = t108 + (-pkin(6) * qJDD(1) - t90) * t125;
t155 = qJDD(1) * t127;
t70 = t125 * qJDD(2) + t127 * t90;
t64 = pkin(6) * t155 + t70;
t13 = -qJD(4) * t34 - t131 * t64 + t134 * t63;
t6 = qJDD(4) * pkin(4) - pkin(7) * t55 + t13;
t159 = qJD(4) * t134;
t163 = t131 * t68;
t12 = -qJD(4) * t163 + t131 * t63 + t134 * t64 + t67 * t159;
t7 = pkin(7) * t56 + t12;
t28 = pkin(7) * t84 + t34;
t165 = t130 * t28;
t33 = t134 * t67 - t163;
t27 = -pkin(7) * t85 + t33;
t26 = qJD(4) * pkin(4) + t27;
t8 = t133 * t26 - t165;
t2 = qJD(5) * t8 + t130 * t6 + t133 * t7;
t20 = -qJD(5) * t50 - t130 * t55 + t133 * t56;
t42 = Ifges(6,4) * t149;
t24 = Ifges(6,1) * t50 + Ifges(6,5) * t123 + t42;
t162 = t133 * t28;
t9 = t130 * t26 + t162;
t3 = -qJD(5) * t9 - t130 * t7 + t133 * t6;
t128 = cos(pkin(8));
t105 = -pkin(1) * t128 - pkin(2);
t116 = t127 * pkin(3);
t190 = t105 - t116;
t83 = qJD(1) * t190 + qJD(3);
t54 = -pkin(4) * t84 + t83;
t201 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t19 + Ifges(6,6) * t20 + Ifges(6,3) * t118 - (Ifges(6,5) * t149 - Ifges(6,6) * t50) * t123 / 0.2e1 + (t149 * t8 + t50 * t9) * mrSges(6,3) - (-Ifges(6,2) * t50 + t24 + t42) * t149 / 0.2e1 - t54 * (mrSges(6,1) * t50 + mrSges(6,2) * t149) - (Ifges(6,1) * t149 - t174) * t50 / 0.2e1;
t200 = t125 ^ 2 + t127 ^ 2;
t122 = pkin(9) + qJ(4);
t111 = sin(t122);
t113 = cos(t122);
t115 = qJ(5) + t122;
t102 = sin(t115);
t103 = cos(t115);
t150 = t103 * mrSges(6,1) - t102 * mrSges(6,2);
t199 = -t113 * mrSges(5,1) + t111 * mrSges(5,2) - t150;
t124 = qJ(1) + pkin(8);
t112 = sin(t124);
t114 = cos(t124);
t187 = g(1) * t114 + g(2) * t112;
t158 = m(4) + m(5) + m(6);
t194 = -m(3) - t158;
t198 = pkin(1) * t194 - mrSges(2,1);
t23 = Ifges(6,2) * t149 + Ifges(6,6) * t123 + t174;
t196 = t23 / 0.2e1;
t69 = -t125 * t90 + t108;
t188 = -t125 * t69 + t127 * t70;
t104 = t116 + pkin(2);
t147 = -mrSges(4,1) * t127 + mrSges(4,2) * t125;
t184 = -m(5) * t104 - mrSges(3,1) - m(6) * (pkin(4) * t113 + t104) - m(4) * pkin(2) + t147 + t199;
t129 = -pkin(6) - qJ(3);
t183 = -m(4) * qJ(3) + m(5) * t129 - m(6) * (pkin(7) - t129) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t180 = t50 / 0.2e1;
t178 = t85 / 0.2e1;
t177 = pkin(4) * t87;
t175 = Ifges(5,4) * t85;
t170 = pkin(6) + t100;
t88 = t170 * t125;
t89 = t170 * t127;
t52 = -t131 * t88 + t134 * t89;
t156 = qJDD(1) * t125;
t152 = -t56 * mrSges(5,1) + t55 * mrSges(5,2);
t151 = -t20 * mrSges(6,1) + t19 * mrSges(6,2);
t51 = -t131 * t89 - t134 * t88;
t148 = -mrSges(4,1) * t155 + mrSges(4,2) * t156;
t145 = mrSges(6,1) * t102 + mrSges(6,2) * t103;
t144 = -t125 * (-t125 * t97 + t110) + t127 * t75;
t38 = -pkin(7) * t93 + t51;
t39 = pkin(7) * t92 + t52;
t21 = -t130 * t39 + t133 * t38;
t22 = t130 * t38 + t133 * t39;
t57 = -t130 * t93 + t133 * t92;
t58 = t130 * t92 + t133 * t93;
t81 = qJDD(1) * t190 + qJDD(3);
t35 = -t88 * t159 + qJD(3) * t160 + (-qJD(3) * t125 - qJD(4) * t89) * t131;
t36 = -qJD(3) * t93 - qJD(4) * t52;
t135 = cos(qJ(1));
t132 = sin(qJ(1));
t95 = qJDD(1) * t105 + qJDD(3);
t80 = Ifges(5,4) * t84;
t72 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t85;
t71 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t84;
t66 = -pkin(4) * t92 + t190;
t46 = t85 * Ifges(5,1) + Ifges(5,5) * qJD(4) + t80;
t45 = t84 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t175;
t44 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t56;
t43 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t55;
t41 = mrSges(6,1) * t123 - mrSges(6,3) * t50;
t40 = -mrSges(6,2) * t123 + mrSges(6,3) * t149;
t37 = -pkin(4) * t56 + t81;
t32 = -pkin(7) * t86 + t36;
t31 = -pkin(7) * t87 + t35;
t30 = -qJD(5) * t58 - t130 * t86 - t133 * t87;
t29 = qJD(5) * t57 - t130 * t87 + t133 * t86;
t25 = -mrSges(6,1) * t149 + mrSges(6,2) * t50;
t15 = -mrSges(6,2) * t118 + mrSges(6,3) * t20;
t14 = mrSges(6,1) * t118 - mrSges(6,3) * t19;
t11 = t133 * t27 - t165;
t10 = -t130 * t27 - t162;
t5 = -qJD(5) * t22 - t130 * t31 + t133 * t32;
t4 = qJD(5) * t21 + t130 * t32 + t133 * t31;
t1 = [(mrSges(2,2) * t132 + t112 * t183 + t114 * t184 + t135 * t198) * g(2) + (mrSges(2,2) * t135 - t112 * t184 + t114 * t183 - t132 * t198) * g(1) + t25 * t177 + (-mrSges(6,1) * t37 + mrSges(6,3) * t2 + Ifges(6,4) * t19 + Ifges(6,2) * t20 + Ifges(6,6) * t118) * t57 + t123 * (Ifges(6,5) * t29 + Ifges(6,6) * t30) / 0.2e1 + t86 * t46 / 0.2e1 - t87 * t45 / 0.2e1 + t35 * t71 + t36 * t72 + t52 * t44 + t54 * (-mrSges(6,1) * t30 + mrSges(6,2) * t29) + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * t128 * mrSges(3,1) - 0.2e1 * t126 * mrSges(3,2) + m(3) * (t126 ^ 2 + t128 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + m(5) * (t12 * t52 + t13 * t51 + t190 * t81 + t33 * t36 + t34 * t35) + t190 * t152 + t51 * t43 + t4 * t40 + t5 * t41 + t29 * t24 / 0.2e1 + (t200 * t90 + t188) * mrSges(4,3) + (-mrSges(5,1) * t81 + mrSges(5,3) * t12 + Ifges(5,4) * t55 + Ifges(5,2) * t56 + Ifges(5,6) * qJDD(4)) * t92 + t30 * t196 + (Ifges(4,4) * t125 + Ifges(4,2) * t127) * t155 + (Ifges(4,1) * t125 + Ifges(4,4) * t127) * t156 + (mrSges(5,2) * t81 - mrSges(5,3) * t13 + Ifges(5,1) * t55 + Ifges(5,4) * t56 + Ifges(5,5) * qJDD(4)) * t93 + t21 * t14 + t22 * t15 + t95 * t147 + t105 * t148 + t66 * t151 + t84 * (Ifges(5,4) * t86 - Ifges(5,2) * t87) / 0.2e1 + t83 * (mrSges(5,1) * t87 + mrSges(5,2) * t86) + qJD(4) * (Ifges(5,5) * t86 - Ifges(5,6) * t87) / 0.2e1 + (-t33 * t86 - t34 * t87) * mrSges(5,3) + (Ifges(5,1) * t86 - Ifges(5,4) * t87) * t178 + m(6) * (t177 * t54 + t2 * t22 + t21 * t3 + t37 * t66 + t4 * t9 + t5 * t8) + (-t8 * t29 + t9 * t30) * mrSges(6,3) + (mrSges(6,2) * t37 - mrSges(6,3) * t3 + Ifges(6,1) * t19 + Ifges(6,4) * t20 + Ifges(6,5) * t118) * t58 + t149 * (Ifges(6,4) * t29 + Ifges(6,2) * t30) / 0.2e1 + (Ifges(6,1) * t29 + Ifges(6,4) * t30) * t180 + m(4) * (t144 * qJD(3) + t100 * t188 + t105 * t95); m(3) * qJDD(2) + t57 * t14 + t58 * t15 + t29 * t40 + t30 * t41 + t92 * t43 + t93 * t44 + t86 * t71 - t87 * t72 + m(4) * (t125 * t70 + t127 * t69) + m(5) * (t12 * t93 + t13 * t92 - t33 * t87 + t34 * t86) + m(6) * (t2 * t58 + t29 * t9 + t3 * t57 + t30 * t8) + t194 * g(3); -t149 * t40 + t50 * t41 - t84 * t71 + t85 * t72 - t200 * qJD(1) ^ 2 * mrSges(4,3) + t148 + t151 + t152 + (-g(1) * t112 + g(2) * t114) * t158 + (-t149 * t9 + t50 * t8 + t37) * m(6) + (t33 * t85 - t34 * t84 + t81) * m(5) + (-qJD(1) * t144 + t95) * m(4); -m(6) * (t10 * t8 + t11 * t9) + t187 * (mrSges(5,1) * t111 + mrSges(5,2) * t113 + t145) - t83 * (mrSges(5,1) * t85 + mrSges(5,2) * t84) - qJD(4) * (Ifges(5,5) * t84 - Ifges(5,6) * t85) / 0.2e1 - t33 * t71 + t34 * t72 + Ifges(5,5) * t55 + Ifges(5,6) * t56 - (-Ifges(5,2) * t85 + t46 + t80) * t84 / 0.2e1 + t50 * t196 - t11 * t40 - t10 * t41 - t12 * mrSges(5,2) + t13 * mrSges(5,1) + t201 + t199 * g(3) + (t33 * t84 + t34 * t85) * mrSges(5,3) - t85 * (Ifges(5,1) * t84 - t175) / 0.2e1 + t45 * t178 + (t130 * t15 + t133 * t14 - t85 * t25 + (-g(3) * t113 + t187 * t111 + t130 * t2 + t133 * t3 - t54 * t85) * m(6) + (-t130 * t41 + t133 * t40 + (-t130 * t8 + t133 * t9) * m(6)) * qJD(5)) * pkin(4) + Ifges(5,3) * qJDD(4); -g(3) * t150 + t187 * t145 + t23 * t180 - t8 * t40 + t9 * t41 + t201;];
tau = t1;
