% Calculate time derivative of joint inertia matrix for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 14:53
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 14:51:12
% EndTime: 2018-11-16 14:51:14
% DurationCPUTime: 1.24s
% Computational Cost: add. (2232->193), mult. (5258->319), div. (0->0), fcn. (5187->8), ass. (0->103)
t77 = sin(qJ(5));
t75 = t77 ^ 2;
t81 = cos(qJ(5));
t76 = t81 ^ 2;
t113 = t75 + t76;
t108 = qJD(5) * t81;
t137 = qJD(2) + qJD(3);
t79 = sin(qJ(3));
t80 = sin(qJ(2));
t83 = cos(qJ(3));
t84 = cos(qJ(2));
t57 = t79 * t80 - t83 * t84;
t44 = t137 * t57;
t93 = t79 * t84 + t83 * t80;
t45 = t137 * t93;
t78 = sin(qJ(4));
t82 = cos(qJ(4));
t94 = t82 * t57 + t78 * t93;
t20 = qJD(4) * t94 + t44 * t82 + t45 * t78;
t40 = t57 * t78 - t82 * t93;
t92 = t40 * t108 + t20 * t77;
t138 = t113 * t82;
t99 = mrSges(6,1) * t77 + mrSges(6,2) * t81;
t59 = t99 * qJD(5);
t70 = -pkin(3) * t82 - pkin(4);
t47 = t70 * t59;
t111 = qJD(4) * t78;
t129 = mrSges(6,1) * t81;
t62 = mrSges(6,2) * t77 - t129;
t55 = pkin(3) * t62 * t111;
t110 = qJD(4) * t82;
t106 = pkin(3) * t110;
t101 = mrSges(6,3) * t106;
t65 = t75 * t101;
t66 = t76 * t101;
t136 = t47 + t55 + t65 + t66;
t21 = qJD(4) * t40 + t44 * t78 - t82 * t45;
t34 = -qJD(2) * t80 * pkin(2) - pkin(3) * t45;
t8 = pkin(4) * t21 - pkin(6) * t20 + t34;
t134 = 0.2e1 * t8;
t133 = 0.2e1 * t34;
t130 = pkin(4) * t59;
t72 = t84 * pkin(2) + pkin(1);
t128 = Ifges(6,4) * t77;
t127 = Ifges(6,4) * t81;
t126 = Ifges(6,6) * t77;
t117 = t78 * t79;
t71 = pkin(2) * t83 + pkin(3);
t36 = t71 * t110 + (-t79 * t111 + (t82 * t83 - t117) * qJD(3)) * pkin(2);
t124 = t36 * mrSges(5,2);
t123 = t36 * t75;
t122 = t36 * t76;
t121 = t40 * t77;
t120 = t40 * t81;
t116 = t79 * t82;
t115 = t81 * t20;
t114 = Ifges(6,5) * t115 + Ifges(6,3) * t21;
t51 = pkin(2) * t116 + t78 * t71;
t112 = pkin(3) * qJD(4);
t109 = qJD(5) * t77;
t107 = 0.2e1 * t84;
t104 = t113 * t36;
t102 = -(2 * Ifges(5,4)) - t126;
t46 = -pkin(3) * t57 + t72;
t100 = -t78 * mrSges(5,1) - t82 * mrSges(5,2);
t98 = Ifges(6,1) * t81 - t128;
t97 = -Ifges(6,2) * t77 + t127;
t96 = Ifges(6,5) * t77 + Ifges(6,6) * t81;
t24 = mrSges(6,2) * t94 - mrSges(6,3) * t121;
t25 = -mrSges(6,1) * t94 - mrSges(6,3) * t120;
t95 = t24 * t81 - t25 * t77;
t50 = -pkin(2) * t117 + t71 * t82;
t91 = t109 * t40 - t115;
t60 = t97 * qJD(5);
t61 = t98 * qJD(5);
t63 = Ifges(6,2) * t81 + t128;
t64 = Ifges(6,1) * t77 + t127;
t90 = t64 * t108 - t109 * t63 + t81 * t60 + t77 * t61;
t89 = (-mrSges(4,1) * t79 - mrSges(4,2) * t83) * qJD(3) * pkin(2);
t6 = mrSges(6,1) * t21 + mrSges(6,3) * t91;
t7 = -mrSges(6,2) * t21 - mrSges(6,3) * t92;
t88 = -t77 * t6 + t81 * t7 + (-t24 * t77 - t25 * t81) * qJD(5);
t14 = -Ifges(6,6) * t94 + t40 * t97;
t15 = -Ifges(6,5) * t94 + t40 * t98;
t3 = -Ifges(6,4) * t91 - Ifges(6,2) * t92 + t21 * Ifges(6,6);
t4 = -Ifges(6,1) * t91 - Ifges(6,4) * t92 + t21 * Ifges(6,5);
t73 = Ifges(6,5) * t108;
t87 = t77 * t4 / 0.2e1 + t64 * t115 / 0.2e1 + t15 * t108 / 0.2e1 + Ifges(5,5) * t20 + t81 * t3 / 0.2e1 - t60 * t121 / 0.2e1 + t61 * t120 / 0.2e1 - t94 * (-Ifges(6,6) * t109 + t73) / 0.2e1 + (t96 / 0.2e1 - Ifges(5,6)) * t21 - t92 * t63 / 0.2e1 - (t40 * t64 + t14) * t109 / 0.2e1;
t37 = t71 * t111 + (t79 * t110 + (t78 * t83 + t116) * qJD(3)) * pkin(2);
t31 = t37 * t62;
t32 = mrSges(6,3) * t123;
t33 = mrSges(6,3) * t122;
t35 = t37 * mrSges(5,1);
t48 = -pkin(4) - t50;
t43 = t48 * t59;
t86 = t31 + t32 + t33 - t35 + t43 + t90 - t124;
t85 = Ifges(4,5) * t44 + Ifges(4,6) * t45 + t87;
t69 = pkin(3) * t78 + pkin(6);
t49 = pkin(6) + t51;
t23 = t99 * t40;
t22 = -pkin(4) * t94 - pkin(6) * t40 + t46;
t5 = mrSges(6,1) * t92 - mrSges(6,2) * t91;
t1 = [0.2e1 * t72 * (-mrSges(4,1) * t45 + mrSges(4,2) * t44) + 0.2e1 * t57 * Ifges(4,2) * t45 - 0.2e1 * t44 * t93 * Ifges(4,1) + (t25 * t134 + t20 * t15) * t81 + (t24 * t134 - t20 * t14) * t77 + 0.2e1 * (t57 * t44 - t45 * t93) * Ifges(4,4) + (m(6) * t113 * t134 + 0.2e1 * qJD(5) * t95 + 0.2e1 * t81 * t6 + 0.2e1 * t77 * t7) * t22 - (mrSges(5,1) * t133 + ((2 * Ifges(5,2)) + Ifges(6,3)) * t21 + t102 * t20 + t114) * t94 + (mrSges(5,2) * t133 + 0.2e1 * Ifges(5,1) * t20 - t77 * t3 + t81 * t4 + (Ifges(6,5) * t81 + t102) * t21 + (-t81 * t14 - t77 * t15 + t94 * t96) * qJD(5)) * t40 + ((-pkin(1) * mrSges(3,2) + Ifges(3,4) * t84) * t107 + (-0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t80 + (Ifges(3,1) - Ifges(3,2)) * t107 - 0.2e1 * (m(4) * t72 - mrSges(4,1) * t57 - mrSges(4,2) * t93) * pkin(2)) * t80) * qJD(2) + 0.2e1 * (m(5) * t34 + mrSges(5,1) * t21 + mrSges(5,2) * t20) * t46; t48 * t5 + t37 * t23 + t88 * t49 + t95 * t36 + (-Ifges(3,5) * t84 + Ifges(3,6) * t80) * qJD(2) + (-t50 * t20 - t51 * t21 + t36 * t94 + t37 * t40) * mrSges(5,3) + t85 + (-t44 * t83 + t45 * t79 + (t57 * t83 - t79 * t93) * qJD(3)) * pkin(2) * mrSges(4,3); -0.2e1 * t124 + 0.2e1 * t31 + 0.2e1 * t32 + 0.2e1 * t33 - 0.2e1 * t35 + 0.2e1 * t43 + 0.2e1 * t89 + 0.2e1 * m(6) * (t104 * t49 + t37 * t48) + 0.2e1 * m(5) * (t36 * t51 - t37 * t50) + t90; t70 * t5 + t85 + t88 * t69 + ((-t20 * t82 - t21 * t78) * mrSges(5,3) + ((mrSges(5,3) * t40 + t23) * t78 + (mrSges(5,3) * t94 + t95) * t82) * qJD(4)) * pkin(3); m(6) * (t37 * t70 + (t122 + t123) * t69) + (m(5) * (t36 * t78 - t37 * t82) + (m(6) * (t138 * t49 + t48 * t78) + m(5) * (-t50 * t78 + t51 * t82) + t100) * qJD(4)) * pkin(3) + t89 + t86 + t136; 0.2e1 * t47 + 0.2e1 * t55 + 0.2e1 * t65 + 0.2e1 * t66 + 0.2e1 * (m(6) * (t138 * t69 + t70 * t78) + t100) * t112 + t90; -pkin(4) * t5 + pkin(6) * t88 + t87; -t130 + m(6) * (-pkin(4) * t37 + pkin(6) * t104) + t86; -t130 + (m(6) * (-pkin(4) * t78 + t138 * pkin(6)) + t100) * t112 + t90 + t136; t90 - 0.2e1 * t130; t8 * t129 + (-mrSges(6,2) * t8 - Ifges(6,6) * t20) * t77 + (-t22 * t99 - t40 * t96) * qJD(5) + t114; t73 - t99 * t36 + (t49 * t62 - t126) * qJD(5); t73 - t99 * t106 + (t62 * t69 - t126) * qJD(5); t73 + (pkin(6) * t62 - t126) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11); t1(2) t1(3) t1(5) t1(8) t1(12); t1(4) t1(5) t1(6) t1(9) t1(13); t1(7) t1(8) t1(9) t1(10) t1(14); t1(11) t1(12) t1(13) t1(14) t1(15);];
Mq  = res;
