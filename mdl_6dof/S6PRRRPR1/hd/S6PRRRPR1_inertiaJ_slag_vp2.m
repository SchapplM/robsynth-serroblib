% Calculate joint inertia matrix for
% S6PRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:22
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:22:20
% EndTime: 2018-11-23 15:22:21
% DurationCPUTime: 0.88s
% Computational Cost: add. (1716->233), mult. (3538->346), div. (0->0), fcn. (4029->12), ass. (0->95)
t86 = sin(qJ(6));
t90 = cos(qJ(6));
t63 = -mrSges(7,1) * t90 + mrSges(7,2) * t86;
t137 = t63 - mrSges(6,1);
t136 = m(6) * pkin(4);
t117 = t86 * mrSges(7,3);
t87 = sin(qJ(4));
t88 = sin(qJ(3));
t91 = cos(qJ(4));
t92 = cos(qJ(3));
t59 = -t87 * t88 + t91 * t92;
t60 = t87 * t92 + t88 * t91;
t82 = sin(pkin(12));
t84 = cos(pkin(12));
t40 = -t84 * t59 + t60 * t82;
t41 = t59 * t82 + t60 * t84;
t23 = -mrSges(7,2) * t40 - t117 * t41;
t122 = t41 * t90;
t24 = mrSges(7,1) * t40 - mrSges(7,3) * t122;
t135 = t90 * t23 - t86 * t24;
t130 = -pkin(9) - pkin(8);
t109 = t130 * t88;
t110 = t130 * t92;
t45 = t87 * t109 - t91 * t110;
t32 = qJ(5) * t59 + t45;
t44 = t109 * t91 + t110 * t87;
t98 = -t60 * qJ(5) + t44;
t15 = t32 * t82 - t84 * t98;
t134 = t15 ^ 2;
t83 = sin(pkin(6));
t89 = sin(qJ(2));
t119 = t83 * t89;
t85 = cos(pkin(6));
t53 = -t119 * t88 + t85 * t92;
t54 = t119 * t92 + t85 * t88;
t33 = t53 * t91 - t54 * t87;
t34 = t53 * t87 + t54 * t91;
t19 = -t84 * t33 + t34 * t82;
t133 = t19 ^ 2;
t81 = t92 ^ 2;
t132 = 0.2e1 * t15;
t131 = 0.2e1 * t63;
t129 = pkin(3) * t87;
t74 = -pkin(3) * t92 - pkin(2);
t46 = -pkin(4) * t59 + t74;
t14 = pkin(5) * t40 - pkin(10) * t41 + t46;
t17 = t84 * t32 + t82 * t98;
t3 = t14 * t86 + t17 * t90;
t128 = t3 * t90;
t127 = Ifges(7,4) * t86;
t126 = Ifges(7,4) * t90;
t93 = cos(qJ(2));
t118 = t83 * t93;
t21 = t33 * t82 + t34 * t84;
t10 = -t118 * t86 + t21 * t90;
t125 = t10 * t90;
t124 = t15 * t19;
t123 = t41 * t86;
t73 = pkin(3) * t91 + pkin(4);
t50 = -t129 * t82 + t73 * t84;
t121 = t50 * mrSges(6,1);
t51 = t84 * t129 + t82 * t73;
t120 = t51 * mrSges(6,2);
t114 = Ifges(7,5) * t122 + Ifges(7,3) * t40;
t113 = Ifges(7,5) * t86 + Ifges(7,6) * t90;
t112 = t86 ^ 2 + t90 ^ 2;
t111 = t88 ^ 2 + t81;
t71 = pkin(4) * t82 + pkin(10);
t108 = t112 * t71;
t25 = t40 * mrSges(6,1) + t41 * mrSges(6,2);
t65 = Ifges(7,2) * t90 + t127;
t66 = Ifges(7,1) * t86 + t126;
t107 = t90 * t65 + t86 * t66 + Ifges(5,3) + Ifges(6,3);
t2 = t14 * t90 - t17 * t86;
t106 = -t2 * t86 + t128;
t9 = -t118 * t90 - t21 * t86;
t105 = -t86 * t9 + t125;
t104 = t84 * mrSges(6,1) - t82 * mrSges(6,2);
t103 = mrSges(7,1) * t86 + mrSges(7,2) * t90;
t102 = -t53 * t88 + t54 * t92;
t101 = 0.2e1 * mrSges(7,3) * t112;
t100 = (mrSges(5,1) * t91 - mrSges(5,2) * t87) * pkin(3);
t99 = t33 * mrSges(5,1) - t34 * mrSges(5,2) - t21 * mrSges(6,2) + mrSges(7,3) * t125 - t117 * t9 + t137 * t19;
t6 = Ifges(7,6) * t40 + (-Ifges(7,2) * t86 + t126) * t41;
t7 = Ifges(7,5) * t40 + (Ifges(7,1) * t90 - t127) * t41;
t97 = -t45 * mrSges(5,2) - t17 * mrSges(6,2) + mrSges(7,3) * t128 - t117 * t2 - t65 * t123 / 0.2e1 + t66 * t122 / 0.2e1 + Ifges(6,5) * t41 + t86 * t7 / 0.2e1 + t44 * mrSges(5,1) + t90 * t6 / 0.2e1 + Ifges(5,6) * t59 + Ifges(5,5) * t60 + (t113 / 0.2e1 - Ifges(6,6)) * t40 + t137 * t15;
t77 = t83 ^ 2;
t72 = -pkin(4) * t84 - pkin(5);
t70 = t77 * t93 ^ 2;
t64 = -mrSges(4,1) * t92 + mrSges(4,2) * t88;
t49 = pkin(10) + t51;
t48 = -pkin(5) - t50;
t42 = -mrSges(5,1) * t59 + mrSges(5,2) * t60;
t22 = t103 * t41;
t1 = [m(2) + m(7) * (t10 ^ 2 + t9 ^ 2 + t133) + m(6) * (t21 ^ 2 + t133 + t70) + m(5) * (t33 ^ 2 + t34 ^ 2 + t70) + m(4) * (t53 ^ 2 + t54 ^ 2 + t70) + m(3) * (t77 * t89 ^ 2 + t85 ^ 2 + t70); t10 * t23 + t19 * t22 + t9 * t24 + (t19 * t41 - t21 * t40) * mrSges(6,3) + (-t33 * t60 + t34 * t59) * mrSges(5,3) + t102 * mrSges(4,3) + (-t89 * mrSges(3,2) + (mrSges(3,1) - t25 - t42 - t64) * t93) * t83 + m(7) * (t10 * t3 + t2 * t9 + t124) + m(6) * (-t118 * t46 + t17 * t21 + t124) + m(5) * (-t118 * t74 + t33 * t44 + t34 * t45) + m(4) * (pkin(2) * t118 + pkin(8) * t102); Ifges(4,2) * t81 - 0.2e1 * pkin(2) * t64 + t22 * t132 + 0.2e1 * t2 * t24 + 0.2e1 * t3 * t23 + 0.2e1 * t46 * t25 + 0.2e1 * t74 * t42 + Ifges(3,3) + (Ifges(4,1) * t88 + 0.2e1 * Ifges(4,4) * t92) * t88 + (-0.2e1 * t44 * mrSges(5,3) + Ifges(5,1) * t60) * t60 + 0.2e1 * t111 * pkin(8) * mrSges(4,3) + (0.2e1 * t45 * mrSges(5,3) + 0.2e1 * Ifges(5,4) * t60 + Ifges(5,2) * t59) * t59 + (-0.2e1 * t17 * mrSges(6,3) + Ifges(6,2) * t40 + t114) * t40 + (mrSges(6,3) * t132 + Ifges(6,1) * t41 - t86 * t6 + t90 * t7 + (-Ifges(7,6) * t86 - (2 * Ifges(6,4))) * t40) * t41 + m(4) * (pkin(8) ^ 2 * t111 + pkin(2) ^ 2) + m(5) * (t44 ^ 2 + t45 ^ 2 + t74 ^ 2) + m(6) * (t17 ^ 2 + t46 ^ 2 + t134) + m(7) * (t2 ^ 2 + t3 ^ 2 + t134); t53 * mrSges(4,1) - t54 * mrSges(4,2) + m(7) * (t105 * t49 + t19 * t48) + m(6) * (-t19 * t50 + t21 * t51) + m(5) * (t33 * t91 + t34 * t87) * pkin(3) + t99; t97 + (m(5) * (t44 * t91 + t45 * t87) + (t87 * t59 - t91 * t60) * mrSges(5,3)) * pkin(3) + m(7) * (t106 * t49 + t15 * t48) + t135 * t49 + Ifges(4,5) * t88 + Ifges(4,6) * t92 + t48 * t22 + m(6) * (-t15 * t50 + t17 * t51) + (-t51 * t40 - t50 * t41) * mrSges(6,3) + (-t88 * mrSges(4,1) - t92 * mrSges(4,2)) * pkin(8); 0.2e1 * t121 - 0.2e1 * t120 + t48 * t131 + Ifges(4,3) + 0.2e1 * t100 + t49 * t101 + m(7) * (t112 * t49 ^ 2 + t48 ^ 2) + m(6) * (t50 ^ 2 + t51 ^ 2) + m(5) * (t87 ^ 2 + t91 ^ 2) * pkin(3) ^ 2 + t107; m(7) * (t105 * t71 + t72 * t19) + (-t19 * t84 + t21 * t82) * t136 + t99; t97 + (m(6) * (-t15 * t84 + t17 * t82) + (-t40 * t82 - t41 * t84) * mrSges(6,3)) * pkin(4) + (m(7) * t15 + t22) * t72 + (m(7) * t106 + t135) * t71; m(7) * (t108 * t49 + t48 * t72) - t120 + t121 + (t48 + t72) * t63 + t100 + (m(6) * (t50 * t84 + t51 * t82) + t104) * pkin(4) + (t112 * t49 + t108) * mrSges(7,3) + t107; t72 * t131 + t71 * t101 + m(7) * (t112 * t71 ^ 2 + t72 ^ 2) + t107 + (0.2e1 * t104 + (t82 ^ 2 + t84 ^ 2) * t136) * pkin(4); m(7) * (t10 * t86 + t9 * t90) - m(6) * t118; t86 * t23 + t90 * t24 + m(7) * (t2 * t90 + t3 * t86) + m(6) * t46 + t25; 0; 0; m(7) * t112 + m(6); mrSges(7,1) * t9 - mrSges(7,2) * t10; mrSges(7,1) * t2 - mrSges(7,2) * t3 - Ifges(7,6) * t123 + t114; -t103 * t49 + t113; -t103 * t71 + t113; -t63; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
