% Calculate joint inertia matrix for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2018-11-23 15:03
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:03:21
% EndTime: 2018-11-23 15:03:22
% DurationCPUTime: 0.86s
% Computational Cost: add. (1030->199), mult. (2257->298), div. (0->0), fcn. (2384->12), ass. (0->86)
t77 = sin(qJ(6));
t81 = cos(qJ(6));
t102 = t77 ^ 2 + t81 ^ 2;
t128 = mrSges(7,3) * t102;
t53 = -t81 * mrSges(7,1) + t77 * mrSges(7,2);
t127 = t53 - mrSges(6,1);
t107 = t77 * mrSges(7,3);
t78 = sin(qJ(5));
t79 = sin(qJ(4));
t82 = cos(qJ(5));
t83 = cos(qJ(4));
t49 = t78 * t79 - t82 * t83;
t51 = t78 * t83 + t82 * t79;
t24 = -t49 * mrSges(7,2) - t51 * t107;
t108 = t51 * t81;
t25 = t49 * mrSges(7,1) - mrSges(7,3) * t108;
t126 = t81 * t24 - t77 * t25;
t73 = sin(pkin(12));
t60 = t73 * pkin(2) + pkin(8);
t115 = pkin(9) + t60;
t41 = t115 * t83;
t99 = t115 * t79;
t19 = t78 * t41 + t82 * t99;
t94 = mrSges(7,1) * t77 + mrSges(7,2) * t81;
t22 = t94 * t51;
t125 = m(7) * t19 + t22;
t118 = pkin(5) * t49;
t75 = cos(pkin(12));
t61 = -t75 * pkin(2) - pkin(3);
t52 = -t83 * pkin(4) + t61;
t17 = -t51 * pkin(10) + t118 + t52;
t21 = t82 * t41 - t78 * t99;
t5 = t81 * t17 - t77 * t21;
t6 = t77 * t17 + t81 * t21;
t95 = -t5 * t77 + t6 * t81;
t124 = m(7) * t95 + t126;
t74 = sin(pkin(6));
t80 = sin(qJ(2));
t84 = cos(qJ(2));
t36 = (t73 * t84 + t75 * t80) * t74;
t76 = cos(pkin(6));
t26 = -t79 * t36 + t76 * t83;
t27 = t83 * t36 + t76 * t79;
t9 = -t82 * t26 + t78 * t27;
t123 = t9 ^ 2;
t122 = t19 ^ 2;
t34 = (t73 * t80 - t75 * t84) * t74;
t33 = t34 ^ 2;
t121 = t49 ^ 2;
t72 = t83 ^ 2;
t120 = 0.2e1 * t19;
t119 = m(6) * pkin(4);
t117 = t19 * t9;
t116 = t49 * t9;
t113 = mrSges(7,3) * t81;
t112 = Ifges(7,4) * t77;
t111 = Ifges(7,4) * t81;
t110 = t49 * t19;
t109 = t51 * t77;
t104 = Ifges(7,5) * t108 + Ifges(7,3) * t49;
t103 = Ifges(7,5) * t77 + Ifges(7,6) * t81;
t101 = t79 ^ 2 + t72;
t55 = Ifges(7,2) * t81 + t112;
t56 = Ifges(7,1) * t77 + t111;
t100 = t81 * t55 + t77 * t56 + Ifges(6,3);
t98 = t102 * pkin(10);
t63 = t78 * pkin(4) + pkin(10);
t97 = t102 * t63;
t28 = t49 * mrSges(6,1) + t51 * mrSges(6,2);
t11 = t78 * t26 + t82 * t27;
t2 = -t77 * t11 + t81 * t34;
t3 = t81 * t11 + t77 * t34;
t96 = -t2 * t77 + t3 * t81;
t54 = -t83 * mrSges(5,1) + t79 * mrSges(5,2);
t93 = -t26 * t79 + t27 * t83;
t92 = 0.2e1 * t128;
t91 = t128 * t51 + t49 * t53 - t28;
t90 = (t82 * mrSges(6,1) - t78 * mrSges(6,2)) * pkin(4);
t89 = -t11 * mrSges(6,2) - t2 * t107 + t3 * t113 + t127 * t9;
t14 = Ifges(7,6) * t49 + (-Ifges(7,2) * t77 + t111) * t51;
t15 = Ifges(7,5) * t49 + (Ifges(7,1) * t81 - t112) * t51;
t88 = -t21 * mrSges(6,2) - t5 * t107 + t6 * t113 + t77 * t15 / 0.2e1 + t81 * t14 / 0.2e1 - t55 * t109 / 0.2e1 + t56 * t108 / 0.2e1 + Ifges(6,5) * t51 + (t103 / 0.2e1 - Ifges(6,6)) * t49 + t127 * t19;
t68 = t76 ^ 2;
t64 = -t82 * pkin(4) - pkin(5);
t46 = t51 ^ 2;
t1 = [m(2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t123) + m(6) * (t11 ^ 2 + t123 + t33) + m(5) * (t26 ^ 2 + t27 ^ 2 + t33) + m(4) * (t36 ^ 2 + t33 + t68) + m(3) * (t68 + (t80 ^ 2 + t84 ^ 2) * t74 ^ 2); -t36 * mrSges(4,2) + t2 * t25 + t9 * t22 + t3 * t24 + (t84 * mrSges(3,1) - t80 * mrSges(3,2)) * t74 + (-t11 * t49 + t9 * t51) * mrSges(6,3) + t93 * mrSges(5,3) + (-mrSges(4,1) + t28 + t54) * t34 + m(7) * (t5 * t2 + t6 * t3 + t117) + m(6) * (t21 * t11 + t52 * t34 + t117) + m(5) * (t61 * t34 + t93 * t60) + m(4) * (-t34 * t75 + t36 * t73) * pkin(2); Ifges(5,2) * t72 + t22 * t120 + 0.2e1 * t6 * t24 + 0.2e1 * t5 * t25 + 0.2e1 * t52 * t28 + 0.2e1 * t61 * t54 + Ifges(3,3) + Ifges(4,3) + (Ifges(5,1) * t79 + 0.2e1 * Ifges(5,4) * t83) * t79 + (-0.2e1 * t21 * mrSges(6,3) + Ifges(6,2) * t49 + t104) * t49 + (mrSges(6,3) * t120 + Ifges(6,1) * t51 - t77 * t14 + t81 * t15 + (-Ifges(7,6) * t77 - (2 * Ifges(6,4))) * t49) * t51 + m(7) * (t5 ^ 2 + t6 ^ 2 + t122) + m(6) * (t21 ^ 2 + t52 ^ 2 + t122) + m(5) * (t101 * t60 ^ 2 + t61 ^ 2) + m(4) * (t73 ^ 2 + t75 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (t75 * mrSges(4,1) - t73 * mrSges(4,2)) * pkin(2) + 0.2e1 * t101 * t60 * mrSges(5,3); m(4) * t76 + m(7) * (t96 * t51 + t116) + m(6) * (t51 * t11 + t116) + m(5) * (t83 * t26 + t79 * t27); t49 * t22 + t126 * t51 + m(7) * (t95 * t51 + t110) + m(6) * (t51 * t21 + t110); m(4) + m(5) * t101 + m(6) * (t46 + t121) + m(7) * (t102 * t46 + t121); t26 * mrSges(5,1) - t27 * mrSges(5,2) + m(7) * (t96 * t63 + t64 * t9) + (t11 * t78 - t82 * t9) * t119 + t89; (m(6) * (-t19 * t82 + t21 * t78) + (-t78 * t49 - t82 * t51) * mrSges(6,3)) * pkin(4) + t88 + (-t79 * mrSges(5,1) - t83 * mrSges(5,2)) * t60 + Ifges(5,5) * t79 + Ifges(5,6) * t83 + t125 * t64 + t124 * t63; m(7) * (t64 * t49 + t51 * t97) + (-t49 * t82 + t51 * t78) * t119 + t91 - t54; 0.2e1 * t64 * t53 + Ifges(5,3) + 0.2e1 * t90 + t63 * t92 + m(7) * (t102 * t63 ^ 2 + t64 ^ 2) + m(6) * (t78 ^ 2 + t82 ^ 2) * pkin(4) ^ 2 + t100; m(7) * (-pkin(5) * t9 + t96 * pkin(10)) + t89; -t125 * pkin(5) + t124 * pkin(10) + t88; m(7) * (t51 * t98 - t118) + t91; m(7) * (-pkin(5) * t64 + pkin(10) * t97) + (-pkin(5) + t64) * t53 + t90 + (t97 + t98) * mrSges(7,3) + t100; -0.2e1 * pkin(5) * t53 + m(7) * (t102 * pkin(10) ^ 2 + pkin(5) ^ 2) + pkin(10) * t92 + t100; t2 * mrSges(7,1) - t3 * mrSges(7,2); t5 * mrSges(7,1) - t6 * mrSges(7,2) - Ifges(7,6) * t109 + t104; -t22; -t94 * t63 + t103; -t94 * pkin(10) + t103; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
