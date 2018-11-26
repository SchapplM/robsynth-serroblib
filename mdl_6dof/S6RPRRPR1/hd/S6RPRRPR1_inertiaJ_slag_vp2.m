% Calculate joint inertia matrix for
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2018-11-23 16:16
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:15:43
% EndTime: 2018-11-23 16:15:44
% DurationCPUTime: 0.82s
% Computational Cost: add. (1665->201), mult. (3024->298), div. (0->0), fcn. (3277->10), ass. (0->83)
t77 = sin(qJ(6));
t80 = cos(qJ(6));
t103 = t77 ^ 2 + t80 ^ 2;
t128 = mrSges(7,3) * t103;
t82 = cos(qJ(3));
t127 = t82 ^ 2;
t126 = m(6) * pkin(4);
t109 = t77 * mrSges(7,3);
t78 = sin(qJ(4));
t79 = sin(qJ(3));
t81 = cos(qJ(4));
t52 = -t78 * t79 + t81 * t82;
t53 = t78 * t82 + t79 * t81;
t73 = sin(pkin(11));
t75 = cos(pkin(11));
t35 = -t75 * t52 + t53 * t73;
t37 = t52 * t73 + t53 * t75;
t15 = -mrSges(7,2) * t35 - t37 * t109;
t112 = t37 * t80;
t16 = mrSges(7,1) * t35 - mrSges(7,3) * t112;
t91 = t80 * t15 - t77 * t16;
t74 = sin(pkin(10));
t63 = pkin(1) * t74 + pkin(7);
t118 = pkin(8) + t63;
t100 = t118 * t82;
t101 = t118 * t79;
t28 = t81 * t100 - t78 * t101;
t18 = qJ(5) * t52 + t28;
t27 = -t78 * t100 - t81 * t101;
t87 = -t53 * qJ(5) + t27;
t10 = t18 * t73 - t75 * t87;
t125 = t10 ^ 2;
t124 = t35 ^ 2;
t123 = t53 ^ 2;
t122 = 0.2e1 * t10;
t55 = -mrSges(7,1) * t80 + mrSges(7,2) * t77;
t121 = 0.2e1 * t55;
t120 = pkin(3) * t78;
t12 = t75 * t18 + t73 * t87;
t76 = cos(pkin(10));
t65 = -pkin(1) * t76 - pkin(2);
t54 = -pkin(3) * t82 + t65;
t38 = -pkin(4) * t52 + t54;
t13 = pkin(5) * t35 - pkin(9) * t37 + t38;
t3 = t12 * t80 + t13 * t77;
t119 = t3 * t80;
t116 = Ifges(7,4) * t77;
t115 = Ifges(7,4) * t80;
t114 = t10 * t35;
t113 = t37 * t77;
t66 = pkin(3) * t81 + pkin(4);
t44 = -t73 * t120 + t66 * t75;
t111 = t44 * mrSges(6,1);
t45 = t75 * t120 + t73 * t66;
t110 = t45 * mrSges(6,2);
t106 = Ifges(7,5) * t112 + Ifges(7,3) * t35;
t105 = t35 * mrSges(6,1) + t37 * mrSges(6,2);
t104 = Ifges(7,5) * t77 + Ifges(7,6) * t80;
t102 = t79 ^ 2 + t127;
t41 = pkin(9) + t45;
t99 = t103 * t41;
t62 = pkin(4) * t73 + pkin(9);
t98 = t103 * t62;
t97 = -t52 * mrSges(5,1) + t53 * mrSges(5,2);
t56 = Ifges(7,2) * t80 + t116;
t57 = Ifges(7,1) * t77 + t115;
t96 = t80 * t56 + t77 * t57 + Ifges(5,3) + Ifges(6,3);
t2 = -t12 * t77 + t13 * t80;
t95 = -t2 * t77 + t119;
t94 = -t82 * mrSges(4,1) + t79 * mrSges(4,2);
t93 = t75 * mrSges(6,1) - t73 * mrSges(6,2);
t92 = mrSges(7,1) * t77 + mrSges(7,2) * t80;
t90 = 0.2e1 * t128;
t89 = (mrSges(5,1) * t81 - mrSges(5,2) * t78) * pkin(3);
t88 = t128 * t37 + t35 * t55 - t105 - t97;
t7 = Ifges(7,6) * t35 + (-Ifges(7,2) * t77 + t115) * t37;
t8 = Ifges(7,5) * t35 + (Ifges(7,1) * t80 - t116) * t37;
t86 = -t28 * mrSges(5,2) - t12 * mrSges(6,2) + mrSges(7,3) * t119 - t2 * t109 - t56 * t113 / 0.2e1 + t57 * t112 / 0.2e1 + t27 * mrSges(5,1) + Ifges(6,5) * t37 + Ifges(5,6) * t52 + Ifges(5,5) * t53 + t77 * t8 / 0.2e1 + t80 * t7 / 0.2e1 + (t104 / 0.2e1 - Ifges(6,6)) * t35 + (t55 - mrSges(6,1)) * t10;
t64 = -pkin(4) * t75 - pkin(5);
t40 = -pkin(5) - t44;
t34 = t37 ^ 2;
t14 = t92 * t37;
t1 = [-0.2e1 * t27 * t53 * mrSges(5,3) + 0.2e1 * t65 * t94 + Ifges(4,2) * t127 + Ifges(5,1) * t123 + 0.2e1 * t54 * t97 + 0.2e1 * t38 * t105 + t14 * t122 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + Ifges(2,3) + Ifges(3,3) + (0.2e1 * t28 * mrSges(5,3) + 0.2e1 * Ifges(5,4) * t53 + Ifges(5,2) * t52) * t52 + (-0.2e1 * t12 * mrSges(6,3) + Ifges(6,2) * t35 + t106) * t35 + (mrSges(6,3) * t122 + Ifges(6,1) * t37 - t77 * t7 + t80 * t8 + (-Ifges(7,6) * t77 - (2 * Ifges(6,4))) * t35) * t37 + m(4) * (t102 * t63 ^ 2 + t65 ^ 2) + m(5) * (t27 ^ 2 + t28 ^ 2 + t54 ^ 2) + m(6) * (t12 ^ 2 + t38 ^ 2 + t125) + m(7) * (t2 ^ 2 + t3 ^ 2 + t125) + m(3) * (t74 ^ 2 + t76 ^ 2) * pkin(1) ^ 2 + (Ifges(4,1) * t79 + 0.2e1 * Ifges(4,4) * t82) * t79 + 0.2e1 * (t76 * mrSges(3,1) - t74 * mrSges(3,2)) * pkin(1) + 0.2e1 * t102 * t63 * mrSges(4,3); t35 * t14 + t91 * t37 + m(6) * (t12 * t37 + t114) + m(7) * (t95 * t37 + t114) + m(5) * (t27 * t52 + t28 * t53); m(3) + m(6) * (t34 + t124) + m(7) * (t103 * t34 + t124) + m(5) * (t52 ^ 2 + t123) + m(4) * t102; t91 * t41 + m(6) * (-t10 * t44 + t12 * t45) + (-mrSges(4,1) * t79 - mrSges(4,2) * t82) * t63 + (-t35 * t45 - t37 * t44) * mrSges(6,3) + (m(5) * (t27 * t81 + t28 * t78) + (t52 * t78 - t53 * t81) * mrSges(5,3)) * pkin(3) + m(7) * (t10 * t40 + t95 * t41) + Ifges(4,6) * t82 + Ifges(4,5) * t79 + t40 * t14 + t86; m(6) * (-t35 * t44 + t37 * t45) + m(7) * (t35 * t40 + t37 * t99) + m(5) * (t52 * t81 + t53 * t78) * pkin(3) + t88 - t94; 0.2e1 * t111 - 0.2e1 * t110 + t40 * t121 + Ifges(4,3) + 0.2e1 * t89 + t41 * t90 + m(7) * (t103 * t41 ^ 2 + t40 ^ 2) + m(6) * (t44 ^ 2 + t45 ^ 2) + m(5) * (t78 ^ 2 + t81 ^ 2) * pkin(3) ^ 2 + t96; (m(6) * (-t10 * t75 + t12 * t73) + (-t35 * t73 - t37 * t75) * mrSges(6,3)) * pkin(4) + t86 + (m(7) * t10 + t14) * t64 + (m(7) * t95 + t91) * t62; m(7) * (t64 * t35 + t37 * t98) + (-t35 * t75 + t37 * t73) * t126 + t88; m(7) * (t40 * t64 + t41 * t98) - t110 + t111 + (t40 + t64) * t55 + t89 + (m(6) * (t44 * t75 + t45 * t73) + t93) * pkin(4) + (t98 + t99) * mrSges(7,3) + t96; t64 * t121 + t62 * t90 + m(7) * (t103 * t62 ^ 2 + t64 ^ 2) + t96 + (0.2e1 * t93 + (t73 ^ 2 + t75 ^ 2) * t126) * pkin(4); t77 * t15 + t80 * t16 + m(7) * (t2 * t80 + t3 * t77) + m(6) * t38 + t105; 0; 0; 0; m(7) * t103 + m(6); mrSges(7,1) * t2 - mrSges(7,2) * t3 - Ifges(7,6) * t113 + t106; -t14; -t92 * t41 + t104; -t92 * t62 + t104; -t55; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
