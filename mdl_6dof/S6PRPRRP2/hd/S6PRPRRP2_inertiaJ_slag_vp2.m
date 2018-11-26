% Calculate joint inertia matrix for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2018-11-23 15:00
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:00:09
% EndTime: 2018-11-23 15:00:09
% DurationCPUTime: 0.74s
% Computational Cost: add. (676->210), mult. (1576->286), div. (0->0), fcn. (1458->10), ass. (0->92)
t70 = sin(qJ(5));
t73 = cos(qJ(5));
t92 = t70 ^ 2 + t73 ^ 2;
t66 = sin(pkin(11));
t67 = sin(pkin(6));
t68 = cos(pkin(11));
t72 = sin(qJ(2));
t75 = cos(qJ(2));
t19 = (t66 * t75 + t68 * t72) * t67;
t69 = cos(pkin(6));
t71 = sin(qJ(4));
t74 = cos(qJ(4));
t14 = t74 * t19 + t69 * t71;
t17 = (t66 * t72 - t68 * t75) * t67;
t4 = t70 * t14 - t73 * t17;
t6 = t73 * t14 + t70 * t17;
t123 = t4 * t70 + t6 * t73;
t52 = t66 * pkin(2) + pkin(8);
t122 = 0.2e1 * t52;
t121 = m(6) + m(7);
t120 = mrSges(7,2) + mrSges(6,3);
t119 = -m(7) * pkin(5) - mrSges(7,1);
t118 = t120 * t92;
t12 = t71 * t19 - t69 * t74;
t11 = t12 ^ 2;
t117 = t17 ^ 2;
t114 = t71 * pkin(9);
t113 = t74 * pkin(4);
t112 = Ifges(6,4) * t70;
t111 = Ifges(6,4) * t73;
t110 = Ifges(7,5) * t70;
t109 = Ifges(7,5) * t73;
t108 = Ifges(6,6) * t74;
t107 = Ifges(7,6) * t74;
t106 = t12 * t71;
t105 = t14 * t74;
t104 = t52 * t74;
t103 = t70 * t71;
t102 = t71 * t73;
t53 = -t68 * pkin(2) - pkin(3);
t28 = -t113 + t53 - t114;
t101 = t73 * t28;
t100 = t74 * t12;
t99 = -Ifges(6,3) - Ifges(7,2);
t82 = t70 * mrSges(7,1) - t73 * mrSges(7,3);
t26 = t82 * t71;
t83 = t70 * mrSges(6,1) + t73 * mrSges(6,2);
t27 = t83 * t71;
t98 = -t26 - t27;
t29 = t74 * mrSges(6,2) - mrSges(6,3) * t103;
t32 = -mrSges(7,2) * t103 - t74 * mrSges(7,3);
t97 = t29 + t32;
t30 = -t74 * mrSges(6,1) - mrSges(6,3) * t102;
t31 = t74 * mrSges(7,1) + mrSges(7,2) * t102;
t96 = -t30 + t31;
t10 = t73 * t104 + t70 * t28;
t36 = -t73 * mrSges(6,1) + t70 * mrSges(6,2);
t95 = t36 - mrSges(5,1);
t94 = t92 * t114;
t93 = t92 * pkin(9) ^ 2;
t63 = t71 ^ 2;
t65 = t74 ^ 2;
t91 = t63 + t65;
t35 = -t73 * mrSges(7,1) - t70 * mrSges(7,3);
t90 = -t35 - t95;
t89 = -Ifges(7,6) * t103 + (-Ifges(7,4) - Ifges(6,5)) * t102;
t87 = t123 * pkin(9);
t9 = -t70 * t104 + t101;
t84 = t10 * t73 - t70 * t9;
t81 = -pkin(5) * t70 + qJ(6) * t73;
t7 = -t74 * qJ(6) + t10;
t8 = -t101 + (t52 * t70 + pkin(5)) * t74;
t79 = m(7) * (t7 * t73 + t70 * t8) + t97 * t73 + t96 * t70;
t78 = m(7) * t81 - t82 - t83;
t61 = t69 ^ 2;
t57 = Ifges(7,4) * t70;
t56 = Ifges(6,5) * t70;
t55 = Ifges(6,6) * t73;
t50 = t52 ^ 2;
t42 = t63 * t50;
t41 = Ifges(6,1) * t70 + t111;
t40 = Ifges(7,1) * t70 - t109;
t39 = Ifges(6,2) * t73 + t112;
t38 = -Ifges(7,3) * t73 + t110;
t37 = -t74 * mrSges(5,1) + t71 * mrSges(5,2);
t34 = -t73 * pkin(5) - t70 * qJ(6) - pkin(4);
t23 = -Ifges(6,5) * t74 + (Ifges(6,1) * t73 - t112) * t71;
t22 = -Ifges(7,4) * t74 + (Ifges(7,1) * t73 + t110) * t71;
t21 = -t108 + (-Ifges(6,2) * t70 + t111) * t71;
t20 = -t107 + (Ifges(7,3) * t70 + t109) * t71;
t15 = (t52 - t81) * t71;
t1 = [m(2) + m(5) * (t14 ^ 2 + t11 + t117) + m(4) * (t19 ^ 2 + t117 + t61) + m(3) * (t61 + (t72 ^ 2 + t75 ^ 2) * t67 ^ 2) + t121 * (t4 ^ 2 + t6 ^ 2 + t11); mrSges(5,3) * t105 - t19 * mrSges(4,2) + (t75 * mrSges(3,1) - t72 * mrSges(3,2)) * t67 + t97 * t6 + t96 * t4 + (-mrSges(4,1) + t37) * t17 + (t71 * mrSges(5,3) - t98) * t12 + m(6) * (t10 * t6 + t52 * t106 - t9 * t4) + m(7) * (t15 * t12 + t8 * t4 + t7 * t6) + m(5) * (t53 * t17 + (t105 + t106) * t52) + m(4) * (-t17 * t68 + t19 * t66) * pkin(2); 0.2e1 * t10 * t29 + 0.2e1 * t15 * t26 + 0.2e1 * t9 * t30 + 0.2e1 * t8 * t31 + 0.2e1 * t7 * t32 + 0.2e1 * t53 * t37 + Ifges(3,3) + Ifges(4,3) + ((Ifges(5,2) - t99) * t74 + t89) * t74 + m(5) * (t65 * t50 + t53 ^ 2 + t42) + m(6) * (t10 ^ 2 + t9 ^ 2 + t42) + m(7) * (t15 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(4) * (t66 ^ 2 + t68 ^ 2) * pkin(2) ^ 2 + (Ifges(5,1) * t71 + 0.2e1 * Ifges(5,4) * t74 + t27 * t122 + (t22 + t23) * t73 + (t20 - t21 + t108) * t70) * t71 + 0.2e1 * (t68 * mrSges(4,1) - t66 * mrSges(4,2)) * pkin(2) + t91 * mrSges(5,3) * t122; m(4) * t69 + m(5) * (t71 * t14 - t100) + t121 * (t6 * t102 + t4 * t103 - t100); (-m(7) * t15 + t98) * t74 + (m(6) * (t84 - t104) + t79) * t71; m(4) + m(5) * t91 + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t92 * t63 + t65); -t14 * mrSges(5,2) - t90 * t12 + m(6) * (-pkin(4) * t12 + t87) + m(7) * (t34 * t12 + t87) + t120 * t123; -pkin(4) * t27 + t34 * t26 + (m(7) * t34 + t35) * t15 + (-t52 * mrSges(5,2) - t57 / 0.2e1 - t56 / 0.2e1 - t55 / 0.2e1 + Ifges(5,6)) * t74 + (t7 * mrSges(7,2) + t10 * mrSges(6,3) - t20 / 0.2e1 + t21 / 0.2e1 + t107 / 0.2e1) * t73 + (t22 / 0.2e1 + t23 / 0.2e1 + t8 * mrSges(7,2) - t9 * mrSges(6,3)) * t70 + (m(6) * t84 + t79) * pkin(9) + (Ifges(5,5) + (t40 / 0.2e1 + t41 / 0.2e1) * t73 + (t38 / 0.2e1 - t39 / 0.2e1) * t70 + (-m(6) * pkin(4) + t95) * t52) * t71; t90 * t74 + m(7) * (-t34 * t74 + t94) + m(6) * (t94 + t113) + (-mrSges(5,2) + t118) * t71; -0.2e1 * pkin(4) * t36 + 0.2e1 * t34 * t35 + Ifges(5,3) + (-t38 + t39) * t73 + (t40 + t41) * t70 + m(7) * (t34 ^ 2 + t93) + m(6) * (pkin(4) ^ 2 + t93) + 0.2e1 * pkin(9) * t118; (m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3)) * t6 + (-mrSges(6,1) + t119) * t4; -Ifges(6,6) * t103 - pkin(5) * t31 + m(7) * (-pkin(5) * t8 + qJ(6) * t7) + qJ(6) * t32 + t7 * mrSges(7,3) - t10 * mrSges(6,2) + t9 * mrSges(6,1) - t8 * mrSges(7,1) + t99 * t74 - t89; t78 * t71; t81 * mrSges(7,2) - Ifges(7,6) * t73 + t78 * pkin(9) + t55 + t56 + t57; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) - t99; m(7) * t4; m(7) * t8 + t31; m(7) * t103; (m(7) * pkin(9) + mrSges(7,2)) * t70; t119; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
