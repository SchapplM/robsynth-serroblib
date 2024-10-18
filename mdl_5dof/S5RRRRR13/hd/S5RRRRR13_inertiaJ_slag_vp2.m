% Calculate joint inertia matrix for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR13_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR13_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR13_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:22
% EndTime: 2024-09-27 17:30:22
% DurationCPUTime: 0.12s
% Computational Cost: add. (1298->175), mult. (2750->241), div. (0->0), fcn. (2572->10), ass. (0->100)
t87 = sin(pkin(5));
t89 = sin(qJ(5));
t90 = sin(qJ(4));
t93 = cos(qJ(5));
t94 = cos(qJ(4));
t53 = (-t89 * t90 + t93 * t94) * t87;
t54 = (t89 * t94 + t90 * t93) * t87;
t132 = Ifges(6,5) * t54 + Ifges(6,6) * t53;
t88 = cos(pkin(5));
t107 = Ifges(6,3) * t88 + t132;
t112 = t87 * t94;
t113 = t87 * t90;
t131 = Ifges(5,5) * t113 + Ifges(5,6) * t112;
t133 = Ifges(5,3) * t88 + t107 + t131;
t86 = t87 ^ 2;
t122 = pkin(10) * t87;
t110 = t88 * t94;
t92 = sin(qJ(2));
t124 = pkin(1) * t92;
t96 = cos(qJ(2));
t80 = t96 * pkin(1) + pkin(2);
t91 = sin(qJ(3));
t95 = cos(qJ(3));
t58 = -t91 * t124 + t95 * t80;
t56 = pkin(3) + t58;
t47 = t56 * t110;
t59 = t95 * t124 + t91 * t80;
t84 = t87 * pkin(9);
t50 = t84 + t59;
t85 = t88 * pkin(4);
t18 = t47 + t85 + (-t50 - t122) * t90;
t111 = t88 * t90;
t27 = t56 * t111 + t94 * t50;
t74 = pkin(10) * t112;
t23 = t74 + t27;
t4 = t93 * t18 - t89 * t23;
t5 = t89 * t18 + t93 * t23;
t130 = t4 * mrSges(6,1) - t5 * mrSges(6,2);
t119 = t95 * pkin(2);
t79 = pkin(3) + t119;
t68 = t79 * t110;
t120 = t91 * pkin(2);
t69 = t84 + t120;
t31 = t68 + t85 + (-t69 - t122) * t90;
t44 = t79 * t111 + t94 * t69;
t32 = t74 + t44;
t11 = t93 * t31 - t89 * t32;
t12 = t89 * t31 + t93 * t32;
t129 = t11 * mrSges(6,1) - t12 * mrSges(6,2);
t77 = pkin(3) * t110;
t38 = t77 + t85 + (-pkin(9) - pkin(10)) * t113;
t62 = pkin(3) * t111 + pkin(9) * t112;
t49 = t74 + t62;
t21 = t93 * t38 - t89 * t49;
t22 = t89 * t38 + t93 * t49;
t128 = t21 * mrSges(6,1) - t22 * mrSges(6,2);
t28 = -t53 * mrSges(6,1) + t54 * mrSges(6,2);
t123 = pkin(4) * t94;
t66 = (-pkin(3) - t123) * t87;
t24 = t66 * t28;
t61 = -pkin(9) * t113 + t77;
t63 = t88 * mrSges(5,1) - mrSges(5,3) * t113;
t36 = t61 * t63;
t64 = -t88 * mrSges(5,2) + mrSges(5,3) * t112;
t37 = t62 * t64;
t42 = t88 * mrSges(6,1) - t54 * mrSges(6,3);
t8 = t21 * t42;
t41 = -t88 * mrSges(6,2) + t53 * mrSges(6,3);
t9 = t22 * t41;
t127 = t24 + t36 + t37 + t8 + t9;
t1 = t4 * t42;
t35 = (-t56 - t123) * t87;
t13 = t35 * t28;
t26 = -t90 * t50 + t47;
t16 = t26 * t63;
t17 = t27 * t64;
t2 = t5 * t41;
t55 = t58 * mrSges(4,1);
t126 = t1 + t13 + t16 + t17 + t2 + t55;
t125 = m(6) * pkin(4);
t116 = t59 * mrSges(4,2);
t115 = t79 * t86;
t114 = t86 * (-mrSges(5,1) * t94 + mrSges(5,2) * t90);
t109 = -0.2e1 * t114;
t108 = mrSges(4,2) * t120;
t105 = (t96 * mrSges(3,1) - t92 * mrSges(3,2)) * pkin(1);
t104 = (t93 * mrSges(6,1) - t89 * mrSges(6,2)) * pkin(4);
t103 = Ifges(6,1) * t54 ^ 2 + Ifges(4,3) + (0.2e1 * Ifges(6,4) * t54 + Ifges(6,2) * t53) * t53 + ((Ifges(5,1) * t90 + Ifges(5,4) * t94) * t113 + (Ifges(5,4) * t90 + Ifges(5,2) * t94) * t112) * t87 + (t131 + t132 + t133) * t88;
t102 = Ifges(3,3) + t103;
t101 = (t41 * t89 + t42 * t93) * pkin(4) + t133;
t57 = (-t79 - t123) * t87;
t19 = t57 * t28;
t43 = -t90 * t69 + t68;
t29 = t43 * t63;
t30 = t44 * t64;
t6 = t11 * t42;
t7 = t12 * t41;
t81 = mrSges(4,1) * t119;
t100 = t19 + t29 + t30 + t6 + t7 + t81 + t103;
t3 = [Ifges(2,3) + t102 + 0.2e1 * t105 + m(3) * (t92 ^ 2 + t96 ^ 2) * pkin(1) ^ 2 + m(6) * (t35 ^ 2 + t4 ^ 2 + t5 ^ 2) + m(5) * (t86 * t56 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(4) * (t58 ^ 2 + t59 ^ 2) + 0.2e1 * t55 + 0.2e1 * t17 + 0.2e1 * t16 + 0.2e1 * t13 + 0.2e1 * t2 + t56 * t109 - 0.2e1 * t116 + 0.2e1 * t1; Ifges(3,3) + t100 + m(5) * (t56 * t115 + t43 * t26 + t44 * t27) + m(6) * (t11 * t4 + t12 * t5 + t57 * t35) + t105 + (-t59 - t120) * mrSges(4,2) + m(4) * (t58 * t95 + t59 * t91) * pkin(2) + (-t56 - t79) * t114 + t126; t102 + m(4) * (t91 ^ 2 + t95 ^ 2) * pkin(2) ^ 2 + m(6) * (t11 ^ 2 + t12 ^ 2 + t57 ^ 2) + m(5) * (t86 * t79 ^ 2 + t43 ^ 2 + t44 ^ 2) + 0.2e1 * t81 + 0.2e1 * t30 + 0.2e1 * t29 + 0.2e1 * t19 + 0.2e1 * t7 + 0.2e1 * t6 - 0.2e1 * t108 + t79 * t109; t103 + m(6) * (t21 * t4 + t22 * t5 + t66 * t35) + m(5) * (t86 * pkin(3) * t56 + t61 * t26 + t62 * t27) + (-pkin(3) - t56) * t114 - t116 + t126 + t127; t100 + m(6) * (t21 * t11 + t22 * t12 + t66 * t57) + m(5) * (pkin(3) * t115 + t61 * t43 + t62 * t44) + (-pkin(3) - t79) * t114 - t108 + t127; pkin(3) * t109 + 0.2e1 * t24 + 0.2e1 * t36 + 0.2e1 * t37 + 0.2e1 * t8 + 0.2e1 * t9 + m(6) * (t21 ^ 2 + t22 ^ 2 + t66 ^ 2) + m(5) * (t86 * pkin(3) ^ 2 + t61 ^ 2 + t62 ^ 2) + t103; t26 * mrSges(5,1) - t27 * mrSges(5,2) + (t4 * t93 + t5 * t89) * t125 + t101 + t130; t43 * mrSges(5,1) - t44 * mrSges(5,2) + (t11 * t93 + t12 * t89) * t125 + t101 + t129; t61 * mrSges(5,1) - t62 * mrSges(5,2) + (t21 * t93 + t22 * t89) * t125 + t101 + t128; Ifges(5,3) + Ifges(6,3) + m(6) * (t89 ^ 2 + t93 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t104; t107 + t130; t107 + t129; t107 + t128; Ifges(6,3) + t104; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
