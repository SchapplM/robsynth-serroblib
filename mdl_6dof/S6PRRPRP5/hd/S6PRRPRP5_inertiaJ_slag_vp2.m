% Calculate joint inertia matrix for
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2018-11-23 15:14
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:13:45
% EndTime: 2018-11-23 15:13:46
% DurationCPUTime: 0.71s
% Computational Cost: add. (619->217), mult. (1284->274), div. (0->0), fcn. (1059->8), ass. (0->88)
t121 = pkin(4) + pkin(8);
t120 = Ifges(7,2) + Ifges(6,3);
t64 = sin(pkin(6));
t71 = cos(qJ(2));
t104 = t64 * t71;
t68 = sin(qJ(2));
t105 = t64 * t68;
t65 = cos(pkin(6));
t67 = sin(qJ(3));
t70 = cos(qJ(3));
t18 = t67 * t105 - t65 * t70;
t66 = sin(qJ(5));
t69 = cos(qJ(5));
t10 = -t69 * t104 + t66 * t18;
t8 = t66 * t104 + t69 * t18;
t115 = t66 * t10 + t69 * t8;
t119 = t67 ^ 2 + t70 ^ 2;
t95 = -t66 ^ 2 - t69 ^ 2;
t118 = m(7) + m(6);
t101 = mrSges(5,1) + mrSges(4,3);
t117 = mrSges(6,3) + mrSges(7,2);
t114 = -m(7) * pkin(5) - mrSges(7,1);
t113 = t117 * t95;
t20 = t70 * t105 + t65 * t67;
t17 = t20 ^ 2;
t72 = -pkin(3) - pkin(9);
t111 = Ifges(6,4) * t66;
t110 = Ifges(6,4) * t69;
t109 = Ifges(7,5) * t66;
t108 = Ifges(7,5) * t69;
t107 = Ifges(6,6) * t66;
t106 = t18 * t67;
t103 = t66 * t70;
t102 = t69 * t70;
t88 = -t67 * qJ(4) - pkin(2);
t24 = t72 * t70 + t88;
t41 = t121 * t67;
t6 = t69 * t24 + t66 * t41;
t27 = t67 * mrSges(6,1) + mrSges(6,3) * t103;
t28 = -t67 * mrSges(7,1) - mrSges(7,2) * t103;
t100 = t27 - t28;
t29 = -t67 * mrSges(6,2) - mrSges(6,3) * t102;
t30 = -mrSges(7,2) * t102 + t67 * mrSges(7,3);
t99 = t29 + t30;
t36 = t66 * mrSges(6,1) + t69 * mrSges(6,2);
t98 = t36 + mrSges(5,3);
t97 = t95 * t72 ^ 2;
t96 = t119 * pkin(8) ^ 2;
t42 = t121 * t70;
t94 = t20 * qJ(4);
t92 = Ifges(7,6) * t102 + t120 * t67;
t91 = -m(5) * pkin(3) + mrSges(5,2);
t87 = t115 * t72;
t1 = t67 * qJ(6) + t6;
t5 = -t66 * t24 + t69 * t41;
t2 = -t67 * pkin(5) - t5;
t84 = t66 * t1 - t69 * t2;
t83 = t69 * t5 + t66 * t6;
t81 = t69 * mrSges(6,1) - t66 * mrSges(6,2);
t80 = t69 * mrSges(7,1) + t66 * mrSges(7,3);
t79 = pkin(5) * t69 + qJ(6) * t66;
t78 = (t20 * t70 + t106) * pkin(8);
t77 = -0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * t95;
t76 = -Ifges(6,6) * t69 + (-Ifges(7,4) - Ifges(6,5)) * t66;
t75 = m(7) * t79 + t80 + t81;
t73 = qJ(4) ^ 2;
t58 = t64 ^ 2;
t53 = Ifges(7,4) * t69;
t52 = Ifges(6,5) * t69;
t50 = Ifges(7,6) * t66;
t45 = t58 * t71 ^ 2;
t40 = Ifges(6,1) * t69 - t111;
t39 = Ifges(7,1) * t69 + t109;
t38 = -Ifges(6,2) * t66 + t110;
t37 = Ifges(7,3) * t66 + t108;
t35 = t66 * mrSges(7,1) - t69 * mrSges(7,3);
t34 = -t70 * mrSges(4,1) + t67 * mrSges(4,2);
t33 = t70 * mrSges(5,2) - t67 * mrSges(5,3);
t32 = -t70 * pkin(3) + t88;
t31 = t66 * pkin(5) - t69 * qJ(6) + qJ(4);
t23 = t81 * t70;
t22 = t80 * t70;
t15 = Ifges(6,5) * t67 + (-Ifges(6,1) * t66 - t110) * t70;
t14 = Ifges(7,4) * t67 + (-Ifges(7,1) * t66 + t108) * t70;
t13 = Ifges(6,6) * t67 + (-Ifges(6,2) * t69 - t111) * t70;
t12 = Ifges(7,6) * t67 + (Ifges(7,3) * t69 - t109) * t70;
t11 = t79 * t70 + t42;
t3 = [m(2) + m(3) * (t58 * t68 ^ 2 + t65 ^ 2 + t45) + t118 * (t10 ^ 2 + t8 ^ 2 + t17) + (m(5) + m(4)) * (t18 ^ 2 + t17 + t45); t100 * t8 + t101 * t106 + t99 * t10 + (-t68 * mrSges(3,2) + (mrSges(3,1) - t33 - t34) * t71) * t64 + (t101 * t70 + t22 + t23) * t20 + m(7) * (t1 * t10 + t11 * t20 - t2 * t8) + m(6) * (t6 * t10 + t42 * t20 + t5 * t8) + m(5) * (-t32 * t104 + t78) + m(4) * (pkin(2) * t104 + t78); -0.2e1 * pkin(2) * t34 + 0.2e1 * t1 * t30 + 0.2e1 * t11 * t22 + 0.2e1 * t2 * t28 + 0.2e1 * t42 * t23 + 0.2e1 * t5 * t27 + 0.2e1 * t6 * t29 + 0.2e1 * t32 * t33 + Ifges(3,3) + ((Ifges(4,1) + Ifges(5,2)) * t67 + t92) * t67 + m(5) * (t32 ^ 2 + t96) + m(4) * (pkin(2) ^ 2 + t96) + m(6) * (t42 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t1 ^ 2 + t11 ^ 2 + t2 ^ 2) + ((Ifges(4,2) + Ifges(5,3)) * t70 + (t12 - t13) * t69 + (-t14 - t15) * t66 + ((2 * Ifges(4,4)) + (2 * Ifges(5,6)) + t76) * t67) * t70 + 0.2e1 * t101 * pkin(8) * t119; (-mrSges(4,1) + mrSges(5,2)) * t18 + (-mrSges(4,2) + t35 + t98) * t20 + m(7) * (t31 * t20 + t87) + m(6) * (t87 + t94) + m(5) * (-pkin(3) * t18 + t94) - t117 * t115; qJ(4) * t23 + t11 * t35 + t31 * t22 + t42 * t36 + (t2 * mrSges(7,2) - t5 * mrSges(6,3) + t14 / 0.2e1 + t15 / 0.2e1 + t100 * t72) * t69 + (-t1 * mrSges(7,2) - t6 * mrSges(6,3) + t12 / 0.2e1 - t13 / 0.2e1 + t99 * t72) * t66 + m(7) * (t31 * t11 + t84 * t72) + m(6) * (qJ(4) * t42 + t83 * t72) + (-pkin(3) * mrSges(5,1) - t107 / 0.2e1 + t52 / 0.2e1 + t53 / 0.2e1 + t50 / 0.2e1 - Ifges(5,4) + Ifges(4,5) + (-mrSges(4,1) + t91) * pkin(8)) * t67 + (qJ(4) * mrSges(5,1) - Ifges(5,5) + Ifges(4,6) + (t37 / 0.2e1 - t38 / 0.2e1) * t69 + (-t39 / 0.2e1 - t40 / 0.2e1) * t66 + (m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3)) * pkin(8)) * t70; -0.2e1 * pkin(3) * mrSges(5,2) + 0.2e1 * t31 * t35 + Ifges(5,1) + Ifges(4,3) + (t40 + t39) * t69 + (t37 - t38) * t66 + 0.2e1 * t98 * qJ(4) + m(6) * (t73 - t97) + m(7) * (t31 ^ 2 - t97) + m(5) * (pkin(3) ^ 2 + t73) + 0.2e1 * t72 * t113; m(5) * t18 + t115 * t118; t100 * t69 + (m(5) * pkin(8) + mrSges(5,1)) * t67 + t99 * t66 + m(7) * t84 + m(6) * t83; t72 * t77 + t113 + t91; m(5) + t77; (mrSges(6,1) - t114) * t8 + (m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3)) * t10; -pkin(5) * t28 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + qJ(6) * t30 + t1 * mrSges(7,3) + t5 * mrSges(6,1) - t2 * mrSges(7,1) - t6 * mrSges(6,2) + t76 * t70 + t92; -t79 * mrSges(7,2) + t75 * t72 - t107 + t50 + t52 + t53; t75; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t120; -m(7) * t8; m(7) * t2 + t28; (-m(7) * t72 + mrSges(7,2)) * t69; -m(7) * t69; t114; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t3(1) t3(2) t3(4) t3(7) t3(11) t3(16); t3(2) t3(3) t3(5) t3(8) t3(12) t3(17); t3(4) t3(5) t3(6) t3(9) t3(13) t3(18); t3(7) t3(8) t3(9) t3(10) t3(14) t3(19); t3(11) t3(12) t3(13) t3(14) t3(15) t3(20); t3(16) t3(17) t3(18) t3(19) t3(20) t3(21);];
Mq  = res;
