% Calculate joint inertia matrix for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2018-11-23 16:13
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:13:36
% EndTime: 2018-11-23 16:13:37
% DurationCPUTime: 0.81s
% Computational Cost: add. (1083->256), mult. (2064->323), div. (0->0), fcn. (2049->6), ass. (0->93)
t120 = Ifges(6,2) + Ifges(5,3);
t77 = sin(qJ(4));
t79 = cos(qJ(4));
t119 = t77 ^ 2 + t79 ^ 2;
t118 = 0.2e1 * t119;
t112 = cos(qJ(3));
t100 = pkin(7) + qJ(2);
t76 = cos(pkin(9));
t46 = t100 * t76;
t78 = sin(qJ(3));
t75 = sin(pkin(9));
t90 = t100 * t75;
t28 = t112 * t90 + t46 * t78;
t117 = t28 ^ 2;
t72 = t76 ^ 2;
t116 = -2 * mrSges(7,3);
t115 = 0.2e1 * t28;
t60 = -pkin(2) * t76 - pkin(1);
t114 = 0.2e1 * t60;
t113 = m(6) + m(7);
t80 = -pkin(4) - pkin(5);
t111 = Ifges(5,4) * t77;
t110 = Ifges(5,4) * t79;
t109 = Ifges(7,4) * t77;
t108 = Ifges(7,4) * t79;
t107 = Ifges(6,5) * t77;
t106 = Ifges(6,5) * t79;
t43 = -t112 * t76 + t75 * t78;
t105 = Ifges(7,5) * t43;
t44 = t112 * t75 + t78 * t76;
t104 = t44 * t77;
t103 = t44 * t79;
t102 = mrSges(6,2) - mrSges(7,3);
t101 = -Ifges(5,6) - Ifges(7,6);
t99 = pkin(8) - qJ(6);
t21 = -mrSges(5,2) * t43 - mrSges(5,3) * t104;
t25 = -mrSges(6,2) * t104 + mrSges(6,3) * t43;
t98 = t21 + t25;
t23 = mrSges(5,1) * t43 - mrSges(5,3) * t103;
t24 = -t43 * mrSges(6,1) + mrSges(6,2) * t103;
t97 = -t23 + t24;
t19 = pkin(3) * t43 - pkin(8) * t44 + t60;
t30 = t112 * t46 - t78 * t90;
t7 = t77 * t19 + t79 * t30;
t62 = t77 * qJ(5);
t96 = t79 * pkin(4) + t62;
t49 = t79 * mrSges(7,1) + t77 * mrSges(7,2);
t95 = t119 * pkin(8) ^ 2;
t94 = t75 ^ 2 + t72;
t92 = qJ(5) * t79;
t91 = qJ(6) * t44;
t45 = -pkin(3) - t96;
t3 = t43 * qJ(5) + t7;
t89 = -t76 * mrSges(3,1) + t75 * mrSges(3,2);
t26 = t77 * t30;
t6 = t19 * t79 - t26;
t17 = -mrSges(7,1) * t104 + mrSges(7,2) * t103;
t22 = -t43 * mrSges(7,1) - mrSges(7,3) * t103;
t87 = Ifges(6,6) * t104 + t120 * t43 + (Ifges(6,4) + Ifges(5,5)) * t103;
t86 = t77 * mrSges(5,1) + t79 * mrSges(5,2);
t85 = t77 * mrSges(6,1) - t79 * mrSges(6,3);
t84 = pkin(4) * t77 - t92;
t81 = qJ(5) ^ 2;
t67 = Ifges(6,4) * t77;
t66 = Ifges(5,5) * t77;
t65 = Ifges(5,6) * t79;
t57 = Ifges(5,1) * t77 + t110;
t56 = Ifges(6,1) * t77 - t106;
t55 = Ifges(7,1) * t77 - t108;
t54 = Ifges(5,2) * t79 + t111;
t53 = -Ifges(7,2) * t79 + t109;
t52 = -Ifges(6,3) * t79 + t107;
t51 = t99 * t79;
t50 = -t79 * mrSges(5,1) + t77 * mrSges(5,2);
t48 = -t79 * mrSges(6,1) - t77 * mrSges(6,3);
t47 = t99 * t77;
t41 = pkin(5) * t79 - t45;
t38 = t44 * mrSges(4,2);
t20 = mrSges(7,2) * t43 + mrSges(7,3) * t104;
t18 = t86 * t44;
t16 = t85 * t44;
t14 = Ifges(5,5) * t43 + (Ifges(5,1) * t79 - t111) * t44;
t13 = Ifges(6,4) * t43 + (Ifges(6,1) * t79 + t107) * t44;
t12 = -t105 + (Ifges(7,1) * t79 + t109) * t44;
t11 = Ifges(5,6) * t43 + (-Ifges(5,2) * t77 + t110) * t44;
t10 = -Ifges(7,6) * t43 + (Ifges(7,2) * t77 + t108) * t44;
t9 = Ifges(6,6) * t43 + (Ifges(6,3) * t77 + t106) * t44;
t8 = t84 * t44 + t28;
t5 = (t80 * t77 + t92) * t44 - t28;
t4 = -pkin(4) * t43 - t6;
t2 = t77 * t91 + t3;
t1 = t26 + (-t19 - t91) * t79 + t80 * t43;
t15 = [Ifges(3,2) * t72 - 0.2e1 * pkin(1) * t89 + t38 * t114 + 0.2e1 * t2 * t20 + 0.2e1 * t7 * t21 + 0.2e1 * t1 * t22 + 0.2e1 * t6 * t23 + 0.2e1 * t4 * t24 + 0.2e1 * t3 * t25 + t18 * t115 + 0.2e1 * t8 * t16 + 0.2e1 * t5 * t17 + Ifges(2,3) + (Ifges(3,1) * t75 + 0.2e1 * Ifges(3,4) * t76) * t75 + 0.2e1 * t94 * qJ(2) * mrSges(3,3) + (mrSges(4,1) * t114 - 0.2e1 * t30 * mrSges(4,3) + (Ifges(7,3) + Ifges(4,2)) * t43 + t87) * t43 + m(3) * (t94 * qJ(2) ^ 2 + pkin(1) ^ 2) + m(4) * (t30 ^ 2 + t60 ^ 2 + t117) + m(5) * (t6 ^ 2 + t7 ^ 2 + t117) + m(6) * (t3 ^ 2 + t4 ^ 2 + t8 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + (mrSges(4,3) * t115 + Ifges(4,1) * t44 - 0.2e1 * Ifges(4,4) * t43 + (t12 + t13 + t14 - t105) * t79 + (t101 * t43 + t10 - t11 + t9) * t77) * t44; -m(3) * pkin(1) + t43 * mrSges(4,1) + t38 + (-t22 - t97) * t79 + (t20 + t98) * t77 + m(6) * (t3 * t77 - t4 * t79) + m(7) * (-t1 * t79 + t2 * t77) + m(5) * (t6 * t79 + t7 * t77) + m(4) * t60 + t89; m(3) + m(4) + (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t118; -t30 * mrSges(4,2) - pkin(3) * t18 + t45 * t16 + t41 * t17 + t51 * t20 + t47 * t22 + t8 * t48 + t5 * t49 + (t67 / 0.2e1 + t66 / 0.2e1 + t65 / 0.2e1 - Ifges(4,6)) * t43 + (t50 - mrSges(4,1)) * t28 + (-t9 / 0.2e1 - t10 / 0.2e1 + t11 / 0.2e1 + t7 * mrSges(5,3) - t2 * mrSges(7,3) + t3 * mrSges(6,2) + (Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t43 + t98 * pkin(8)) * t79 + (t12 / 0.2e1 + t13 / 0.2e1 + t14 / 0.2e1 - t105 / 0.2e1 - t6 * mrSges(5,3) - t1 * mrSges(7,3) + t4 * mrSges(6,2) + t97 * pkin(8)) * t77 + m(5) * (-pkin(3) * t28 + (-t6 * t77 + t7 * t79) * pkin(8)) + m(6) * (t45 * t8 + (t3 * t79 + t4 * t77) * pkin(8)) + m(7) * (t1 * t47 + t2 * t51 + t41 * t5) + (Ifges(4,5) + (t55 / 0.2e1 + t56 / 0.2e1 + t57 / 0.2e1) * t79 + (t52 / 0.2e1 + t53 / 0.2e1 - t54 / 0.2e1) * t77) * t44; m(7) * (-t47 * t79 + t51 * t77); -0.2e1 * pkin(3) * t50 + 0.2e1 * t41 * t49 + 0.2e1 * t45 * t48 + Ifges(4,3) + m(6) * (t45 ^ 2 + t95) + m(7) * (t41 ^ 2 + t47 ^ 2 + t51 ^ 2) + m(5) * (pkin(3) ^ 2 + t95) + (t51 * t116 - t52 - t53 + t54) * t79 + (t47 * t116 + t55 + t56 + t57) * t77 + (mrSges(6,2) + mrSges(5,3)) * pkin(8) * t118; t6 * mrSges(5,1) - t4 * mrSges(6,1) - t1 * mrSges(7,1) - t7 * mrSges(5,2) + t2 * mrSges(7,2) + t3 * mrSges(6,3) + Ifges(7,3) * t43 - pkin(4) * t24 + t80 * t22 + (t20 + t25) * qJ(5) + m(6) * (-pkin(4) * t4 + qJ(5) * t3) + m(7) * (qJ(5) * t2 + t1 * t80) + (-Ifges(7,5) * t79 + t101 * t77) * t44 + t87; (mrSges(5,1) + mrSges(6,1)) * t79 + (-mrSges(5,2) + mrSges(6,3)) * t77 + m(6) * t96 + m(7) * (-t79 * t80 + t62) + t49; m(7) * (qJ(5) * t51 + t47 * t80) - t47 * mrSges(7,1) + t51 * mrSges(7,2) + t67 + t66 + t65 + (-pkin(4) * mrSges(6,2) - t80 * mrSges(7,3) - Ifges(7,5)) * t77 + (t102 * qJ(5) - Ifges(6,6) + Ifges(7,6)) * t79 + (-m(6) * t84 - t85 - t86) * pkin(8); 0.2e1 * pkin(4) * mrSges(6,1) - 0.2e1 * t80 * mrSges(7,1) + Ifges(7,3) + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * qJ(5) + m(6) * (pkin(4) ^ 2 + t81) + m(7) * (t80 ^ 2 + t81) + t120; m(6) * t4 + m(7) * t1 + t22 + t24; -t113 * t79; m(7) * t47 + (m(6) * pkin(8) + t102) * t77; -m(6) * pkin(4) + m(7) * t80 - mrSges(6,1) - mrSges(7,1); t113; m(7) * t5 + t17; 0; m(7) * t41 + t49; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
