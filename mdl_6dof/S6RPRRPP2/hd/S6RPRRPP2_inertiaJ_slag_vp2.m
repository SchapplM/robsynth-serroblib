% Calculate joint inertia matrix for
% S6RPRRPP2
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
% Datum: 2018-11-23 16:12
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:12:03
% EndTime: 2018-11-23 16:12:04
% DurationCPUTime: 0.89s
% Computational Cost: add. (687->256), mult. (1330->315), div. (0->0), fcn. (977->6), ass. (0->98)
t70 = sin(qJ(4));
t72 = cos(qJ(4));
t90 = t70 ^ 2 + t72 ^ 2;
t68 = sin(pkin(9));
t52 = pkin(1) * t68 + pkin(7);
t121 = 0.2e1 * t52;
t71 = sin(qJ(3));
t79 = t70 * mrSges(6,1) - t72 * mrSges(6,3);
t17 = t79 * t71;
t111 = pkin(4) * t70;
t88 = qJ(5) * t72;
t43 = t71 * t88;
t8 = -t43 + (t52 + t111) * t71;
t120 = m(6) * t8 + t17;
t101 = t70 * t71;
t100 = t71 * t72;
t44 = mrSges(7,2) * t100;
t18 = -mrSges(7,1) * t101 + t44;
t114 = pkin(4) + pkin(5);
t5 = t43 + (-t114 * t70 - t52) * t71;
t119 = m(7) * t5 + t18;
t118 = -2 * mrSges(7,3);
t115 = m(6) + m(7);
t73 = cos(qJ(3));
t113 = m(7) * t73;
t112 = pkin(3) * t73;
t110 = pkin(8) * t71;
t109 = Ifges(5,4) * t70;
t108 = Ifges(5,4) * t72;
t107 = Ifges(7,4) * t70;
t106 = Ifges(7,4) * t72;
t105 = Ifges(6,5) * t70;
t104 = Ifges(6,5) * t72;
t103 = Ifges(7,5) * t73;
t102 = t52 * t73;
t34 = -mrSges(5,1) * t72 + mrSges(5,2) * t70;
t99 = mrSges(4,1) - t34;
t98 = -mrSges(6,1) - mrSges(7,1);
t97 = mrSges(6,2) - mrSges(7,3);
t96 = Ifges(5,6) + Ifges(7,6);
t95 = pkin(8) - qJ(6);
t23 = mrSges(5,2) * t73 - mrSges(5,3) * t101;
t27 = -mrSges(6,2) * t101 - mrSges(6,3) * t73;
t94 = t23 + t27;
t25 = -mrSges(5,1) * t73 - mrSges(5,3) * t100;
t26 = t73 * mrSges(6,1) + mrSges(6,2) * t100;
t93 = -t25 + t26;
t69 = cos(pkin(9));
t53 = -pkin(1) * t69 - pkin(2);
t20 = -t110 + t53 - t112;
t7 = t72 * t102 + t70 * t20;
t92 = t90 * t110;
t33 = t72 * mrSges(7,1) + t70 * mrSges(7,2);
t91 = t90 * pkin(8) ^ 2;
t65 = t71 ^ 2;
t67 = t73 ^ 2;
t89 = t65 + t67;
t87 = qJ(6) * t71;
t85 = Ifges(5,3) + Ifges(7,3) + Ifges(6,2);
t84 = -Ifges(6,6) * t101 + (-Ifges(6,4) - Ifges(5,5)) * t100;
t83 = qJ(5) * t70 + pkin(3);
t28 = t70 * t102;
t6 = t20 * t72 - t28;
t24 = t73 * mrSges(7,1) - mrSges(7,3) * t100;
t3 = -qJ(5) * t73 + t7;
t81 = -t6 * t70 + t7 * t72;
t80 = t70 * mrSges(5,1) + t72 * mrSges(5,2);
t63 = t73 * pkin(4);
t4 = -t6 + t63;
t78 = m(6) * (t3 * t72 + t4 * t70);
t75 = qJ(5) ^ 2;
t60 = Ifges(6,4) * t70;
t59 = Ifges(5,5) * t70;
t58 = Ifges(5,6) * t72;
t51 = t52 ^ 2;
t42 = t65 * t51;
t41 = Ifges(5,1) * t70 + t108;
t40 = Ifges(6,1) * t70 - t104;
t39 = Ifges(7,1) * t70 - t106;
t38 = Ifges(5,2) * t72 + t109;
t37 = -Ifges(7,2) * t72 + t107;
t36 = -Ifges(6,3) * t72 + t105;
t35 = t95 * t72;
t32 = -mrSges(6,1) * t72 - mrSges(6,3) * t70;
t31 = t95 * t70;
t30 = -pkin(4) * t72 - t83;
t22 = -mrSges(7,2) * t73 + mrSges(7,3) * t101;
t21 = t114 * t72 + t83;
t19 = t80 * t71;
t14 = -Ifges(5,5) * t73 + (Ifges(5,1) * t72 - t109) * t71;
t13 = -Ifges(6,4) * t73 + (Ifges(6,1) * t72 + t105) * t71;
t12 = t103 + (Ifges(7,1) * t72 + t107) * t71;
t11 = -Ifges(5,6) * t73 + (-Ifges(5,2) * t70 + t108) * t71;
t10 = Ifges(7,6) * t73 + (Ifges(7,2) * t70 + t106) * t71;
t9 = -Ifges(6,6) * t73 + (Ifges(6,3) * t70 + t104) * t71;
t2 = t70 * t87 + t3;
t1 = pkin(5) * t73 + t28 + t63 + (-t20 - t87) * t72;
t15 = [0.2e1 * t1 * t24 + 0.2e1 * t8 * t17 + 0.2e1 * t5 * t18 + 0.2e1 * t2 * t22 + 0.2e1 * t7 * t23 + 0.2e1 * t6 * t25 + 0.2e1 * t4 * t26 + 0.2e1 * t3 * t27 + Ifges(2,3) + Ifges(3,3) + m(4) * (t51 * t67 + t53 ^ 2 + t42) + m(5) * (t6 ^ 2 + t7 ^ 2 + t42) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2 + t8 ^ 2) + m(3) * (t68 ^ 2 + t69 ^ 2) * pkin(1) ^ 2 + (-0.2e1 * t53 * mrSges(4,1) + (Ifges(4,2) + t85) * t73 + t84) * t73 + (0.2e1 * t53 * mrSges(4,2) + Ifges(4,1) * t71 + 0.2e1 * Ifges(4,4) * t73 + t19 * t121 + (t12 + t13 + t14 + t103) * t72 + (t73 * t96 + t10 - t11 + t9) * t70) * t71 + 0.2e1 * (mrSges(3,1) * t69 - mrSges(3,2) * t68) * pkin(1) + t89 * mrSges(4,3) * t121; (-t19 + t119 - t120) * t73 + ((t22 + t94) * t72 + (t24 + t93) * t70 + m(5) * (t81 - t102) + t78 + m(7) * (t1 * t70 + t2 * t72)) * t71; m(3) + m(4) * t89 + 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t65 * t90 + t67); -pkin(3) * t19 + t21 * t18 + t35 * t22 + t31 * t24 + t8 * t32 + t5 * t33 + m(7) * (t1 * t31 + t2 * t35 + t21 * t5) + (-t52 * mrSges(4,2) - t60 / 0.2e1 - t59 / 0.2e1 - t58 / 0.2e1 + Ifges(4,6)) * t73 + (t7 * mrSges(5,3) - t2 * mrSges(7,3) + t3 * mrSges(6,2) - t9 / 0.2e1 - t10 / 0.2e1 + t11 / 0.2e1 + (-Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t73) * t72 + (t103 / 0.2e1 + t12 / 0.2e1 + t13 / 0.2e1 + t14 / 0.2e1 + t4 * mrSges(6,2) - t6 * mrSges(5,3) - t1 * mrSges(7,3)) * t70 + (m(5) * t81 + t70 * t93 + t72 * t94 + t78) * pkin(8) + (Ifges(4,5) + (t39 / 0.2e1 + t40 / 0.2e1 + t41 / 0.2e1) * t72 + (t36 / 0.2e1 + t37 / 0.2e1 - t38 / 0.2e1) * t70 + (-m(5) * pkin(3) - t99) * t52) * t71 + t120 * t30; (-t32 + t33 + t99) * t73 + m(5) * (t92 + t112) + m(6) * (-t30 * t73 + t92) + t21 * t113 + (m(7) * (t31 * t70 + t35 * t72) - mrSges(4,2) + t90 * (mrSges(5,3) + t97)) * t71; -0.2e1 * pkin(3) * t34 + 0.2e1 * t21 * t33 + 0.2e1 * t30 * t32 + Ifges(4,3) + m(6) * (t30 ^ 2 + t91) + m(7) * (t21 ^ 2 + t31 ^ 2 + t35 ^ 2) + m(5) * (pkin(3) ^ 2 + t91) + (t118 * t35 - t36 - t37 + t38) * t72 + (t118 * t31 + t39 + t40 + t41) * t70 + 0.2e1 * (mrSges(6,2) + mrSges(5,3)) * pkin(8) * t90; t6 * mrSges(5,1) - t4 * mrSges(6,1) - t1 * mrSges(7,1) - t7 * mrSges(5,2) + t2 * mrSges(7,2) + t3 * mrSges(6,3) - pkin(4) * t26 - t114 * t24 + (t22 + t27) * qJ(5) + m(7) * (qJ(5) * t2 - t1 * t114) + m(6) * (-pkin(4) * t4 + qJ(5) * t3) - t85 * t73 + (-Ifges(7,5) * t72 - t70 * t96) * t71 - t84; t44 + ((-mrSges(5,2) + mrSges(6,3)) * t72 + (-mrSges(5,1) + t98) * t70) * t71 + m(6) * (-pkin(4) * t101 + t43) + m(7) * (-t101 * t114 + t43); -t31 * mrSges(7,1) + t35 * mrSges(7,2) + t59 + t58 + m(7) * (qJ(5) * t35 - t114 * t31) + t60 + (-pkin(4) * mrSges(6,2) + mrSges(7,3) * t114 - Ifges(7,5)) * t70 + (qJ(5) * t97 - Ifges(6,6) + Ifges(7,6)) * t72 + (m(6) * (t88 - t111) - t79 - t80) * pkin(8); 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * t114 * mrSges(7,1) + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * qJ(5) + m(6) * (pkin(4) ^ 2 + t75) + m(7) * (t114 ^ 2 + t75) + t85; m(6) * t4 + m(7) * t1 + t24 + t26; t115 * t101; m(7) * t31 + (m(6) * pkin(8) + t97) * t70; -m(6) * pkin(4) - m(7) * t114 + t98; t115; t119; t113; m(7) * t21 + t33; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
