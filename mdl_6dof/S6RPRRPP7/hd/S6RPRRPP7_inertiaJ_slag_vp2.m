% Calculate joint inertia matrix for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2018-11-23 16:15
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:14:51
% EndTime: 2018-11-23 16:14:52
% DurationCPUTime: 0.91s
% Computational Cost: add. (651->263), mult. (1229->322), div. (0->0), fcn. (868->4), ass. (0->99)
t67 = sin(qJ(4));
t69 = cos(qJ(4));
t90 = t67 ^ 2 + t69 ^ 2;
t120 = pkin(8) * t90;
t119 = Ifges(6,2) + Ifges(5,3);
t72 = (-pkin(1) - pkin(7));
t118 = -2 * t72;
t117 = -2 * mrSges(7,3);
t116 = 2 * qJ(2);
t115 = m(6) + m(7);
t114 = pkin(4) + pkin(5);
t70 = cos(qJ(3));
t113 = m(7) * t70;
t111 = Ifges(5,4) * t67;
t110 = Ifges(5,4) * t69;
t109 = Ifges(7,4) * t67;
t108 = Ifges(7,4) * t69;
t107 = Ifges(6,5) * t67;
t106 = Ifges(6,5) * t69;
t68 = sin(qJ(3));
t105 = Ifges(7,5) * t68;
t104 = t67 * t68;
t103 = t67 * t70;
t102 = t68 * t72;
t101 = t69 * t70;
t32 = -mrSges(5,1) * t69 + mrSges(5,2) * t67;
t100 = mrSges(4,1) - t32;
t99 = -mrSges(6,1) - mrSges(7,1);
t98 = mrSges(6,2) - mrSges(7,3);
t97 = mrSges(6,3) + mrSges(7,2);
t96 = -Ifges(5,6) - Ifges(7,6);
t95 = pkin(8) - qJ(6);
t22 = -mrSges(5,2) * t68 - mrSges(5,3) * t103;
t26 = -mrSges(6,2) * t103 + mrSges(6,3) * t68;
t94 = t22 + t26;
t24 = mrSges(5,1) * t68 - mrSges(5,3) * t101;
t44 = mrSges(6,2) * t101;
t25 = -t68 * mrSges(6,1) + t44;
t93 = -t24 + t25;
t27 = pkin(3) * t68 - pkin(8) * t70 + qJ(2);
t7 = t69 * t102 + t67 * t27;
t92 = t68 * t120;
t31 = t69 * mrSges(7,1) + t67 * mrSges(7,2);
t91 = t90 * pkin(8) ^ 2;
t63 = t68 ^ 2;
t65 = t70 ^ 2;
t89 = t65 + t63;
t88 = qJ(5) * t69;
t87 = qJ(6) * t70;
t52 = t68 * qJ(5);
t86 = mrSges(7,3) * t101;
t4 = t52 + t7;
t84 = t89 * mrSges(4,3);
t83 = t67 * qJ(5) + pkin(3);
t40 = t67 * t102;
t6 = t69 * t27 - t40;
t17 = -mrSges(7,1) * t103 + mrSges(7,2) * t101;
t81 = Ifges(6,6) * t103 + t119 * t68 + (Ifges(6,4) + Ifges(5,5)) * t101;
t5 = -t68 * pkin(4) - t6;
t80 = t4 * t69 + t5 * t67;
t79 = -t6 * t67 + t69 * t7;
t78 = t67 * mrSges(5,1) + t69 * mrSges(5,2);
t77 = t67 * mrSges(6,1) - t69 * mrSges(6,3);
t76 = -pkin(4) * t67 + t88;
t74 = qJ(2) ^ 2;
t73 = qJ(5) ^ 2;
t66 = t72 ^ 2;
t59 = Ifges(6,4) * t67;
t58 = Ifges(5,5) * t67;
t56 = Ifges(5,6) * t69;
t51 = t65 * t72;
t50 = t65 * t66;
t42 = t69 * t52;
t39 = Ifges(5,1) * t67 + t110;
t38 = Ifges(6,1) * t67 - t106;
t37 = Ifges(7,1) * t67 - t108;
t36 = Ifges(5,2) * t69 + t111;
t35 = -Ifges(7,2) * t69 + t109;
t34 = -Ifges(6,3) * t69 + t107;
t33 = t95 * t69;
t30 = -mrSges(6,1) * t69 - mrSges(6,3) * t67;
t29 = t95 * t67;
t28 = -pkin(4) * t69 - t83;
t23 = -t68 * mrSges(7,1) - t86;
t21 = mrSges(7,2) * t68 + mrSges(7,3) * t103;
t19 = t114 * t69 + t83;
t18 = t78 * t70;
t16 = t77 * t70;
t14 = Ifges(5,5) * t68 + (Ifges(5,1) * t69 - t111) * t70;
t13 = Ifges(6,4) * t68 + (Ifges(6,1) * t69 + t107) * t70;
t12 = -t105 + (Ifges(7,1) * t69 + t109) * t70;
t11 = Ifges(5,6) * t68 + (-Ifges(5,2) * t67 + t110) * t70;
t10 = -Ifges(7,6) * t68 + (Ifges(7,2) * t67 + t108) * t70;
t9 = Ifges(6,6) * t68 + (Ifges(6,3) * t67 + t106) * t70;
t8 = (-t72 - t76) * t70;
t3 = (-t114 * t67 + t72 + t88) * t70;
t2 = t67 * t87 + t4;
t1 = t40 + (-t27 - t87) * t69 - t114 * t68;
t15 = [-(2 * pkin(1) * mrSges(3,2)) + (mrSges(3,3) * t116) + 0.2e1 * t1 * t23 + 0.2e1 * t8 * t16 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t21 + 0.2e1 * t7 * t22 + 0.2e1 * t6 * t24 + 0.2e1 * t5 * t25 + 0.2e1 * t4 * t26 + Ifges(3,1) + Ifges(2,3) + t84 * t118 + (mrSges(4,1) * t116 + (Ifges(7,3) + Ifges(4,2)) * t68 + t81) * t68 + m(4) * (t63 * t66 + t50 + t74) + (m(3) * (pkin(1) ^ 2 + t74)) + m(5) * (t6 ^ 2 + t7 ^ 2 + t50) + m(6) * (t4 ^ 2 + t5 ^ 2 + t8 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + ((mrSges(4,2) * t116) + Ifges(4,1) * t70 - 0.2e1 * Ifges(4,4) * t68 + t18 * t118 + (t12 + t13 + t14 - t105) * t69 + (t96 * t68 + t10 - t11 + t9) * t67) * t70; -(m(3) * pkin(1)) + mrSges(3,2) - t84 + (-t16 + t17 - t18) * t70 + ((t21 + t94) * t69 + (t23 + t93) * t67) * t68 + m(6) * (t80 * t68 - t70 * t8) + m(7) * (t70 * t3 + (t1 * t67 + t2 * t69) * t68) + m(5) * (t79 * t68 + t51) + m(4) * (t63 * t72 + t51); m(3) + m(4) * t89 + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1 + m(5) / 0.2e1) * (t90 * t63 + t65); -pkin(3) * t18 + t19 * t17 + t33 * t21 + t29 * t23 + t3 * t31 + t8 * t30 + m(7) * (t1 * t29 + t19 * t3 + t2 * t33) + (-(t72 * mrSges(4,2)) + t59 / 0.2e1 + t58 / 0.2e1 + t56 / 0.2e1 - Ifges(4,6)) * t68 + (t7 * mrSges(5,3) - t2 * mrSges(7,3) + t4 * mrSges(6,2) - t9 / 0.2e1 - t10 / 0.2e1 + t11 / 0.2e1 + (Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t68) * t69 + (t5 * mrSges(6,2) - t6 * mrSges(5,3) - t1 * mrSges(7,3) - t105 / 0.2e1 + t12 / 0.2e1 + t13 / 0.2e1 + t14 / 0.2e1) * t67 + (m(5) * t79 + m(6) * t80 + t93 * t67 + t94 * t69) * pkin(8) + (Ifges(4,5) + (m(5) * pkin(3) + t100) * t72 + (t37 / 0.2e1 + t38 / 0.2e1 + t39 / 0.2e1) * t69 + (t34 / 0.2e1 + t35 / 0.2e1 - t36 / 0.2e1) * t67) * t70 + (m(6) * t8 + t16) * t28; (-t30 + t31 + t100) * t70 + m(6) * (-t28 * t70 + t92) + t19 * t113 + m(5) * (pkin(3) * t70 + t92) + (m(7) * (t29 * t67 + t33 * t69) - mrSges(4,2) + t90 * (mrSges(5,3) + t98)) * t68; -0.2e1 * pkin(3) * t32 + 0.2e1 * t19 * t31 + 0.2e1 * t28 * t30 + Ifges(4,3) + m(6) * (t28 ^ 2 + t91) + m(7) * (t19 ^ 2 + t29 ^ 2 + t33 ^ 2) + m(5) * (pkin(3) ^ 2 + t91) + (t33 * t117 - t34 - t35 + t36) * t69 + (t29 * t117 + t37 + t38 + t39) * t67 + 0.2e1 * (mrSges(6,2) + mrSges(5,3)) * t120; t6 * mrSges(5,1) - t5 * mrSges(6,1) - t1 * mrSges(7,1) - t7 * mrSges(5,2) + t2 * mrSges(7,2) + t4 * mrSges(6,3) + Ifges(7,3) * t68 - pkin(4) * t25 - t114 * t23 + (t21 + t26) * qJ(5) + m(6) * (-pkin(4) * t5 + qJ(5) * t4) + m(7) * (qJ(5) * t2 - t1 * t114) + (-Ifges(7,5) * t69 + t96 * t67) * t70 + t81; (-mrSges(5,1) + t99) * t104 + m(6) * (-pkin(4) * t104 + t42) + m(7) * (-t104 * t114 + t42) + (-mrSges(5,2) + t97) * t68 * t69; m(7) * (qJ(5) * t33 - t114 * t29) + t33 * mrSges(7,2) - t29 * mrSges(7,1) + t59 + t58 + t56 + (-mrSges(6,2) * pkin(4) + mrSges(7,3) * t114 - Ifges(7,5)) * t67 + (t98 * qJ(5) - Ifges(6,6) + Ifges(7,6)) * t69 + (m(6) * t76 - t77 - t78) * pkin(8); 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * t114 * mrSges(7,1) + Ifges(7,3) + 0.2e1 * t97 * qJ(5) + m(6) * (pkin(4) ^ 2 + t73) + m(7) * (t114 ^ 2 + t73) + t119; m(6) * t5 + m(7) * t1 + t99 * t68 + t44 - t86; t115 * t104; m(7) * t29 + (m(6) * pkin(8) + t98) * t67; -m(6) * pkin(4) - m(7) * t114 + t99; t115; m(7) * t3 + t17; t113; m(7) * t19 + t31; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
