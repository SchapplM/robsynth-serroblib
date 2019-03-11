% Calculate joint inertia matrix for
% S6RPRRPP3
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:35:09
% EndTime: 2019-03-09 04:35:11
% DurationCPUTime: 0.91s
% Computational Cost: add. (688->253), mult. (1329->316), div. (0->0), fcn. (977->6), ass. (0->97)
t119 = -Ifges(7,4) - Ifges(6,5);
t71 = sin(qJ(4));
t73 = cos(qJ(4));
t90 = t71 ^ 2 + t73 ^ 2;
t68 = sin(pkin(9));
t55 = pkin(1) * t68 + pkin(7);
t118 = 0.2e1 * t55;
t72 = sin(qJ(3));
t79 = -t71 * mrSges(6,2) - t73 * mrSges(6,3);
t19 = t79 * t72;
t101 = t71 * t72;
t88 = qJ(5) * t73;
t44 = t72 * t88;
t92 = -pkin(4) * t101 + t44;
t8 = t72 * t55 - t92;
t117 = m(6) * t8 + t19;
t116 = 2 * mrSges(7,1);
t114 = m(6) + m(7);
t113 = pkin(5) + pkin(8);
t74 = cos(qJ(3));
t112 = pkin(3) * t74;
t111 = pkin(8) * t72;
t110 = Ifges(5,4) * t71;
t109 = Ifges(5,4) * t73;
t108 = Ifges(6,4) * t74;
t107 = Ifges(5,6) * t74;
t106 = Ifges(6,6) * t71;
t105 = Ifges(6,6) * t73;
t104 = Ifges(7,6) * t71;
t103 = Ifges(7,6) * t73;
t102 = t55 * t74;
t100 = t72 * t73;
t99 = t73 * mrSges(7,1);
t32 = -mrSges(5,1) * t73 + mrSges(5,2) * t71;
t98 = mrSges(4,1) - t32;
t97 = mrSges(6,1) + mrSges(7,1);
t96 = mrSges(6,3) + mrSges(7,2);
t70 = -pkin(4) - qJ(6);
t23 = -mrSges(5,1) * t74 - mrSges(5,3) * t100;
t27 = mrSges(6,1) * t100 - t74 * mrSges(6,2);
t95 = -t23 + t27;
t25 = mrSges(6,1) * t101 + mrSges(6,3) * t74;
t26 = -mrSges(7,1) * t101 - t74 * mrSges(7,2);
t94 = -t25 + t26;
t69 = cos(pkin(9));
t56 = -pkin(1) * t69 - pkin(2);
t21 = -t111 + t56 - t112;
t7 = t73 * t102 + t71 * t21;
t93 = t90 * t111;
t24 = t74 * mrSges(7,3) + t72 * t99;
t91 = t90 * pkin(8) ^ 2;
t65 = t72 ^ 2;
t67 = t74 ^ 2;
t89 = t65 + t67;
t86 = Ifges(7,1) + Ifges(6,1) + Ifges(5,3);
t85 = -qJ(5) * t71 - pkin(3);
t28 = t71 * t102;
t6 = t21 * t73 - t28;
t84 = t119 * t101 + (-Ifges(5,5) - Ifges(7,5)) * t100;
t82 = m(7) * t70 + mrSges(6,2) - mrSges(7,3);
t3 = qJ(5) * t74 - t7;
t81 = -t6 * t71 + t7 * t73;
t80 = t71 * mrSges(5,1) + t73 * mrSges(5,2);
t63 = t74 * pkin(4);
t4 = -t6 + t63;
t78 = m(6) * (-t3 * t73 + t4 * t71);
t75 = qJ(5) ^ 2;
t60 = Ifges(5,5) * t71;
t59 = Ifges(7,5) * t71;
t58 = Ifges(5,6) * t73;
t54 = t55 ^ 2;
t43 = t113 * t73;
t42 = t113 * t71;
t40 = t65 * t54;
t39 = Ifges(5,1) * t71 + t109;
t38 = Ifges(5,2) * t73 + t110;
t37 = -Ifges(6,2) * t71 - t105;
t36 = -Ifges(7,2) * t73 + t104;
t35 = -Ifges(6,3) * t73 - t106;
t34 = Ifges(7,3) * t71 - t103;
t33 = -mrSges(7,2) * t71 - mrSges(7,3) * t73;
t31 = mrSges(6,2) * t73 - mrSges(6,3) * t71;
t30 = -pkin(4) * t73 + t85;
t22 = mrSges(5,2) * t74 - mrSges(5,3) * t101;
t20 = t80 * t72;
t18 = (-t73 * mrSges(7,2) + t71 * mrSges(7,3)) * t72;
t17 = t70 * t73 + t85;
t14 = -t108 + (-Ifges(6,2) * t73 + t106) * t72;
t13 = -Ifges(7,4) * t74 + (Ifges(7,2) * t71 + t103) * t72;
t12 = -Ifges(6,5) * t74 + (Ifges(6,3) * t71 - t105) * t72;
t11 = -Ifges(7,5) * t74 + (Ifges(7,3) * t73 + t104) * t72;
t10 = -Ifges(5,5) * t74 + (Ifges(5,1) * t73 - t110) * t72;
t9 = -t107 + (-Ifges(5,2) * t71 + t109) * t72;
t5 = qJ(6) * t101 + t8;
t2 = -pkin(5) * t101 - t3;
t1 = qJ(6) * t74 + t28 + t63 + (pkin(5) * t72 - t21) * t73;
t15 = [0.2e1 * t1 * t24 + 0.2e1 * t5 * t18 + 0.2e1 * t8 * t19 + 0.2e1 * t2 * t26 + 0.2e1 * t7 * t22 + 0.2e1 * t6 * t23 + 0.2e1 * t3 * t25 + 0.2e1 * t4 * t27 + Ifges(2,3) + Ifges(3,3) + m(4) * (t54 * t67 + t56 ^ 2 + t40) + m(5) * (t6 ^ 2 + t7 ^ 2 + t40) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2 + t8 ^ 2) + m(3) * (t68 ^ 2 + t69 ^ 2) * pkin(1) ^ 2 + (-0.2e1 * t56 * mrSges(4,1) + (Ifges(4,2) + t86) * t74 + t84) * t74 + (0.2e1 * t56 * mrSges(4,2) + Ifges(4,1) * t72 + 0.2e1 * Ifges(4,4) * t74 + t20 * t118 + (t10 + t11 - t14 + t108) * t73 + (t12 + t13 - t9 + t107) * t71) * t72 + 0.2e1 * (mrSges(3,1) * t69 - mrSges(3,2) * t68) * pkin(1) + t89 * mrSges(4,3) * t118; (-m(7) * t5 - t117 - t18 - t20) * t74 + ((t22 + t94) * t73 + (t24 + t95) * t71 + t78 + m(7) * (t1 * t71 + t2 * t73) + m(5) * (t81 - t102)) * t72; m(3) + m(4) * t89 + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1 + m(5) / 0.2e1) * (t65 * t90 + t67); -pkin(3) * t20 + t17 * t18 + t42 * t24 + t43 * t26 + t8 * t31 + t5 * t33 + m(7) * (t1 * t42 + t17 * t5 + t2 * t43) + (-t55 * mrSges(4,2) - t59 / 0.2e1 - t60 / 0.2e1 - t58 / 0.2e1 + Ifges(4,6)) * t74 + (t7 * mrSges(5,3) + t2 * mrSges(7,1) - t3 * mrSges(6,1) + t9 / 0.2e1 - t12 / 0.2e1 - t13 / 0.2e1 + (Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1) * t74) * t73 + (t108 / 0.2e1 + t11 / 0.2e1 - t14 / 0.2e1 + t10 / 0.2e1 + t4 * mrSges(6,1) - t6 * mrSges(5,3) + t1 * mrSges(7,1)) * t71 + ((t22 - t25) * t73 + t95 * t71 + m(5) * t81 + t78) * pkin(8) + (Ifges(4,5) + (t34 / 0.2e1 - t37 / 0.2e1 + t39 / 0.2e1) * t73 + (t35 / 0.2e1 + t36 / 0.2e1 - t38 / 0.2e1) * t71 + (-m(5) * pkin(3) - t98) * t55) * t72 + t117 * t30; m(6) * t93 + m(5) * (t93 + t112) + (m(7) * (t42 * t71 + t43 * t73) - mrSges(4,2) + t90 * (mrSges(5,3) + t97)) * t72 + (-m(6) * t30 - m(7) * t17 - t31 - t33 + t98) * t74; -0.2e1 * pkin(3) * t32 + 0.2e1 * t17 * t33 + 0.2e1 * t30 * t31 + Ifges(4,3) + m(6) * (t30 ^ 2 + t91) + m(7) * (t17 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (pkin(3) ^ 2 + t91) + (t116 * t43 - t35 - t36 + t38) * t73 + (t116 * t42 + t34 - t37 + t39) * t71 + 0.2e1 * (mrSges(6,1) + mrSges(5,3)) * pkin(8) * t90; t6 * mrSges(5,1) - t7 * mrSges(5,2) + t4 * mrSges(6,2) + t2 * mrSges(7,2) - t3 * mrSges(6,3) - t1 * mrSges(7,3) - pkin(4) * t27 + t70 * t24 + (-Ifges(6,4) * t73 - Ifges(5,6) * t71) * t72 + t94 * qJ(5) + m(6) * (-pkin(4) * t4 - qJ(5) * t3) + m(7) * (qJ(5) * t2 + t1 * t70) - t86 * t74 - t84; m(6) * t92 + m(7) * t44 + ((-mrSges(5,2) + t96) * t73 + (-mrSges(5,1) + t82) * t71) * t72; t59 + m(7) * (qJ(5) * t43 + t42 * t70) + t58 + t60 + t43 * mrSges(7,2) - t42 * mrSges(7,3) + (-pkin(4) * mrSges(6,1) + t70 * mrSges(7,1) - Ifges(6,4)) * t71 + (qJ(5) * t97 + t119) * t73 + (m(6) * (-pkin(4) * t71 + t88) - t79 - t80) * pkin(8); -0.2e1 * pkin(4) * mrSges(6,2) - 0.2e1 * t70 * mrSges(7,3) + 0.2e1 * t96 * qJ(5) + m(6) * (pkin(4) ^ 2 + t75) + m(7) * (t70 ^ 2 + t75) + t86; m(6) * t4 + m(7) * t1 + t24 + t27; t114 * t101; m(7) * t42 + (m(6) * pkin(8) + t97) * t71; -m(6) * pkin(4) + t82; t114; m(7) * t2 + t26; m(7) * t100; m(7) * t43 + t99; m(7) * qJ(5) + mrSges(7,2); 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
