% Calculate joint inertia matrix for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
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
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:03:32
% EndTime: 2019-03-09 10:03:34
% DurationCPUTime: 0.86s
% Computational Cost: add. (710->258), mult. (1252->313), div. (0->0), fcn. (890->4), ass. (0->93)
t118 = pkin(3) + pkin(7);
t117 = Ifges(6,2) + Ifges(5,3);
t68 = sin(qJ(2));
t70 = cos(qJ(2));
t116 = t68 ^ 2 + t70 ^ 2;
t67 = sin(qJ(4));
t69 = cos(qJ(4));
t93 = -t67 ^ 2 - t69 ^ 2;
t115 = 2 * mrSges(7,3);
t51 = qJ(5) * t67;
t77 = pkin(4) * t69 + t51;
t114 = m(6) * t77;
t72 = -pkin(2) - pkin(8);
t113 = t72 * t93;
t71 = -pkin(4) - pkin(5);
t112 = t69 * t71 - t51;
t111 = m(6) + m(7);
t109 = Ifges(5,4) * t67;
t108 = Ifges(5,4) * t69;
t107 = Ifges(7,4) * t67;
t106 = Ifges(7,4) * t69;
t105 = Ifges(6,5) * t67;
t104 = Ifges(6,5) * t69;
t103 = t67 * t70;
t102 = t69 * t70;
t101 = -mrSges(7,1) - mrSges(6,1);
t100 = -mrSges(6,2) + mrSges(7,3);
t99 = mrSges(7,2) + mrSges(6,3);
t98 = -Ifges(5,6) - Ifges(7,6);
t85 = -qJ(3) * t68 - pkin(1);
t21 = t70 * t72 + t85;
t43 = t118 * t68;
t6 = t69 * t21 + t67 * t43;
t25 = mrSges(5,1) * t68 + mrSges(5,3) * t103;
t90 = mrSges(6,2) * t103;
t26 = -t68 * mrSges(6,1) - t90;
t97 = t25 - t26;
t28 = -mrSges(5,2) * t68 - mrSges(5,3) * t102;
t29 = -mrSges(6,2) * t102 + mrSges(6,3) * t68;
t96 = t28 + t29;
t95 = t93 * t72 ^ 2;
t94 = t116 * pkin(7) ^ 2;
t44 = t118 * t70;
t92 = qJ(6) * t70;
t91 = qJ(6) + t72;
t89 = m(6) / 0.2e1 + m(5) / 0.2e1;
t3 = t68 * qJ(5) + t6;
t87 = Ifges(6,6) * t102 + t117 * t68;
t86 = -m(4) * pkin(2) + mrSges(4,2);
t35 = -t67 * mrSges(7,1) + t69 * mrSges(7,2);
t5 = -t67 * t21 + t43 * t69;
t84 = qJ(5) * t69 - qJ(3);
t4 = -pkin(4) * t68 - t5;
t81 = t3 * t67 - t4 * t69;
t80 = t5 * t69 + t6 * t67;
t79 = t69 * mrSges(5,1) - t67 * mrSges(5,2);
t78 = t69 * mrSges(6,1) + t67 * mrSges(6,3);
t19 = (-t69 * mrSges(7,1) - t67 * mrSges(7,2)) * t70;
t76 = t98 * t69 + (-Ifges(6,4) - Ifges(5,5) + Ifges(7,5)) * t67;
t74 = qJ(3) ^ 2;
t73 = qJ(5) ^ 2;
t57 = Ifges(6,4) * t69;
t56 = Ifges(5,5) * t69;
t54 = Ifges(6,6) * t67;
t46 = mrSges(7,3) * t103;
t42 = Ifges(5,1) * t69 - t109;
t41 = Ifges(6,1) * t69 + t105;
t40 = Ifges(7,1) * t69 + t107;
t39 = -Ifges(5,2) * t67 + t108;
t38 = Ifges(7,2) * t67 + t106;
t37 = Ifges(6,3) * t67 + t104;
t36 = mrSges(5,1) * t67 + mrSges(5,2) * t69;
t34 = mrSges(6,1) * t67 - mrSges(6,3) * t69;
t33 = -pkin(2) * t70 + t85;
t32 = pkin(4) * t67 - t84;
t31 = t91 * t69;
t30 = t91 * t67;
t27 = mrSges(7,2) * t68 + mrSges(7,3) * t102;
t24 = -t68 * mrSges(7,1) + t46;
t20 = t79 * t70;
t18 = t78 * t70;
t17 = t67 * t71 + t84;
t14 = Ifges(5,5) * t68 + (-Ifges(5,1) * t67 - t108) * t70;
t13 = Ifges(6,4) * t68 + (-Ifges(6,1) * t67 + t104) * t70;
t12 = -Ifges(7,5) * t68 + (-Ifges(7,1) * t67 + t106) * t70;
t11 = Ifges(5,6) * t68 + (-Ifges(5,2) * t69 - t109) * t70;
t10 = -Ifges(7,6) * t68 + (Ifges(7,2) * t69 - t107) * t70;
t9 = Ifges(6,6) * t68 + (Ifges(6,3) * t69 - t105) * t70;
t8 = t70 * t77 + t44;
t7 = t112 * t70 - t44;
t2 = t69 * t92 + t3;
t1 = t67 * t92 + t71 * t68 - t5;
t15 = [0.2e1 * t1 * t24 + 0.2e1 * t8 * t18 + 0.2e1 * t7 * t19 + 0.2e1 * t2 * t27 + 0.2e1 * t44 * t20 + 0.2e1 * t5 * t25 + 0.2e1 * t4 * t26 + 0.2e1 * t6 * t28 + 0.2e1 * t3 * t29 + Ifges(2,3) + m(4) * (t33 ^ 2 + t94) + m(3) * (pkin(1) ^ 2 + t94) + m(5) * (t44 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2 + t8 ^ 2) + (-0.2e1 * pkin(1) * mrSges(3,2) - 0.2e1 * t33 * mrSges(4,3) + (Ifges(4,2) + Ifges(3,1) + Ifges(7,3)) * t68 + t87) * t68 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * t33 * mrSges(4,2) + (Ifges(4,3) + Ifges(3,2)) * t70 + (t10 - t11 + t9) * t69 + (-t12 - t13 - t14) * t67 + ((2 * Ifges(3,4)) + (2 * Ifges(4,6)) + t76) * t68) * t70 + 0.2e1 * (mrSges(4,1) + mrSges(3,3)) * pkin(7) * t116; qJ(3) * t20 + t17 * t19 + t32 * t18 - t31 * t24 + t30 * t27 + t8 * t34 + t7 * t35 + t44 * t36 + (t12 / 0.2e1 + t13 / 0.2e1 + t14 / 0.2e1 + t4 * mrSges(6,2) - t5 * mrSges(5,3) - t1 * mrSges(7,3) + t97 * t72) * t69 + (-t3 * mrSges(6,2) - t6 * mrSges(5,3) + t2 * mrSges(7,3) + t9 / 0.2e1 + t10 / 0.2e1 - t11 / 0.2e1 + t96 * t72) * t67 + m(6) * (t32 * t8 + t72 * t81) + m(5) * (qJ(3) * t44 + t72 * t80) + m(7) * (-t1 * t31 + t17 * t7 + t2 * t30) + (-Ifges(7,5) * t69 / 0.2e1 - pkin(2) * mrSges(4,1) + t56 / 0.2e1 + t57 / 0.2e1 + t54 / 0.2e1 - Ifges(4,4) + Ifges(3,5) + (-Ifges(7,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * t67 + (-mrSges(3,1) + t86) * pkin(7)) * t68 + (qJ(3) * mrSges(4,1) - Ifges(4,5) + Ifges(3,6) + (t37 / 0.2e1 + t38 / 0.2e1 - t39 / 0.2e1) * t69 + (-t40 / 0.2e1 - t41 / 0.2e1 - t42 / 0.2e1) * t67 + (m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * pkin(7)) * t70; -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t17 * t35 + 0.2e1 * t32 * t34 + Ifges(4,1) + Ifges(3,3) + 0.2e1 * (t36 + mrSges(4,3)) * qJ(3) + (t31 * t115 + t40 + t41 + t42) * t69 + (t30 * t115 + t37 + t38 - t39) * t67 + m(6) * (t32 ^ 2 - t95) + m(7) * (t17 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(5) * (t74 - t95) + m(4) * (pkin(2) ^ 2 + t74) + 0.2e1 * (mrSges(6,2) + mrSges(5,3)) * t113; (m(4) * pkin(7) + mrSges(4,1)) * t68 + (-t24 + t97) * t69 + (t27 + t96) * t67 + m(7) * (-t1 * t69 + t2 * t67) + m(6) * t81 + m(5) * t80; m(7) * (t30 * t67 + t31 * t69) - 0.2e1 * t89 * t113 + t86 - t93 * (-mrSges(5,3) + t100); m(4) - 0.2e1 * (m(7) / 0.2e1 + t89) * t93; t5 * mrSges(5,1) - t4 * mrSges(6,1) - t1 * mrSges(7,1) - t6 * mrSges(5,2) + t2 * mrSges(7,2) + t3 * mrSges(6,3) + Ifges(7,3) * t68 - pkin(4) * t26 + t71 * t24 + (t27 + t29) * qJ(5) + m(7) * (qJ(5) * t2 + t1 * t71) + m(6) * (-pkin(4) * t4 + qJ(5) * t3) + t76 * t70 + t87; t57 + t54 + m(7) * (qJ(5) * t30 - t31 * t71) + t31 * mrSges(7,1) + t30 * mrSges(7,2) + t56 + (-pkin(4) * mrSges(6,2) - t71 * mrSges(7,3) - Ifges(7,5)) * t69 + (qJ(5) * t100 + t98) * t67 + (t78 + t79 + t114) * t72; (mrSges(5,1) - t101) * t69 + t114 - m(7) * t112 + (-mrSges(5,2) + t99) * t67; 0.2e1 * pkin(4) * mrSges(6,1) - 0.2e1 * t71 * mrSges(7,1) + Ifges(7,3) + 0.2e1 * t99 * qJ(5) + m(7) * (t71 ^ 2 + t73) + m(6) * (pkin(4) ^ 2 + t73) + t117; m(6) * t4 + m(7) * t1 + t101 * t68 + t46 - t90; -m(7) * t31 + (-m(6) * t72 - t100) * t69; -t111 * t69; -m(6) * pkin(4) + m(7) * t71 + t101; t111; m(7) * t7 + t19; m(7) * t17 + t35; 0; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
