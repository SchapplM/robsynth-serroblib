% Calculate joint inertia matrix for
% S6PRRPRP4
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:40:18
% EndTime: 2019-03-08 21:40:19
% DurationCPUTime: 0.75s
% Computational Cost: add. (608->219), mult. (1258->276), div. (0->0), fcn. (1056->8), ass. (0->90)
t116 = pkin(4) + pkin(8);
t115 = Ifges(6,3) + Ifges(7,3);
t66 = sin(qJ(3));
t69 = cos(qJ(3));
t114 = t66 ^ 2 + t69 ^ 2;
t71 = -pkin(3) - pkin(9);
t107 = m(6) * t71;
t113 = -mrSges(7,3) + t107;
t112 = m(6) + m(7);
t96 = mrSges(5,1) + mrSges(4,3);
t110 = m(5) * pkin(8) + mrSges(5,1);
t63 = sin(pkin(6));
t67 = sin(qJ(2));
t101 = t63 * t67;
t64 = cos(pkin(6));
t18 = t69 * t101 + t64 * t66;
t14 = t18 ^ 2;
t109 = -2 * mrSges(7,3);
t108 = m(7) * pkin(5);
t65 = sin(qJ(5));
t106 = Ifges(6,4) * t65;
t68 = cos(qJ(5));
t105 = Ifges(6,4) * t68;
t104 = Ifges(7,4) * t65;
t103 = Ifges(7,4) * t68;
t16 = t66 * t101 - t64 * t69;
t102 = t16 * t66;
t70 = cos(qJ(2));
t100 = t63 * t70;
t99 = t65 * t69;
t98 = t68 * mrSges(6,1);
t97 = t68 * t69;
t95 = -mrSges(6,2) - mrSges(7,2);
t94 = -Ifges(6,6) - Ifges(7,6);
t79 = -t66 * qJ(4) - pkin(2);
t22 = t71 * t69 + t79;
t40 = t116 * t66;
t6 = t68 * t22 + t65 * t40;
t25 = t66 * mrSges(7,1) + mrSges(7,3) * t99;
t26 = t66 * mrSges(6,1) + mrSges(6,3) * t99;
t93 = t25 + t26;
t27 = -t66 * mrSges(7,2) - mrSges(7,3) * t97;
t28 = -t66 * mrSges(6,2) - mrSges(6,3) * t97;
t92 = t27 + t28;
t35 = t65 * mrSges(6,1) + t68 * mrSges(6,2);
t91 = t35 + mrSges(5,3);
t34 = t65 * mrSges(7,1) + t68 * mrSges(7,2);
t90 = t115 * t66;
t89 = t114 * pkin(8) ^ 2;
t41 = t116 * t69;
t88 = t65 ^ 2 + t68 ^ 2;
t87 = qJ(6) * t69;
t86 = t18 * qJ(4);
t85 = -qJ(6) + t71;
t82 = -m(5) * pkin(3) + mrSges(5,2);
t81 = -mrSges(6,3) + t107;
t80 = t88 * mrSges(6,3);
t78 = mrSges(6,1) + mrSges(7,1) + t108;
t20 = mrSges(7,1) * t97 - mrSges(7,2) * t99;
t7 = t65 * t100 + t16 * t68;
t8 = -t68 * t100 + t16 * t65;
t76 = t65 * t8 + t68 * t7;
t75 = (t18 * t69 + t102) * pkin(8);
t74 = t94 * t68 + (-Ifges(6,5) - Ifges(7,5)) * t65;
t72 = qJ(4) ^ 2;
t57 = t63 ^ 2;
t52 = Ifges(6,5) * t68;
t51 = Ifges(7,5) * t68;
t46 = t65 * pkin(5) + qJ(4);
t44 = t57 * t70 ^ 2;
t39 = Ifges(6,1) * t68 - t106;
t38 = Ifges(7,1) * t68 - t104;
t37 = -Ifges(6,2) * t65 + t105;
t36 = -Ifges(7,2) * t65 + t103;
t33 = -t69 * mrSges(4,1) + t66 * mrSges(4,2);
t32 = t69 * mrSges(5,2) - t66 * mrSges(5,3);
t31 = -t69 * pkin(3) + t79;
t30 = t85 * t68;
t29 = t85 * t65;
t24 = t68 * t40;
t21 = (-t65 * mrSges(6,2) + t98) * t69;
t15 = pkin(5) * t97 + t41;
t12 = Ifges(6,5) * t66 + (-Ifges(6,1) * t65 - t105) * t69;
t11 = Ifges(7,5) * t66 + (-Ifges(7,1) * t65 - t103) * t69;
t10 = Ifges(6,6) * t66 + (-Ifges(6,2) * t68 - t106) * t69;
t9 = Ifges(7,6) * t66 + (-Ifges(7,2) * t68 - t104) * t69;
t5 = -t65 * t22 + t24;
t4 = -t68 * t87 + t6;
t3 = t66 * pkin(5) + t24 + (-t22 + t87) * t65;
t1 = [m(2) + m(3) * (t57 * t67 ^ 2 + t64 ^ 2 + t44) + (t7 ^ 2 + t8 ^ 2 + t14) * t112 + (m(5) + m(4)) * (t16 ^ 2 + t14 + t44); t92 * t8 + t93 * t7 + t96 * t102 + (-t67 * mrSges(3,2) + (mrSges(3,1) - t32 - t33) * t70) * t63 + (t96 * t69 + t20 + t21) * t18 + m(7) * (t15 * t18 + t3 * t7 + t4 * t8) + m(6) * (t41 * t18 + t5 * t7 + t6 * t8) + m(5) * (-t31 * t100 + t75) + m(4) * (pkin(2) * t100 + t75); -0.2e1 * pkin(2) * t33 + 0.2e1 * t15 * t20 + 0.2e1 * t41 * t21 + 0.2e1 * t3 * t25 + 0.2e1 * t5 * t26 + 0.2e1 * t4 * t27 + 0.2e1 * t6 * t28 + 0.2e1 * t31 * t32 + Ifges(3,3) + ((Ifges(4,1) + Ifges(5,2)) * t66 + t90) * t66 + m(5) * (t31 ^ 2 + t89) + m(4) * (pkin(2) ^ 2 + t89) + m(6) * (t41 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t15 ^ 2 + t3 ^ 2 + t4 ^ 2) + ((Ifges(5,3) + Ifges(4,2)) * t69 + (-t9 - t10) * t68 + (-t11 - t12) * t65 + ((2 * Ifges(4,4)) + (2 * Ifges(5,6)) + t74) * t66) * t69 + 0.2e1 * t96 * pkin(8) * t114; (-mrSges(4,1) + mrSges(5,2)) * t16 + (-mrSges(4,2) + t34 + t91) * t18 + m(7) * (t46 * t18 + t29 * t8 + t30 * t7) + m(6) * t86 + m(5) * (-pkin(3) * t16 + t86) + (-mrSges(6,3) + t113) * t76; t15 * t34 + t46 * t20 + t30 * t25 + t29 * t27 + t41 * t35 + m(7) * (t46 * t15 + t29 * t4 + t30 * t3) + (t71 * t26 - t3 * mrSges(7,3) + t11 / 0.2e1 + t12 / 0.2e1 + t81 * t5) * t68 + (t71 * t28 - t4 * mrSges(7,3) - t9 / 0.2e1 - t10 / 0.2e1 + t81 * t6) * t65 + (-pkin(3) * mrSges(5,1) + t51 / 0.2e1 + t52 / 0.2e1 + Ifges(4,5) - Ifges(5,4) + (-Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t65 + (-mrSges(4,1) + t82) * pkin(8)) * t66 + (-Ifges(5,5) + Ifges(4,6) + (-t36 / 0.2e1 - t37 / 0.2e1) * t68 + (-t38 / 0.2e1 - t39 / 0.2e1) * t65 + (mrSges(5,3) - mrSges(4,2)) * pkin(8)) * t69 + (m(6) * t41 + t110 * t69 + t21) * qJ(4); -0.2e1 * pkin(3) * mrSges(5,2) + 0.2e1 * t46 * t34 + Ifges(5,1) + Ifges(4,3) + (t30 * t109 + t38 + t39) * t68 + (t29 * t109 - t36 - t37) * t65 + m(7) * (t29 ^ 2 + t30 ^ 2 + t46 ^ 2) + m(6) * (t88 * t71 ^ 2 + t72) + m(5) * (pkin(3) ^ 2 + t72) + 0.2e1 * t91 * qJ(4) - 0.2e1 * t71 * t80; m(5) * t16 + t112 * t76; t93 * t68 + t110 * t66 + t92 * t65 + m(7) * (t68 * t3 + t65 * t4) + m(6) * (t68 * t5 + t65 * t6); -t80 + m(7) * (t65 * t29 + t68 * t30) + t82 + t113 * t88; t112 * t88 + m(5); t78 * t7 + t95 * t8; t5 * mrSges(6,1) + t3 * mrSges(7,1) - t6 * mrSges(6,2) - t4 * mrSges(7,2) + (m(7) * t3 + t25) * pkin(5) + t74 * t69 + t90; t71 * t98 + t30 * mrSges(7,1) - t29 * mrSges(7,2) + t51 + t52 + (m(7) * t30 - t68 * mrSges(7,3)) * pkin(5) + (-mrSges(6,2) * t71 + t94) * t65; t95 * t65 + t78 * t68; (0.2e1 * mrSges(7,1) + t108) * pkin(5) + t115; m(7) * t18; m(7) * t15 + t20; m(7) * t46 + t34; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
