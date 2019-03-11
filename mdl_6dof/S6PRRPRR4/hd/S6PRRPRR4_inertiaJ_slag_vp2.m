% Calculate joint inertia matrix for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:10:21
% EndTime: 2019-03-08 22:10:23
% DurationCPUTime: 0.84s
% Computational Cost: add. (929->227), mult. (1890->305), div. (0->0), fcn. (1835->10), ass. (0->91)
t68 = cos(qJ(3));
t59 = t68 ^ 2;
t64 = sin(qJ(3));
t121 = t64 ^ 2 + t59;
t61 = cos(pkin(6));
t60 = sin(pkin(6));
t65 = sin(qJ(2));
t94 = t60 * t65;
t22 = -t61 * t68 + t64 * t94;
t24 = t61 * t64 + t68 * t94;
t120 = t22 * t64 + t24 * t68;
t62 = sin(qJ(6));
t66 = cos(qJ(6));
t87 = t62 ^ 2 + t66 ^ 2;
t84 = t87 * mrSges(7,3);
t36 = -t66 * mrSges(7,1) + mrSges(7,2) * t62;
t90 = mrSges(6,1) - t36;
t119 = -mrSges(4,2) + mrSges(5,3);
t118 = mrSges(4,3) + mrSges(5,2);
t63 = sin(qJ(5));
t67 = cos(qJ(5));
t27 = t63 * t64 + t67 * t68;
t28 = -t63 * t68 + t64 * t67;
t98 = t28 * t62;
t12 = -mrSges(7,2) * t27 - mrSges(7,3) * t98;
t97 = t28 * t66;
t13 = mrSges(7,1) * t27 - mrSges(7,3) * t97;
t117 = t66 * t12 - t62 * t13;
t116 = -m(5) * pkin(3) - mrSges(5,1);
t115 = m(7) * pkin(10) + mrSges(7,3);
t114 = m(7) * pkin(5) + t90;
t8 = -t22 * t67 + t24 * t63;
t113 = t8 ^ 2;
t108 = pkin(8) - pkin(9);
t42 = t108 * t68;
t85 = t108 * t64;
t15 = t42 * t63 - t67 * t85;
t112 = t15 ^ 2;
t111 = 0.2e1 * t15;
t110 = -0.2e1 * t36;
t109 = t66 / 0.2e1;
t107 = t15 * t8;
t106 = t67 * t8;
t105 = Ifges(7,4) * t62;
t104 = Ifges(7,4) * t66;
t10 = t22 * t63 + t24 * t67;
t103 = t10 * mrSges(6,2);
t102 = t15 * t67;
t17 = t67 * t42 + t63 * t85;
t101 = t17 * mrSges(6,2);
t70 = -pkin(3) - pkin(4);
t33 = -qJ(4) * t63 + t67 * t70;
t96 = t33 * mrSges(6,1);
t34 = t67 * qJ(4) + t63 * t70;
t95 = t34 * mrSges(6,2);
t69 = cos(qJ(2));
t93 = t60 * t69;
t89 = Ifges(7,5) * t97 + Ifges(7,3) * t27;
t88 = t121 * pkin(8) ^ 2;
t39 = Ifges(7,5) * t62 + Ifges(7,6) * t66;
t86 = t39 / 0.2e1 - Ifges(6,6);
t35 = -t68 * pkin(3) - t64 * qJ(4) - pkin(2);
t32 = -pkin(10) + t34;
t83 = t87 * t32;
t82 = t87 * t63;
t25 = t68 * pkin(4) - t35;
t79 = t120 * pkin(8);
t7 = pkin(5) * t27 - pkin(10) * t28 + t25;
t1 = -t17 * t62 + t66 * t7;
t2 = t17 * t66 + t62 * t7;
t78 = -t1 * t62 + t2 * t66;
t3 = -t10 * t62 + t66 * t93;
t4 = t10 * t66 + t62 * t93;
t77 = t3 * t62 - t4 * t66;
t76 = mrSges(7,1) * t62 + mrSges(7,2) * t66;
t40 = Ifges(7,2) * t66 + t105;
t41 = Ifges(7,1) * t62 + t104;
t74 = t66 * t40 + t62 * t41 + Ifges(6,3);
t73 = -t62 * t40 / 0.2e1 + t41 * t109 + Ifges(6,5);
t58 = t67 ^ 2;
t55 = t63 ^ 2;
t53 = t60 ^ 2;
t43 = t53 * t69 ^ 2;
t38 = -mrSges(4,1) * t68 + mrSges(4,2) * t64;
t37 = -mrSges(5,1) * t68 - mrSges(5,3) * t64;
t31 = pkin(5) - t33;
t14 = mrSges(6,1) * t27 + mrSges(6,2) * t28;
t11 = t76 * t28;
t6 = Ifges(7,5) * t27 + (Ifges(7,1) * t66 - t105) * t28;
t5 = Ifges(7,6) * t27 + (-Ifges(7,2) * t62 + t104) * t28;
t9 = [m(2) + m(7) * (t3 ^ 2 + t4 ^ 2 + t113) + m(6) * (t10 ^ 2 + t113 + t43) + m(3) * (t53 * t65 ^ 2 + t61 ^ 2 + t43) + (m(5) + m(4)) * (t22 ^ 2 + t24 ^ 2 + t43); t8 * t11 + t4 * t12 + t3 * t13 + (-t10 * t27 + t28 * t8) * mrSges(6,3) + (-t65 * mrSges(3,2) + (mrSges(3,1) + t14 - t37 - t38) * t69) * t60 + m(7) * (t1 * t3 + t2 * t4 + t107) + m(6) * (t10 * t17 + t25 * t93 + t107) + m(5) * (-t35 * t93 + t79) + m(4) * (pkin(2) * t93 + t79) + t118 * t120; -0.2e1 * pkin(2) * t38 + 0.2e1 * t1 * t13 + t11 * t111 + 0.2e1 * t2 * t12 + 0.2e1 * t25 * t14 + 0.2e1 * t35 * t37 + Ifges(3,3) + (Ifges(5,3) + Ifges(4,2)) * t59 + (-0.2e1 * mrSges(6,3) * t17 + Ifges(6,2) * t27 + t89) * t27 + ((Ifges(4,1) + Ifges(5,1)) * t64 + 0.2e1 * (Ifges(4,4) - Ifges(5,5)) * t68) * t64 + (mrSges(6,3) * t111 + Ifges(6,1) * t28 - t62 * t5 + t66 * t6 + (-Ifges(7,6) * t62 - (2 * Ifges(6,4))) * t27) * t28 + m(5) * (t35 ^ 2 + t88) + m(4) * (pkin(2) ^ 2 + t88) + m(6) * (t17 ^ 2 + t25 ^ 2 + t112) + m(7) * (t1 ^ 2 + t2 ^ 2 + t112) + 0.2e1 * t118 * pkin(8) * t121; t103 + t90 * t8 + t119 * t24 + (-mrSges(4,1) - mrSges(5,1)) * t22 + t77 * mrSges(7,3) + m(7) * (t31 * t8 - t32 * t77) + m(6) * (t10 * t34 - t33 * t8) + m(5) * (-pkin(3) * t22 + qJ(4) * t24); t101 + t31 * t11 + t90 * t15 + (t32 * t12 - t2 * mrSges(7,3) - t5 / 0.2e1) * t66 + (t1 * mrSges(7,3) - t32 * t13 - t6 / 0.2e1) * t62 + m(7) * (t15 * t31 + t32 * t78) + m(6) * (-t15 * t33 + t17 * t34) + (qJ(4) * mrSges(5,2) + Ifges(4,6) - Ifges(5,6)) * t68 + (-pkin(3) * mrSges(5,2) + Ifges(5,4) + Ifges(4,5)) * t64 + (-t34 * mrSges(6,3) - t86) * t27 + (-t33 * mrSges(6,3) - t73) * t28 + ((m(5) * qJ(4) + t119) * t68 + (-mrSges(4,1) + t116) * t64) * pkin(8); 0.2e1 * pkin(3) * mrSges(5,1) - 0.2e1 * t96 + 0.2e1 * t95 + 0.2e1 * qJ(4) * mrSges(5,3) + t31 * t110 + Ifges(5,2) + Ifges(4,3) - 0.2e1 * t32 * t84 + m(7) * (t32 ^ 2 * t87 + t31 ^ 2) + m(6) * (t33 ^ 2 + t34 ^ 2) + m(5) * (pkin(3) ^ 2 + qJ(4) ^ 2) + t74; m(7) * (-t63 * t77 - t106) + m(6) * (t10 * t63 - t106) + m(5) * t22; (-t28 * mrSges(6,3) - t11) * t67 + (m(5) * pkin(8) + mrSges(5,2)) * t64 + (-t27 * mrSges(6,3) + t117) * t63 + m(7) * (t63 * t78 - t102) + m(6) * (t17 * t63 - t102); -t90 * t67 + (mrSges(6,2) - t84) * t63 + m(7) * (-t31 * t67 + t32 * t82) + m(6) * (t33 * t67 + t34 * t63) + t116; m(5) + m(6) * (t55 + t58) + m(7) * (t55 * t87 + t58); -t114 * t8 - t115 * t77 - t103; t62 * t6 / 0.2e1 + t5 * t109 - pkin(5) * t11 - t101 + t86 * t27 + t78 * mrSges(7,3) + t73 * t28 - t114 * t15 + (m(7) * t78 + t117) * pkin(10); m(7) * (-pkin(5) * t31 + pkin(10) * t83) + t96 - t95 + (pkin(5) + t31) * t36 + (-pkin(10) * t87 + t83) * mrSges(7,3) - t74; -t63 * mrSges(6,2) + t114 * t67 + t115 * t82; pkin(5) * t110 + m(7) * (t87 * pkin(10) ^ 2 + pkin(5) ^ 2) + 0.2e1 * pkin(10) * t84 + t74; mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t1 - mrSges(7,2) * t2 - Ifges(7,6) * t98 + t89; -t32 * t76 - t39; -t76 * t63; -pkin(10) * t76 + t39; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t9(1) t9(2) t9(4) t9(7) t9(11) t9(16); t9(2) t9(3) t9(5) t9(8) t9(12) t9(17); t9(4) t9(5) t9(6) t9(9) t9(13) t9(18); t9(7) t9(8) t9(9) t9(10) t9(14) t9(19); t9(11) t9(12) t9(13) t9(14) t9(15) t9(20); t9(16) t9(17) t9(18) t9(19) t9(20) t9(21);];
Mq  = res;
