% Calculate joint inertia matrix for
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:40:39
% EndTime: 2019-03-09 03:40:41
% DurationCPUTime: 0.97s
% Computational Cost: add. (1678->260), mult. (3319->374), div. (0->0), fcn. (3420->10), ass. (0->105)
t100 = cos(qJ(5));
t92 = sin(pkin(11));
t94 = cos(pkin(11));
t97 = sin(qJ(5));
t69 = t100 * t92 + t94 * t97;
t98 = sin(qJ(3));
t55 = t69 * t98;
t68 = t100 * t94 - t92 * t97;
t56 = t68 * t98;
t32 = t55 * mrSges(6,1) + mrSges(6,2) * t56;
t96 = sin(qJ(6));
t99 = cos(qJ(6));
t26 = -t55 * t99 - t56 * t96;
t27 = -t55 * t96 + t56 * t99;
t9 = t26 * mrSges(7,1) - mrSges(7,2) * t27;
t105 = -t32 + t9;
t121 = t94 * t98;
t122 = t92 * t98;
t60 = mrSges(5,1) * t122 + mrSges(5,2) * t121;
t131 = -t60 + t105;
t93 = sin(pkin(10));
t84 = pkin(1) * t93 + pkin(7);
t130 = 0.2e1 * t84;
t38 = t68 * t99 - t69 * t96;
t39 = t68 * t96 + t69 * t99;
t14 = -t38 * mrSges(7,1) + mrSges(7,2) * t39;
t42 = -t68 * mrSges(6,1) + mrSges(6,2) * t69;
t129 = -t14 - t42;
t128 = -Ifges(6,5) * t56 + Ifges(6,6) * t55;
t127 = -t92 / 0.2e1;
t126 = t94 / 0.2e1;
t125 = Ifges(5,4) * t92;
t124 = Ifges(5,4) * t94;
t101 = cos(qJ(3));
t123 = pkin(3) * t101;
t120 = t98 * mrSges(5,3);
t79 = t98 * t84;
t119 = -Ifges(7,3) - Ifges(6,3);
t118 = pkin(8) + qJ(4);
t117 = -Ifges(7,5) * t27 - Ifges(7,6) * t26;
t95 = cos(pkin(10));
t85 = -pkin(1) * t95 - pkin(2);
t64 = -qJ(4) * t98 - t123 + t85;
t58 = t94 * t64;
t28 = -pkin(8) * t121 + t58 + (-t84 * t92 - pkin(4)) * t101;
t113 = t101 * t84;
t41 = t113 * t94 + t64 * t92;
t33 = -pkin(8) * t122 + t41;
t13 = t100 * t33 + t28 * t97;
t73 = t118 * t92;
t75 = t118 * t94;
t46 = t100 * t75 - t73 * t97;
t74 = -t94 * mrSges(5,1) + mrSges(5,2) * t92;
t116 = t74 - mrSges(4,1);
t59 = pkin(4) * t122 + t79;
t115 = t92 ^ 2 + t94 ^ 2;
t90 = t98 ^ 2;
t91 = t101 ^ 2;
t114 = t90 + t91;
t112 = m(5) + m(6) + m(7);
t86 = -pkin(4) * t94 - pkin(3);
t12 = t100 * t28 - t33 * t97;
t45 = -t100 * t73 - t75 * t97;
t111 = qJ(4) * t115;
t40 = -t113 * t92 + t58;
t110 = -t40 * t92 + t41 * t94;
t70 = mrSges(5,2) * t101 - t120 * t92;
t71 = -mrSges(5,1) * t101 - t120 * t94;
t109 = t94 * t70 - t92 * t71;
t4 = -pkin(5) * t101 - pkin(9) * t56 + t12;
t5 = -pkin(9) * t55 + t13;
t2 = t4 * t99 - t5 * t96;
t3 = t4 * t96 + t5 * t99;
t108 = mrSges(7,1) * t2 - t3 * mrSges(7,2) - t117;
t29 = -pkin(9) * t69 + t45;
t30 = pkin(9) * t68 + t46;
t10 = t29 * t99 - t30 * t96;
t11 = t29 * t96 + t30 * t99;
t36 = Ifges(7,6) * t38;
t37 = Ifges(7,5) * t39;
t107 = mrSges(7,1) * t10 - t11 * mrSges(7,2) + t36 + t37;
t106 = (mrSges(7,1) * t99 - mrSges(7,2) * t96) * pkin(5);
t83 = t84 ^ 2;
t78 = t90 * t83;
t77 = Ifges(5,1) * t92 + t124;
t76 = Ifges(5,2) * t94 + t125;
t63 = Ifges(6,5) * t69;
t62 = Ifges(6,6) * t68;
t54 = -Ifges(5,5) * t101 + (Ifges(5,1) * t94 - t125) * t98;
t53 = -Ifges(5,6) * t101 + (-Ifges(5,2) * t92 + t124) * t98;
t49 = -pkin(5) * t68 + t86;
t48 = -mrSges(6,1) * t101 - mrSges(6,3) * t56;
t47 = mrSges(6,2) * t101 - mrSges(6,3) * t55;
t44 = Ifges(6,1) * t69 + Ifges(6,4) * t68;
t43 = Ifges(6,4) * t69 + Ifges(6,2) * t68;
t34 = pkin(5) * t55 + t59;
t25 = Ifges(6,1) * t56 - Ifges(6,4) * t55 - Ifges(6,5) * t101;
t24 = Ifges(6,4) * t56 - Ifges(6,2) * t55 - Ifges(6,6) * t101;
t18 = -mrSges(7,1) * t101 - mrSges(7,3) * t27;
t17 = mrSges(7,2) * t101 + mrSges(7,3) * t26;
t16 = Ifges(7,1) * t39 + Ifges(7,4) * t38;
t15 = Ifges(7,4) * t39 + Ifges(7,2) * t38;
t7 = Ifges(7,1) * t27 + Ifges(7,4) * t26 - Ifges(7,5) * t101;
t6 = Ifges(7,4) * t27 + Ifges(7,2) * t26 - Ifges(7,6) * t101;
t1 = [0.2e1 * t12 * t48 + 0.2e1 * t13 * t47 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 - t55 * t24 + t56 * t25 + t26 * t6 + t27 * t7 + 0.2e1 * t59 * t32 - 0.2e1 * t34 * t9 + 0.2e1 * t40 * t71 + 0.2e1 * t41 * t70 + Ifges(2,3) + Ifges(3,3) + (0.2e1 * t85 * mrSges(4,2) + Ifges(4,1) * t98 + t60 * t130 - t92 * t53 + t94 * t54) * t98 + (-0.2e1 * t85 * mrSges(4,1) + (Ifges(5,3) + Ifges(4,2) - t119) * t101 + (-Ifges(5,5) * t94 + Ifges(5,6) * t92 + (2 * Ifges(4,4))) * t98 + t117 + t128) * t101 + m(4) * (t83 * t91 + t85 ^ 2 + t78) + m(6) * (t12 ^ 2 + t13 ^ 2 + t59 ^ 2) + m(5) * (t40 ^ 2 + t41 ^ 2 + t78) + m(7) * (t2 ^ 2 + t3 ^ 2 + t34 ^ 2) + m(3) * (t93 ^ 2 + t95 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (mrSges(3,1) * t95 - mrSges(3,2) * t93) * pkin(1) + t114 * mrSges(4,3) * t130; t27 * t17 + t26 * t18 + t56 * t47 - t55 * t48 + t131 * t101 + m(7) * (-t101 * t34 + t2 * t26 + t27 * t3) + m(6) * (-t101 * t59 - t12 * t55 + t13 * t56) + (m(5) * (t110 - t113) + t109) * t98; m(3) + m(7) * (t26 ^ 2 + t27 ^ 2 + t91) + m(6) * (t55 ^ 2 + t56 ^ 2 + t91) + m(5) * (t115 * t90 + t91) + m(4) * t114; (-t2 * t39 + t3 * t38) * mrSges(7,3) + (-t12 * t69 + t13 * t68) * mrSges(6,3) + m(5) * (-pkin(3) * t79 + qJ(4) * t110) + t109 * qJ(4) + t110 * mrSges(5,3) + t92 * t54 / 0.2e1 + t69 * t25 / 0.2e1 + t86 * t32 + t59 * t42 - pkin(3) * t60 + t68 * t24 / 0.2e1 + t46 * t47 + t45 * t48 - t49 * t9 - t55 * t43 / 0.2e1 + t56 * t44 / 0.2e1 + t34 * t14 + t38 * t6 / 0.2e1 + t39 * t7 / 0.2e1 + t26 * t15 / 0.2e1 + t27 * t16 / 0.2e1 + t11 * t17 + t10 * t18 + (t116 * t84 + t126 * t77 + t127 * t76 + Ifges(4,5)) * t98 + (-t84 * mrSges(4,2) + Ifges(5,5) * t127 - Ifges(5,6) * t94 / 0.2e1 - t63 / 0.2e1 - t62 / 0.2e1 - t37 / 0.2e1 - t36 / 0.2e1 + Ifges(4,6)) * t101 + t53 * t126 + m(6) * (t12 * t45 + t13 * t46 + t59 * t86) + m(7) * (t10 * t2 + t11 * t3 + t34 * t49); (-t26 * t39 + t27 * t38) * mrSges(7,3) + (t55 * t69 + t56 * t68) * mrSges(6,3) + (mrSges(5,3) * t115 - mrSges(4,2)) * t98 + (-t116 + t129) * t101 + m(7) * (t10 * t26 - t101 * t49 + t11 * t27) + m(6) * (-t101 * t86 - t45 * t55 + t46 * t56) + m(5) * (t111 * t98 + t123); -0.2e1 * pkin(3) * t74 + 0.2e1 * t49 * t14 + t38 * t15 + t39 * t16 + 0.2e1 * t86 * t42 + t68 * t43 + t69 * t44 + t94 * t76 + t92 * t77 + Ifges(4,3) + m(7) * (t10 ^ 2 + t11 ^ 2 + t49 ^ 2) + m(6) * (t45 ^ 2 + t46 ^ 2 + t86 ^ 2) + m(5) * (qJ(4) ^ 2 * t115 + pkin(3) ^ 2) + 0.2e1 * (-t10 * t39 + t11 * t38) * mrSges(7,3) + 0.2e1 * (-t45 * t69 + t46 * t68) * mrSges(6,3) + 0.2e1 * mrSges(5,3) * t111; m(5) * t79 + m(6) * t59 + m(7) * t34 - t131; -t112 * t101; -m(5) * pkin(3) + m(6) * t86 + m(7) * t49 - t129 + t74; t112; t12 * mrSges(6,1) - t13 * mrSges(6,2) + t119 * t101 + (m(7) * (t2 * t99 + t3 * t96) + t99 * t18 + t96 * t17) * pkin(5) + t108 - t128; m(7) * (t26 * t99 + t27 * t96) * pkin(5) + t105; t45 * mrSges(6,1) - t46 * mrSges(6,2) + t62 + t63 + (m(7) * (t10 * t99 + t11 * t96) + (t38 * t96 - t39 * t99) * mrSges(7,3)) * pkin(5) + t107; 0; m(7) * (t96 ^ 2 + t99 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t106 - t119; -Ifges(7,3) * t101 + t108; t9; t107; 0; Ifges(7,3) + t106; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
