% Calculate joint inertia matrix for
% S6RPRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:12:18
% EndTime: 2019-03-09 07:12:20
% DurationCPUTime: 1.06s
% Computational Cost: add. (3093->290), mult. (5954->429), div. (0->0), fcn. (6835->10), ass. (0->111)
t115 = cos(qJ(4));
t107 = sin(pkin(11));
t108 = cos(pkin(11));
t112 = sin(qJ(3));
t145 = cos(qJ(3));
t82 = t145 * t107 + t112 * t108;
t134 = t115 * t82;
t81 = t107 * t112 - t145 * t108;
t152 = Ifges(5,5) * t134 + Ifges(5,3) * t81;
t110 = sin(qJ(5));
t111 = sin(qJ(4));
t114 = cos(qJ(5));
t84 = t110 * t115 + t111 * t114;
t44 = t84 * t82;
t83 = -t110 * t111 + t114 * t115;
t45 = t83 * t82;
t151 = Ifges(6,5) * t45 - Ifges(6,6) * t44 + Ifges(6,3) * t81;
t141 = pkin(7) + qJ(2);
t88 = t141 * t107;
t89 = t141 * t108;
t61 = t112 * t89 + t145 * t88;
t150 = t61 ^ 2;
t149 = 0.2e1 * t61;
t96 = -pkin(2) * t108 - pkin(1);
t148 = 0.2e1 * t96;
t104 = t108 ^ 2;
t146 = -pkin(9) - pkin(8);
t144 = pkin(4) * t110;
t109 = sin(qJ(6));
t113 = cos(qJ(6));
t97 = pkin(4) * t114 + pkin(5);
t72 = t109 * t97 + t113 * t144;
t143 = t72 * mrSges(7,2);
t142 = Ifges(6,3) + Ifges(7,3);
t49 = pkin(3) * t81 - pkin(8) * t82 + t96;
t63 = -t112 * t88 + t145 * t89;
t30 = -t111 * t63 + t115 * t49;
t14 = pkin(4) * t81 - pkin(9) * t134 + t30;
t135 = t111 * t82;
t31 = t111 * t49 + t115 * t63;
t21 = -pkin(9) * t135 + t31;
t8 = t110 * t14 + t114 * t21;
t55 = -t109 * t84 + t113 * t83;
t56 = t109 * t83 + t113 * t84;
t140 = Ifges(7,5) * t56 + Ifges(7,6) * t55;
t139 = Ifges(6,5) * t84 + Ifges(6,6) * t83;
t93 = t146 * t111;
t94 = t146 * t115;
t66 = t110 * t93 - t114 * t94;
t138 = Ifges(5,4) * t111;
t137 = Ifges(5,4) * t115;
t136 = t109 * mrSges(7,2);
t133 = Ifges(5,5) * t111 + Ifges(5,6) * t115;
t132 = t107 ^ 2 + t104;
t131 = t111 ^ 2 + t115 ^ 2;
t130 = pkin(5) * t136;
t27 = -t109 * t45 - t113 * t44;
t28 = -t109 * t44 + t113 * t45;
t129 = Ifges(7,5) * t28 + Ifges(7,6) * t27 + Ifges(7,3) * t81;
t98 = -pkin(4) * t115 - pkin(3);
t58 = -t83 * mrSges(6,1) + t84 * mrSges(6,2);
t32 = -t55 * mrSges(7,1) + t56 * mrSges(7,2);
t7 = -t110 * t21 + t114 * t14;
t65 = t110 * t94 + t114 * t93;
t128 = -t108 * mrSges(3,1) + t107 * mrSges(3,2);
t71 = -t109 * t144 + t113 * t97;
t68 = t71 * mrSges(7,1);
t127 = Ifges(7,3) + t68 - t143;
t39 = pkin(4) * t135 + t61;
t90 = -t115 * mrSges(5,1) + t111 * mrSges(5,2);
t126 = t111 * mrSges(5,1) + t115 * mrSges(5,2);
t42 = -pkin(10) * t84 + t65;
t43 = pkin(10) * t83 + t66;
t24 = -t109 * t43 + t113 * t42;
t25 = t109 * t42 + t113 * t43;
t125 = t24 * mrSges(7,1) - t25 * mrSges(7,2) + t140;
t4 = pkin(5) * t81 - pkin(10) * t45 + t7;
t5 = -pkin(10) * t44 + t8;
t2 = -t109 * t5 + t113 * t4;
t3 = t109 * t4 + t113 * t5;
t124 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t129;
t123 = (mrSges(6,1) * t114 - mrSges(6,2) * t110) * pkin(4);
t122 = -t32 - t58;
t121 = t65 * mrSges(6,1) - t66 * mrSges(6,2) + t125 + t139;
t120 = t7 * mrSges(6,1) - t8 * mrSges(6,2) + t124 + t151;
t99 = t113 * pkin(5) * mrSges(7,1);
t92 = Ifges(5,1) * t111 + t137;
t91 = Ifges(5,2) * t115 + t138;
t73 = t82 * mrSges(4,2);
t67 = -pkin(5) * t83 + t98;
t60 = Ifges(6,1) * t84 + Ifges(6,4) * t83;
t59 = Ifges(6,4) * t84 + Ifges(6,2) * t83;
t54 = mrSges(5,1) * t81 - mrSges(5,3) * t134;
t53 = -mrSges(5,2) * t81 - mrSges(5,3) * t135;
t48 = t126 * t82;
t38 = Ifges(5,5) * t81 + (Ifges(5,1) * t115 - t138) * t82;
t37 = Ifges(5,6) * t81 + (-Ifges(5,2) * t111 + t137) * t82;
t36 = mrSges(6,1) * t81 - mrSges(6,3) * t45;
t35 = -mrSges(6,2) * t81 - mrSges(6,3) * t44;
t34 = Ifges(7,1) * t56 + Ifges(7,4) * t55;
t33 = Ifges(7,4) * t56 + Ifges(7,2) * t55;
t29 = mrSges(6,1) * t44 + mrSges(6,2) * t45;
t26 = pkin(5) * t44 + t39;
t18 = Ifges(6,1) * t45 - Ifges(6,4) * t44 + Ifges(6,5) * t81;
t17 = Ifges(6,4) * t45 - Ifges(6,2) * t44 + Ifges(6,6) * t81;
t16 = mrSges(7,1) * t81 - mrSges(7,3) * t28;
t15 = -mrSges(7,2) * t81 + mrSges(7,3) * t27;
t11 = -mrSges(7,1) * t27 + mrSges(7,2) * t28;
t10 = Ifges(7,1) * t28 + Ifges(7,4) * t27 + Ifges(7,5) * t81;
t9 = Ifges(7,4) * t28 + Ifges(7,2) * t27 + Ifges(7,6) * t81;
t1 = [m(6) * (t39 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(7) * (t2 ^ 2 + t26 ^ 2 + t3 ^ 2) + (Ifges(3,1) * t107 + 0.2e1 * Ifges(3,4) * t108) * t107 + 0.2e1 * t132 * qJ(2) * mrSges(3,3) + 0.2e1 * t30 * t54 - t44 * t17 + t45 * t18 + 0.2e1 * t31 * t53 + 0.2e1 * t8 * t35 + 0.2e1 * t7 * t36 + 0.2e1 * t39 * t29 + 0.2e1 * t26 * t11 + t27 * t9 + t28 * t10 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + (mrSges(4,3) * t149 + Ifges(4,1) * t82 - t111 * t37 + t115 * t38) * t82 + m(4) * (t63 ^ 2 + t96 ^ 2 + t150) + m(5) * (t30 ^ 2 + t31 ^ 2 + t150) + Ifges(3,2) * t104 + m(3) * (t132 * qJ(2) ^ 2 + pkin(1) ^ 2) - 0.2e1 * pkin(1) * t128 + (mrSges(4,1) * t148 - 0.2e1 * t63 * mrSges(4,3) + Ifges(4,2) * t81 + (-Ifges(5,6) * t111 - (2 * Ifges(4,4))) * t82 + t129 + t151 + t152) * t81 + Ifges(2,3) + t73 * t148 + t48 * t149; -m(3) * pkin(1) + t81 * mrSges(4,1) + t111 * t53 + t115 * t54 + t56 * t15 + t55 * t16 + t84 * t35 + t83 * t36 + t73 + m(7) * (t2 * t55 + t3 * t56) + m(6) * (t7 * t83 + t8 * t84) + m(5) * (t111 * t31 + t115 * t30) + m(4) * t96 + t128; m(3) + m(4) + m(5) * t131 + m(6) * (t83 ^ 2 + t84 ^ 2) + m(7) * (t55 ^ 2 + t56 ^ 2); m(5) * (-pkin(3) * t61 + (-t30 * t111 + t31 * t115) * pkin(8)) + t98 * t29 + t83 * t17 / 0.2e1 + t84 * t18 / 0.2e1 - Ifges(4,6) * t81 + Ifges(4,5) * t82 + t66 * t35 + t67 * t11 + t55 * t9 / 0.2e1 + t56 * t10 / 0.2e1 + t39 * t58 - t44 * t59 / 0.2e1 + t45 * t60 / 0.2e1 - t63 * mrSges(4,2) + t65 * t36 - pkin(3) * t48 + t26 * t32 + t27 * t33 / 0.2e1 + t28 * t34 / 0.2e1 + t24 * t16 + t25 * t15 + (t140 + t139 + t133) * t81 / 0.2e1 + (t90 - mrSges(4,1)) * t61 + (t37 / 0.2e1 + t82 * t92 / 0.2e1 + pkin(8) * t53 + t31 * mrSges(5,3)) * t115 + (t38 / 0.2e1 - t82 * t91 / 0.2e1 - pkin(8) * t54 - t30 * mrSges(5,3)) * t111 + m(6) * (t39 * t98 + t65 * t7 + t66 * t8) + m(7) * (t2 * t24 + t25 * t3 + t26 * t67) + (-t2 * t56 + t3 * t55) * mrSges(7,3) + (-t7 * t84 + t8 * t83) * mrSges(6,3); m(6) * (t65 * t83 + t66 * t84) + m(7) * (t24 * t55 + t25 * t56); -0.2e1 * pkin(3) * t90 + t111 * t92 + t115 * t91 + 0.2e1 * t67 * t32 + t55 * t33 + t56 * t34 + 0.2e1 * t98 * t58 + t83 * t59 + t84 * t60 + Ifges(4,3) + m(7) * (t24 ^ 2 + t25 ^ 2 + t67 ^ 2) + m(6) * (t65 ^ 2 + t66 ^ 2 + t98 ^ 2) + m(5) * (t131 * pkin(8) ^ 2 + pkin(3) ^ 2) + 0.2e1 * (-t24 * t56 + t25 * t55) * mrSges(7,3) + 0.2e1 * (-t65 * t84 + t66 * t83) * mrSges(6,3) + 0.2e1 * t131 * pkin(8) * mrSges(5,3); (m(6) * (t110 * t8 + t114 * t7) + t114 * t36 + t110 * t35) * pkin(4) + t71 * t16 + t72 * t15 + t30 * mrSges(5,1) - t31 * mrSges(5,2) - Ifges(5,6) * t135 + t120 + m(7) * (t2 * t71 + t3 * t72) + t152; m(7) * (t55 * t71 + t56 * t72) + m(6) * (t110 * t84 + t114 * t83) * pkin(4) + t122 - t90; m(7) * (t24 * t71 + t25 * t72) - t126 * pkin(8) + (t72 * t55 - t71 * t56) * mrSges(7,3) + (m(6) * (t110 * t66 + t114 * t65) + (t110 * t83 - t114 * t84) * mrSges(6,3)) * pkin(4) + t121 + t133; -0.2e1 * t143 + Ifges(5,3) + 0.2e1 * t68 + 0.2e1 * t123 + m(7) * (t71 ^ 2 + t72 ^ 2) + m(6) * (t110 ^ 2 + t114 ^ 2) * pkin(4) ^ 2 + t142; (m(7) * (t109 * t3 + t113 * t2) + t109 * t15 + t113 * t16) * pkin(5) + t120; m(7) * (t109 * t56 + t113 * t55) * pkin(5) + t122; (m(7) * (t109 * t25 + t113 * t24) + (t109 * t55 - t113 * t56) * mrSges(7,3)) * pkin(5) + t121; Ifges(6,3) + t99 + t123 + (m(7) * (t109 * t72 + t113 * t71) - t136) * pkin(5) + t127; -0.2e1 * t130 + 0.2e1 * t99 + m(7) * (t109 ^ 2 + t113 ^ 2) * pkin(5) ^ 2 + t142; t124; -t32; t125; t127; Ifges(7,3) + t99 - t130; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
