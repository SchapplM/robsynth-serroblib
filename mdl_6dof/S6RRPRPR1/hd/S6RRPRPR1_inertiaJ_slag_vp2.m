% Calculate joint inertia matrix for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:07:56
% EndTime: 2019-03-09 10:07:58
% DurationCPUTime: 0.95s
% Computational Cost: add. (2859->259), mult. (5384->372), div. (0->0), fcn. (6287->10), ass. (0->99)
t112 = sin(qJ(4));
t108 = sin(pkin(10));
t110 = cos(pkin(10));
t113 = sin(qJ(2));
t139 = -qJ(3) - pkin(7);
t93 = t139 * t113;
t115 = cos(qJ(2));
t94 = t139 * t115;
t66 = t108 * t94 + t110 * t93;
t84 = t108 * t115 + t110 * t113;
t120 = -pkin(8) * t84 + t66;
t144 = cos(qJ(4));
t67 = t108 * t93 - t110 * t94;
t82 = -t108 * t113 + t110 * t115;
t50 = pkin(8) * t82 + t67;
t30 = t112 * t50 - t144 * t120;
t150 = t30 ^ 2;
t149 = 0.2e1 * t30;
t107 = sin(pkin(11));
t109 = cos(pkin(11));
t111 = sin(qJ(6));
t114 = cos(qJ(6));
t83 = -t107 * t111 + t109 * t114;
t85 = t107 * t114 + t109 * t111;
t59 = -t83 * mrSges(7,1) + t85 * mrSges(7,2);
t148 = 0.2e1 * t59;
t100 = -pkin(2) * t115 - pkin(1);
t69 = -pkin(3) * t82 + t100;
t147 = 0.2e1 * t69;
t146 = 0.2e1 * t82;
t143 = pkin(2) * t108;
t142 = t109 * pkin(5);
t98 = pkin(2) * t110 + pkin(3);
t73 = -t112 * t143 + t144 * t98;
t141 = t73 * mrSges(5,1);
t74 = t112 * t98 + t144 * t143;
t140 = t74 * mrSges(5,2);
t57 = t112 * t84 - t144 * t82;
t58 = t112 * t82 + t144 * t84;
t29 = pkin(4) * t57 - qJ(5) * t58 + t69;
t32 = t112 * t120 + t144 * t50;
t12 = t107 * t29 + t109 * t32;
t132 = t109 * t58;
t134 = t107 * t58;
t38 = mrSges(6,1) * t134 + mrSges(6,2) * t132;
t138 = Ifges(7,5) * t85 + Ifges(7,6) * t83;
t137 = Ifges(6,4) * t107;
t136 = Ifges(6,4) * t109;
t11 = -t107 * t32 + t109 * t29;
t135 = t107 * t11;
t133 = t109 * t12;
t71 = qJ(5) + t74;
t131 = t109 * t71;
t130 = t107 ^ 2 + t109 ^ 2;
t129 = t113 ^ 2 + t115 ^ 2;
t128 = 2 * mrSges(7,3);
t36 = t85 * t58;
t37 = t83 * t58;
t127 = Ifges(7,5) * t37 - Ifges(7,6) * t36 + Ifges(7,3) * t57;
t126 = -t82 * mrSges(4,1) + t84 * mrSges(4,2);
t13 = t36 * mrSges(7,1) + t37 * mrSges(7,2);
t89 = -t109 * mrSges(6,1) + t107 * mrSges(6,2);
t125 = t130 * qJ(5);
t60 = Ifges(7,4) * t85 + Ifges(7,2) * t83;
t61 = Ifges(7,1) * t85 + Ifges(7,4) * t83;
t91 = Ifges(6,2) * t109 + t137;
t92 = Ifges(6,1) * t107 + t136;
t124 = t107 * t92 + t109 * t91 + t83 * t60 + t85 * t61 + Ifges(5,3);
t123 = t133 - t135;
t122 = 0.2e1 * t130 * mrSges(6,3);
t72 = -pkin(4) - t73;
t121 = t59 + t89;
t15 = pkin(5) * t134 + t30;
t4 = pkin(5) * t57 - pkin(9) * t132 + t11;
t5 = -pkin(9) * t134 + t12;
t2 = -t111 * t5 + t114 * t4;
t22 = t57 * Ifges(6,6) + (-Ifges(6,2) * t107 + t136) * t58;
t23 = t57 * Ifges(6,5) + (Ifges(6,1) * t109 - t137) * t58;
t3 = t111 * t4 + t114 * t5;
t8 = Ifges(7,4) * t37 - Ifges(7,2) * t36 + Ifges(7,6) * t57;
t9 = Ifges(7,1) * t37 - Ifges(7,4) * t36 + Ifges(7,5) * t57;
t119 = -t32 * mrSges(5,2) + mrSges(6,3) * t133 + t15 * t59 + t107 * t23 / 0.2e1 + t109 * t22 / 0.2e1 - t36 * t60 / 0.2e1 + t37 * t61 / 0.2e1 - t91 * t134 / 0.2e1 + t92 * t132 / 0.2e1 - Ifges(5,6) * t57 + Ifges(5,5) * t58 + t83 * t8 / 0.2e1 + t85 * t9 / 0.2e1 + (t89 - mrSges(5,1)) * t30 + (Ifges(6,5) * t107 + Ifges(6,6) * t109 + t138) * t57 / 0.2e1 + (-t2 * t85 + t3 * t83) * mrSges(7,3);
t102 = t109 * pkin(9);
t99 = -pkin(4) - t142;
t90 = qJ(5) * t109 + t102;
t88 = (-pkin(9) - qJ(5)) * t107;
t68 = t72 - t142;
t65 = t102 + t131;
t64 = (-pkin(9) - t71) * t107;
t63 = t111 * t88 + t114 * t90;
t62 = -t111 * t90 + t114 * t88;
t53 = t58 * mrSges(5,2);
t43 = t111 * t64 + t114 * t65;
t42 = -t111 * t65 + t114 * t64;
t40 = mrSges(6,1) * t57 - mrSges(6,3) * t132;
t39 = -mrSges(6,2) * t57 - mrSges(6,3) * t134;
t17 = mrSges(7,1) * t57 - mrSges(7,3) * t37;
t16 = -mrSges(7,2) * t57 - mrSges(7,3) * t36;
t1 = [m(4) * (t100 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(7) * (t15 ^ 2 + t2 ^ 2 + t3 ^ 2) - 0.2e1 * pkin(1) * (-mrSges(3,1) * t115 + mrSges(3,2) * t113) + t115 * (Ifges(3,4) * t113 + Ifges(3,2) * t115) + t113 * (Ifges(3,1) * t113 + Ifges(3,4) * t115) - t36 * t8 + t37 * t9 + 0.2e1 * t12 * t39 + 0.2e1 * t11 * t40 + 0.2e1 * t15 * t13 + 0.2e1 * t3 * t16 + 0.2e1 * t2 * t17 + (mrSges(5,3) * t149 + Ifges(5,1) * t58 - t107 * t22 + t109 * t23 + (Ifges(6,5) * t109 - Ifges(6,6) * t107 - (2 * Ifges(5,4))) * t57) * t58 + m(5) * (t32 ^ 2 + t69 ^ 2 + t150) + m(6) * (t11 ^ 2 + t12 ^ 2 + t150) + (-0.2e1 * t66 * mrSges(4,3) + Ifges(4,1) * t84 + Ifges(4,4) * t146) * t84 + (mrSges(5,1) * t147 - 0.2e1 * t32 * mrSges(5,3) + (Ifges(6,3) + Ifges(5,2)) * t57 + t127) * t57 + m(3) * (t129 * pkin(7) ^ 2 + pkin(1) ^ 2) + 0.2e1 * t100 * t126 + 0.2e1 * t129 * pkin(7) * mrSges(3,3) + Ifges(2,3) + t67 * mrSges(4,3) * t146 + t53 * t147 + t38 * t149 + Ifges(4,2) * t82 ^ 2; t39 * t131 + t119 + m(6) * (t123 * t71 + t30 * t72) + Ifges(3,6) * t115 + Ifges(3,5) * t113 + Ifges(4,6) * t82 + Ifges(4,5) * t84 + t66 * mrSges(4,1) - t67 * mrSges(4,2) + t68 * t13 + t72 * t38 + t42 * t17 + t43 * t16 + (-t11 * mrSges(6,3) - t71 * t40) * t107 + m(7) * (t15 * t68 + t2 * t42 + t3 * t43) + m(5) * (-t30 * t73 + t32 * t74) + (-t113 * mrSges(3,1) - t115 * mrSges(3,2)) * pkin(7) + (-t74 * t57 - t73 * t58) * mrSges(5,3) + (m(4) * (t108 * t67 + t110 * t66) + (t108 * t82 - t110 * t84) * mrSges(4,3)) * pkin(2); 0.2e1 * t141 - 0.2e1 * t140 + t68 * t148 + 0.2e1 * t72 * t89 + Ifges(3,3) + Ifges(4,3) + (-t42 * t85 + t43 * t83) * t128 + t71 * t122 + m(7) * (t42 ^ 2 + t43 ^ 2 + t68 ^ 2) + m(6) * (t130 * t71 ^ 2 + t72 ^ 2) + m(5) * (t73 ^ 2 + t74 ^ 2) + t124 + (0.2e1 * mrSges(4,1) * t110 - 0.2e1 * mrSges(4,2) * t108 + m(4) * (t108 ^ 2 + t110 ^ 2) * pkin(2)) * pkin(2); t57 * mrSges(5,1) + t107 * t39 + t109 * t40 + t85 * t16 + t83 * t17 + t53 + m(7) * (t2 * t83 + t3 * t85) + m(6) * (t107 * t12 + t109 * t11) + m(5) * t69 + m(4) * t100 + t126; m(7) * (t42 * t83 + t43 * t85); m(4) + m(5) + m(6) * t130 + m(7) * (t83 ^ 2 + t85 ^ 2); -mrSges(6,3) * t135 + t119 + m(6) * (-pkin(4) * t30 + qJ(5) * t123) + t99 * t13 + t62 * t17 + t63 * t16 - pkin(4) * t38 + m(7) * (t15 * t99 + t2 * t62 + t3 * t63) + (-t107 * t40 + t109 * t39) * qJ(5); t141 - t140 + (t72 - pkin(4)) * t89 + (t68 + t99) * t59 + m(7) * (t42 * t62 + t43 * t63 + t68 * t99) + m(6) * (-pkin(4) * t72 + t125 * t71) + ((-t42 - t62) * t85 + (t43 + t63) * t83) * mrSges(7,3) + (t130 * t71 + t125) * mrSges(6,3) + t124; m(7) * (t62 * t83 + t63 * t85); -0.2e1 * pkin(4) * t89 + t99 * t148 + (-t62 * t85 + t63 * t83) * t128 + qJ(5) * t122 + m(7) * (t62 ^ 2 + t63 ^ 2 + t99 ^ 2) + m(6) * (t130 * qJ(5) ^ 2 + pkin(4) ^ 2) + t124; m(6) * t30 + m(7) * t15 + t13 + t38; m(6) * t72 + m(7) * t68 + t121; 0; -m(6) * pkin(4) + m(7) * t99 + t121; m(6) + m(7); mrSges(7,1) * t2 - mrSges(7,2) * t3 + t127; mrSges(7,1) * t42 - mrSges(7,2) * t43 + t138; -t59; mrSges(7,1) * t62 - t63 * mrSges(7,2) + t138; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
