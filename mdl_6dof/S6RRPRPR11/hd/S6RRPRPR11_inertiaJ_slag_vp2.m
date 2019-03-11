% Calculate joint inertia matrix for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR11_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR11_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR11_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR11_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:11:46
% EndTime: 2019-03-09 11:11:48
% DurationCPUTime: 1.16s
% Computational Cost: add. (1848->282), mult. (3336->391), div. (0->0), fcn. (3377->8), ass. (0->113)
t162 = pkin(3) + pkin(7);
t161 = Ifges(5,3) + Ifges(6,3);
t112 = sin(qJ(2));
t115 = cos(qJ(2));
t160 = t112 ^ 2 + t115 ^ 2;
t110 = sin(qJ(6));
t113 = cos(qJ(6));
t108 = sin(pkin(10));
t109 = cos(pkin(10));
t111 = sin(qJ(4));
t114 = cos(qJ(4));
t121 = t108 * t111 - t109 * t114;
t74 = -t108 * t114 - t109 * t111;
t131 = t110 * t74 - t113 * t121;
t41 = -t110 * t121 - t113 * t74;
t159 = t131 ^ 2 + t41 ^ 2;
t116 = -pkin(2) - pkin(8);
t132 = -qJ(3) * t112 - pkin(1);
t72 = t116 * t115 + t132;
t87 = t162 * t112;
t78 = t114 * t87;
t30 = pkin(4) * t112 + t78 + (qJ(5) * t115 - t72) * t111;
t139 = t114 * t115;
t44 = t111 * t87 + t114 * t72;
t33 = -qJ(5) * t139 + t44;
t12 = -t108 * t33 + t109 * t30;
t59 = t74 * t115;
t4 = pkin(5) * t112 - pkin(9) * t59 + t12;
t13 = t108 * t30 + t109 * t33;
t58 = t121 * t115;
t5 = pkin(9) * t58 + t13;
t2 = -t110 * t5 + t113 * t4;
t3 = t110 * t4 + t113 * t5;
t158 = t131 * t2 + t3 * t41;
t138 = -qJ(5) + t116;
t81 = t138 * t111;
t82 = t138 * t114;
t48 = -t108 * t81 + t109 * t82;
t26 = pkin(9) * t121 + t48;
t49 = t108 * t82 + t109 * t81;
t27 = pkin(9) * t74 + t49;
t10 = t110 * t26 + t113 * t27;
t9 = -t110 * t27 + t113 * t26;
t157 = t10 * t41 + t131 * t9;
t146 = pkin(4) * t108;
t92 = pkin(4) * t109 + pkin(5);
t61 = -t110 * t146 + t113 * t92;
t62 = t110 * t92 + t113 * t146;
t156 = t131 * t61 + t41 * t62;
t154 = m(6) * pkin(4);
t151 = t2 * mrSges(7,1) - t3 * mrSges(7,2);
t149 = -m(4) * pkin(2) + mrSges(4,2);
t147 = -t111 / 0.2e1;
t145 = t61 * mrSges(7,1);
t144 = t62 * mrSges(7,2);
t143 = Ifges(5,4) * t111;
t142 = Ifges(5,4) * t114;
t141 = t114 * mrSges(5,1);
t93 = t111 * pkin(4) + qJ(3);
t140 = t160 * pkin(7) ^ 2;
t88 = t162 * t115;
t137 = t111 ^ 2 + t114 ^ 2;
t136 = t121 ^ 2 + t74 ^ 2;
t28 = -t110 * t59 + t113 * t58;
t29 = t110 * t58 + t113 * t59;
t135 = Ifges(7,5) * t29 + Ifges(7,6) * t28 + Ifges(7,3) * t112;
t64 = pkin(4) * t139 + t88;
t32 = -t58 * mrSges(6,1) + t59 * mrSges(6,2);
t45 = -t74 * mrSges(6,1) - mrSges(6,2) * t121;
t11 = -t28 * mrSges(7,1) + t29 * mrSges(7,2);
t14 = mrSges(7,1) * t41 + mrSges(7,2) * t131;
t134 = mrSges(7,1) * t131 - t41 * mrSges(7,2);
t133 = m(5) * t137;
t130 = t137 * mrSges(5,3);
t128 = -t12 * t121 - t13 * t74;
t127 = -t121 * t48 - t49 * t74;
t37 = Ifges(7,6) * t41;
t38 = Ifges(7,5) * t131;
t126 = t9 * mrSges(7,1) - t10 * mrSges(7,2) - t37 + t38;
t125 = -t111 * mrSges(5,2) + t141;
t124 = -Ifges(5,5) * t111 - Ifges(5,6) * t114;
t123 = t108 * t74 + t109 * t121;
t43 = -t111 * t72 + t78;
t122 = t44 * t111 + t43 * t114;
t120 = Ifges(6,5) * t59 + Ifges(6,6) * t58 + t112 * t161 + t135;
t117 = qJ(3) ^ 2;
t97 = Ifges(5,5) * t114;
t86 = Ifges(5,1) * t114 - t143;
t85 = -Ifges(5,2) * t111 + t142;
t84 = mrSges(5,1) * t111 + mrSges(5,2) * t114;
t83 = -pkin(2) * t115 + t132;
t80 = -mrSges(5,2) * t112 - mrSges(5,3) * t139;
t79 = mrSges(5,3) * t111 * t115 + mrSges(5,1) * t112;
t68 = t125 * t115;
t67 = Ifges(6,5) * t121;
t66 = Ifges(6,6) * t74;
t57 = Ifges(5,5) * t112 + (-Ifges(5,1) * t111 - t142) * t115;
t56 = t112 * Ifges(5,6) + (-Ifges(5,2) * t114 - t143) * t115;
t52 = -pkin(5) * t74 + t93;
t51 = mrSges(6,1) * t112 - mrSges(6,3) * t59;
t50 = -mrSges(6,2) * t112 + mrSges(6,3) * t58;
t47 = -Ifges(6,1) * t121 + Ifges(6,4) * t74;
t46 = -Ifges(6,4) * t121 + Ifges(6,2) * t74;
t34 = -pkin(5) * t58 + t64;
t25 = Ifges(6,1) * t59 + Ifges(6,4) * t58 + Ifges(6,5) * t112;
t24 = Ifges(6,4) * t59 + Ifges(6,2) * t58 + Ifges(6,6) * t112;
t18 = mrSges(7,1) * t112 - mrSges(7,3) * t29;
t17 = -mrSges(7,2) * t112 + mrSges(7,3) * t28;
t16 = Ifges(7,1) * t131 - Ifges(7,4) * t41;
t15 = Ifges(7,4) * t131 - Ifges(7,2) * t41;
t8 = Ifges(7,1) * t29 + Ifges(7,4) * t28 + Ifges(7,5) * t112;
t7 = Ifges(7,4) * t29 + Ifges(7,2) * t28 + Ifges(7,6) * t112;
t1 = [0.2e1 * t34 * t11 + 0.2e1 * t12 * t51 + 0.2e1 * t13 * t50 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 + t58 * t24 + t59 * t25 + t28 * t7 + t29 * t8 + 0.2e1 * t64 * t32 + 0.2e1 * t43 * t79 + 0.2e1 * t44 * t80 + 0.2e1 * t88 * t68 + Ifges(2,3) + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * t83 * mrSges(4,2) - t111 * t57 - t114 * t56 + (Ifges(4,3) + Ifges(3,2)) * t115) * t115 + (-0.2e1 * pkin(1) * mrSges(3,2) - 0.2e1 * t83 * mrSges(4,3) + (Ifges(3,1) + Ifges(4,2)) * t112 + ((2 * Ifges(3,4)) + (2 * Ifges(4,6)) + t124) * t115 + t120) * t112 + m(5) * (t43 ^ 2 + t44 ^ 2 + t88 ^ 2) + m(4) * (t83 ^ 2 + t140) + m(3) * (pkin(1) ^ 2 + t140) + m(6) * (t12 ^ 2 + t13 ^ 2 + t64 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t34 ^ 2) + 0.2e1 * (mrSges(4,1) + mrSges(3,3)) * pkin(7) * t160; -t121 * t25 / 0.2e1 + (-pkin(2) * mrSges(4,1) + Ifges(5,6) * t147 + t97 / 0.2e1 - t67 / 0.2e1 + t66 / 0.2e1 + t38 / 0.2e1 - t37 / 0.2e1 - Ifges(4,4) + Ifges(3,5)) * t112 - t41 * t7 / 0.2e1 + m(6) * (t12 * t48 + t13 * t49 + t64 * t93) + m(7) * (t10 * t3 + t2 * t9 + t34 * t52) + t88 * t84 + t93 * t32 + t64 * t45 + qJ(3) * t68 + t74 * t24 / 0.2e1 + t52 * t11 + t58 * t46 / 0.2e1 + t59 * t47 / 0.2e1 + t49 * t50 + t48 * t51 + t29 * t16 / 0.2e1 + t34 * t14 + t10 * t17 + t9 * t18 + t28 * t15 / 0.2e1 - t128 * mrSges(6,3) - t158 * mrSges(7,3) + m(5) * (qJ(3) * t88 + t122 * t116) + ((m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * t115 + (-mrSges(3,1) + t149) * t112) * pkin(7) + (t116 * t80 - t44 * mrSges(5,3) - t56 / 0.2e1) * t111 + (-t43 * mrSges(5,3) + t116 * t79 + t57 / 0.2e1) * t114 + t131 * t8 / 0.2e1 + (-t114 * t85 / 0.2e1 + t86 * t147 + qJ(3) * mrSges(4,1) - Ifges(4,5) + Ifges(3,6)) * t115; -0.2e1 * pkin(2) * mrSges(4,2) - t111 * t85 + t114 * t86 + 0.2e1 * t52 * t14 - t41 * t15 + t131 * t16 + 0.2e1 * t93 * t45 + t74 * t46 - t121 * t47 + Ifges(4,1) + Ifges(3,3) + m(7) * (t10 ^ 2 + t52 ^ 2 + t9 ^ 2) + m(6) * (t48 ^ 2 + t49 ^ 2 + t93 ^ 2) + m(5) * (t137 * t116 ^ 2 + t117) + m(4) * (pkin(2) ^ 2 + t117) + 0.2e1 * (t84 + mrSges(4,3)) * qJ(3) - 0.2e1 * t157 * mrSges(7,3) - 0.2e1 * t127 * mrSges(6,3) - 0.2e1 * t116 * t130; t111 * t80 + t114 * t79 + t41 * t17 + t131 * t18 - t74 * t50 - t121 * t51 + (m(4) * pkin(7) + mrSges(4,1)) * t112 + m(7) * t158 + m(6) * t128 + m(5) * t122; m(6) * t127 + m(7) * t157 - t136 * mrSges(6,3) - mrSges(7,3) * t159 + t116 * t133 - t130 + t149; m(6) * t136 + m(7) * t159 + m(4) + t133; m(7) * (t2 * t61 + t3 * t62) + t124 * t115 + (t109 * t51 + t108 * t50 + m(6) * (t108 * t13 + t109 * t12)) * pkin(4) + t120 + t62 * t17 + t61 * t18 + t43 * mrSges(5,1) - t44 * mrSges(5,2) + t12 * mrSges(6,1) - t13 * mrSges(6,2) + t151; m(7) * (t10 * t62 + t61 * t9) - t49 * mrSges(6,2) + t48 * mrSges(6,1) + t66 + t97 - t67 + t116 * t141 + (-mrSges(5,2) * t116 - Ifges(5,6)) * t111 - t156 * mrSges(7,3) + (m(6) * (t108 * t49 + t109 * t48) + t123 * mrSges(6,3)) * pkin(4) + t126; m(7) * t156 - t121 * mrSges(6,1) + t74 * mrSges(6,2) - t123 * t154 + t125 + t134; 0.2e1 * t145 - 0.2e1 * t144 + Ifges(7,3) + m(7) * (t61 ^ 2 + t62 ^ 2) + (0.2e1 * mrSges(6,1) * t109 - 0.2e1 * mrSges(6,2) * t108 + (t108 ^ 2 + t109 ^ 2) * t154) * pkin(4) + t161; m(6) * t64 + m(7) * t34 + t11 + t32; m(6) * t93 + m(7) * t52 + t14 + t45; 0; 0; m(6) + m(7); t135 + t151; t126; t134; Ifges(7,3) - t144 + t145; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
