% Calculate joint inertia matrix for
% S6RRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2018-11-23 17:22
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:21:46
% EndTime: 2018-11-23 17:21:47
% DurationCPUTime: 1.16s
% Computational Cost: add. (3244->303), mult. (6137->452), div. (0->0), fcn. (6999->10), ass. (0->111)
t116 = cos(qJ(4));
t108 = sin(pkin(11));
t109 = cos(pkin(11));
t113 = sin(qJ(2));
t117 = cos(qJ(2));
t86 = t108 * t117 + t109 * t113;
t136 = t116 * t86;
t85 = t108 * t113 - t109 * t117;
t153 = Ifges(5,5) * t136 + Ifges(5,3) * t85;
t111 = sin(qJ(5));
t112 = sin(qJ(4));
t115 = cos(qJ(5));
t89 = t111 * t116 + t112 * t115;
t44 = t89 * t86;
t88 = -t111 * t112 + t115 * t116;
t45 = t88 * t86;
t152 = Ifges(6,5) * t45 - Ifges(6,6) * t44 + Ifges(6,3) * t85;
t143 = -qJ(3) - pkin(7);
t130 = t143 * t113;
t92 = t143 * t117;
t64 = -t108 * t92 - t109 * t130;
t151 = t64 ^ 2;
t150 = 0.2e1 * t64;
t100 = -pkin(2) * t117 - pkin(1);
t149 = 0.2e1 * t100;
t97 = pkin(2) * t108 + pkin(8);
t147 = pkin(9) + t97;
t146 = pkin(4) * t111;
t110 = sin(qJ(6));
t114 = cos(qJ(6));
t99 = pkin(4) * t115 + pkin(5);
t75 = t110 * t99 + t114 * t146;
t145 = t75 * mrSges(7,2);
t144 = Ifges(6,3) + Ifges(7,3);
t52 = pkin(3) * t85 - pkin(8) * t86 + t100;
t66 = t108 * t130 - t109 * t92;
t30 = -t112 * t66 + t116 * t52;
t19 = pkin(4) * t85 - pkin(9) * t136 + t30;
t137 = t112 * t86;
t31 = t112 * t52 + t116 * t66;
t25 = -pkin(9) * t137 + t31;
t8 = t111 * t19 + t115 * t25;
t58 = -t110 * t89 + t114 * t88;
t59 = t110 * t88 + t114 * t89;
t142 = Ifges(7,5) * t59 + Ifges(7,6) * t58;
t80 = t147 * t112;
t81 = t147 * t116;
t51 = -t111 * t80 + t115 * t81;
t141 = Ifges(6,5) * t89 + Ifges(6,6) * t88;
t140 = Ifges(5,4) * t112;
t139 = Ifges(5,4) * t116;
t138 = t110 * mrSges(7,2);
t135 = Ifges(5,5) * t112 + Ifges(5,6) * t116;
t134 = t112 ^ 2 + t116 ^ 2;
t133 = t113 ^ 2 + t117 ^ 2;
t132 = pkin(5) * t138;
t26 = -t110 * t45 - t114 * t44;
t27 = -t110 * t44 + t114 * t45;
t131 = Ifges(7,5) * t27 + Ifges(7,6) * t26 + Ifges(7,3) * t85;
t98 = -pkin(2) * t109 - pkin(3);
t61 = -t88 * mrSges(6,1) + t89 * mrSges(6,2);
t32 = -t58 * mrSges(7,1) + t59 * mrSges(7,2);
t7 = -t111 * t25 + t115 * t19;
t50 = -t111 * t81 - t115 * t80;
t74 = -t110 * t146 + t114 * t99;
t68 = t74 * mrSges(7,1);
t129 = Ifges(7,3) + t68 - t145;
t41 = pkin(4) * t137 + t64;
t91 = -t116 * mrSges(5,1) + t112 * mrSges(5,2);
t128 = t112 * mrSges(5,1) + t116 * mrSges(5,2);
t39 = -pkin(10) * t89 + t50;
t40 = pkin(10) * t88 + t51;
t13 = -t110 * t40 + t114 * t39;
t14 = t110 * t39 + t114 * t40;
t127 = t13 * mrSges(7,1) - t14 * mrSges(7,2) + t142;
t90 = -pkin(4) * t116 + t98;
t4 = pkin(5) * t85 - pkin(10) * t45 + t7;
t5 = -pkin(10) * t44 + t8;
t2 = -t110 * t5 + t114 * t4;
t3 = t110 * t4 + t114 * t5;
t126 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t131;
t125 = (mrSges(6,1) * t115 - mrSges(6,2) * t111) * pkin(4);
t124 = -t32 - t61;
t123 = t50 * mrSges(6,1) - t51 * mrSges(6,2) + t127 + t141;
t122 = t7 * mrSges(6,1) - t8 * mrSges(6,2) + t126 + t152;
t101 = t114 * pkin(5) * mrSges(7,1);
t94 = Ifges(5,1) * t112 + t139;
t93 = Ifges(5,2) * t116 + t140;
t76 = t86 * mrSges(4,2);
t67 = -pkin(5) * t88 + t90;
t63 = Ifges(6,1) * t89 + Ifges(6,4) * t88;
t62 = Ifges(6,4) * t89 + Ifges(6,2) * t88;
t57 = mrSges(5,1) * t85 - mrSges(5,3) * t136;
t56 = -mrSges(5,2) * t85 - mrSges(5,3) * t137;
t49 = t128 * t86;
t38 = Ifges(5,5) * t85 + (Ifges(5,1) * t116 - t140) * t86;
t37 = Ifges(5,6) * t85 + (-Ifges(5,2) * t112 + t139) * t86;
t36 = mrSges(6,1) * t85 - mrSges(6,3) * t45;
t35 = -mrSges(6,2) * t85 - mrSges(6,3) * t44;
t34 = Ifges(7,1) * t59 + Ifges(7,4) * t58;
t33 = Ifges(7,4) * t59 + Ifges(7,2) * t58;
t29 = mrSges(6,1) * t44 + mrSges(6,2) * t45;
t28 = pkin(5) * t44 + t41;
t21 = Ifges(6,1) * t45 - Ifges(6,4) * t44 + Ifges(6,5) * t85;
t20 = Ifges(6,4) * t45 - Ifges(6,2) * t44 + Ifges(6,6) * t85;
t18 = mrSges(7,1) * t85 - mrSges(7,3) * t27;
t17 = -mrSges(7,2) * t85 + mrSges(7,3) * t26;
t11 = -mrSges(7,1) * t26 + mrSges(7,2) * t27;
t10 = Ifges(7,1) * t27 + Ifges(7,4) * t26 + Ifges(7,5) * t85;
t9 = Ifges(7,4) * t27 + Ifges(7,2) * t26 + Ifges(7,6) * t85;
t1 = [t117 * (Ifges(3,4) * t113 + Ifges(3,2) * t117) - 0.2e1 * pkin(1) * (-mrSges(3,1) * t117 + mrSges(3,2) * t113) + t113 * (Ifges(3,1) * t113 + Ifges(3,4) * t117) + (mrSges(4,1) * t149 - 0.2e1 * t66 * mrSges(4,3) + Ifges(4,2) * t85 + (-Ifges(5,6) * t112 - (2 * Ifges(4,4))) * t86 + t131 + t152 + t153) * t85 + t76 * t149 + t49 * t150 - t44 * t20 + t45 * t21 + 0.2e1 * t31 * t56 + 0.2e1 * t30 * t57 + 0.2e1 * t8 * t35 + 0.2e1 * t7 * t36 + 0.2e1 * t41 * t29 + t26 * t9 + t27 * t10 + 0.2e1 * t28 * t11 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 + m(6) * (t41 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(7) * (t2 ^ 2 + t28 ^ 2 + t3 ^ 2) + (mrSges(4,3) * t150 + Ifges(4,1) * t86 - t112 * t37 + t116 * t38) * t86 + m(4) * (t100 ^ 2 + t66 ^ 2 + t151) + m(5) * (t30 ^ 2 + t31 ^ 2 + t151) + m(3) * (pkin(7) ^ 2 * t133 + pkin(1) ^ 2) + 0.2e1 * t133 * pkin(7) * mrSges(3,3) + Ifges(2,3); Ifges(3,6) * t117 + Ifges(3,5) * t113 + t98 * t49 + t88 * t20 / 0.2e1 + t89 * t21 / 0.2e1 + t90 * t29 + m(5) * (t64 * t98 + (-t30 * t112 + t31 * t116) * t97) + (t91 - mrSges(4,1)) * t64 + (t37 / 0.2e1 + t86 * t94 / 0.2e1 + t97 * t56 + t31 * mrSges(5,3)) * t116 + (t38 / 0.2e1 - t86 * t93 / 0.2e1 - t97 * t57 - t30 * mrSges(5,3)) * t112 - Ifges(4,6) * t85 + Ifges(4,5) * t86 + t58 * t9 / 0.2e1 + t59 * t10 / 0.2e1 + t41 * t61 - t44 * t62 / 0.2e1 + t45 * t63 / 0.2e1 - t66 * mrSges(4,2) + t67 * t11 + t50 * t36 + t51 * t35 + t28 * t32 + t26 * t33 / 0.2e1 + t27 * t34 / 0.2e1 + t14 * t17 + t13 * t18 + m(6) * (t41 * t90 + t50 * t7 + t51 * t8) + m(7) * (t13 * t2 + t14 * t3 + t28 * t67) + (-t113 * mrSges(3,1) - t117 * mrSges(3,2)) * pkin(7) + (-t2 * t59 + t3 * t58) * mrSges(7,3) + (-t7 * t89 + t8 * t88) * mrSges(6,3) + (m(4) * (t108 * t66 - t109 * t64) + (-t108 * t85 - t109 * t86) * mrSges(4,3)) * pkin(2) + (t142 + t141 + t135) * t85 / 0.2e1; t112 * t94 + t116 * t93 + 0.2e1 * t67 * t32 + t58 * t33 + t59 * t34 + 0.2e1 * t90 * t61 + t88 * t62 + t89 * t63 + 0.2e1 * t98 * t91 + Ifges(3,3) + Ifges(4,3) + m(7) * (t13 ^ 2 + t14 ^ 2 + t67 ^ 2) + m(6) * (t50 ^ 2 + t51 ^ 2 + t90 ^ 2) + m(5) * (t134 * t97 ^ 2 + t98 ^ 2) + m(4) * (t108 ^ 2 + t109 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (t109 * mrSges(4,1) - t108 * mrSges(4,2)) * pkin(2) + 0.2e1 * (-t13 * t59 + t14 * t58) * mrSges(7,3) + 0.2e1 * (-t50 * t89 + t51 * t88) * mrSges(6,3) + 0.2e1 * t134 * t97 * mrSges(5,3); t85 * mrSges(4,1) + t112 * t56 + t116 * t57 + t59 * t17 + t58 * t18 + t89 * t35 + t88 * t36 + t76 + m(7) * (t2 * t58 + t3 * t59) + m(6) * (t7 * t88 + t8 * t89) + m(5) * (t112 * t31 + t116 * t30) + m(4) * t100; m(7) * (t13 * t58 + t14 * t59) + m(6) * (t50 * t88 + t51 * t89); m(4) + m(5) * t134 + m(6) * (t88 ^ 2 + t89 ^ 2) + m(7) * (t58 ^ 2 + t59 ^ 2); -Ifges(5,6) * t137 + (m(6) * (t111 * t8 + t115 * t7) + t115 * t36 + t111 * t35) * pkin(4) + t74 * t18 + t75 * t17 + t30 * mrSges(5,1) - t31 * mrSges(5,2) + t122 + m(7) * (t2 * t74 + t3 * t75) + t153; m(7) * (t13 * t74 + t14 * t75) - t128 * t97 + (t75 * t58 - t74 * t59) * mrSges(7,3) + (m(6) * (t111 * t51 + t115 * t50) + (t111 * t88 - t115 * t89) * mrSges(6,3)) * pkin(4) + t123 + t135; m(7) * (t58 * t74 + t59 * t75) + m(6) * (t111 * t89 + t115 * t88) * pkin(4) + t124 - t91; -0.2e1 * t145 + Ifges(5,3) + 0.2e1 * t68 + 0.2e1 * t125 + m(7) * (t74 ^ 2 + t75 ^ 2) + m(6) * (t111 ^ 2 + t115 ^ 2) * pkin(4) ^ 2 + t144; (m(7) * (t110 * t3 + t114 * t2) + t110 * t17 + t114 * t18) * pkin(5) + t122; (m(7) * (t110 * t14 + t114 * t13) + (t110 * t58 - t114 * t59) * mrSges(7,3)) * pkin(5) + t123; m(7) * (t110 * t59 + t114 * t58) * pkin(5) + t124; Ifges(6,3) + t101 + t125 + (-t138 + m(7) * (t110 * t75 + t114 * t74)) * pkin(5) + t129; -0.2e1 * t132 + 0.2e1 * t101 + m(7) * (t110 ^ 2 + t114 ^ 2) * pkin(5) ^ 2 + t144; t126; t127; -t32; t129; Ifges(7,3) + t101 - t132; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
