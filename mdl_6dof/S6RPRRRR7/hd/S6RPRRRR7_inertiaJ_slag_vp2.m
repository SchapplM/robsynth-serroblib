% Calculate joint inertia matrix for
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 16:37
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:36:41
% EndTime: 2018-11-23 16:36:42
% DurationCPUTime: 0.90s
% Computational Cost: add. (1936->219), mult. (3337->315), div. (0->0), fcn. (3599->8), ass. (0->106)
t88 = sin(qJ(6));
t92 = cos(qJ(6));
t65 = -mrSges(7,1) * t92 + mrSges(7,2) * t88;
t159 = t65 - mrSges(6,1);
t90 = sin(qJ(4));
t91 = sin(qJ(3));
t94 = cos(qJ(4));
t95 = cos(qJ(3));
t60 = -t90 * t95 - t94 * t91;
t62 = -t90 * t91 + t94 * t95;
t89 = sin(qJ(5));
t93 = cos(qJ(5));
t109 = -t93 * t60 + t62 * t89;
t35 = t60 * t89 + t93 * t62;
t158 = t109 * t89 + t35 * t93;
t144 = pkin(3) * t90;
t76 = pkin(3) * t94 + pkin(4);
t47 = -t144 * t89 + t76 * t93;
t48 = t93 * t144 + t89 * t76;
t157 = t109 * t48 + t35 * t47;
t156 = t95 ^ 2;
t131 = t88 * mrSges(7,3);
t15 = -mrSges(7,2) * t109 - t131 * t35;
t133 = t35 * t92;
t16 = mrSges(7,1) * t109 - mrSges(7,3) * t133;
t111 = t92 * t15 - t88 * t16;
t142 = pkin(5) * t65;
t83 = t88 ^ 2;
t139 = mrSges(7,3) * t83;
t77 = pkin(10) * t139;
t85 = t92 ^ 2;
t138 = mrSges(7,3) * t85;
t78 = pkin(10) * t138;
t153 = t77 + t78 - t142;
t96 = -pkin(1) - pkin(7);
t140 = -pkin(8) + t96;
t118 = t140 * t95;
t119 = t140 * t91;
t37 = t118 * t94 - t119 * t90;
t103 = -t62 * pkin(9) + t37;
t38 = t90 * t118 + t94 * t119;
t21 = pkin(9) * t60 + t38;
t10 = -t103 * t93 + t21 * t89;
t112 = mrSges(7,1) * t88 + mrSges(7,2) * t92;
t14 = t112 * t35;
t152 = m(7) * t10 + t14;
t143 = pkin(4) * t93;
t75 = -pkin(5) - t143;
t51 = t75 * t65;
t74 = pkin(4) * t89 + pkin(10);
t63 = t74 * t139;
t64 = t74 * t138;
t79 = mrSges(6,1) * t143;
t151 = t51 + t63 + t64 + t79;
t12 = t103 * t89 + t93 * t21;
t71 = t91 * pkin(3) + qJ(2);
t42 = -pkin(4) * t60 + t71;
t13 = pkin(5) * t109 - pkin(10) * t35 + t42;
t3 = t12 * t92 + t13 * t88;
t141 = t3 * t92;
t2 = -t12 * t88 + t13 * t92;
t113 = -t2 * t88 + t141;
t150 = m(7) * t113 + t111;
t149 = t10 ^ 2;
t148 = t35 ^ 2;
t147 = 0.2e1 * t10;
t146 = 0.2e1 * t42;
t137 = Ifges(7,4) * t88;
t136 = Ifges(7,4) * t92;
t135 = t10 * t35;
t134 = t35 * t88;
t132 = t48 * mrSges(6,2);
t129 = t89 * mrSges(6,2);
t127 = Ifges(7,5) * t133 + Ifges(7,3) * t109;
t126 = Ifges(7,5) * t88 + Ifges(7,6) * t92;
t125 = t83 + t85;
t124 = t91 ^ 2 + t156;
t123 = t60 ^ 2 + t62 ^ 2;
t122 = pkin(4) * t129;
t66 = Ifges(7,2) * t92 + t137;
t67 = Ifges(7,1) * t88 + t136;
t121 = t92 * t66 + t88 * t67 + Ifges(6,3);
t120 = m(4) * t124;
t46 = pkin(10) + t48;
t117 = t125 * t46;
t116 = t125 * t74;
t115 = t124 * mrSges(4,3);
t114 = Ifges(5,3) + t121;
t110 = t37 * t62 - t38 * t60;
t108 = t60 * t90 - t62 * t94;
t107 = -t159 * t35 + (-mrSges(6,2) + t138 + t139) * t109;
t106 = (mrSges(5,1) * t94 - mrSges(5,2) * t90) * pkin(3);
t45 = -pkin(5) - t47;
t39 = t45 * t65;
t40 = t46 * t139;
t41 = t46 * t138;
t43 = t47 * mrSges(6,1);
t105 = t121 + t39 + t40 + t41 + t43 - t132;
t104 = t62 * mrSges(5,1) + t60 * mrSges(5,2) + t107;
t6 = Ifges(7,6) * t109 + (-Ifges(7,2) * t88 + t136) * t35;
t7 = Ifges(7,5) * t109 + (Ifges(7,1) * t92 - t137) * t35;
t102 = -t12 * mrSges(6,2) + mrSges(7,3) * t141 - t131 * t2 - t66 * t134 / 0.2e1 + t67 * t133 / 0.2e1 + Ifges(6,5) * t35 + t88 * t7 / 0.2e1 + t92 * t6 / 0.2e1 + (t126 / 0.2e1 - Ifges(6,6)) * t109 + t159 * t10;
t101 = t37 * mrSges(5,1) - t38 * mrSges(5,2) + Ifges(5,5) * t62 + Ifges(5,6) * t60 + t102;
t97 = qJ(2) ^ 2;
t30 = t109 ^ 2;
t1 = [Ifges(4,1) * t156 + 0.2e1 * t71 * (-t60 * mrSges(5,1) + t62 * mrSges(5,2)) + t62 * (Ifges(5,1) * t62 + Ifges(5,4) * t60) + t60 * (Ifges(5,4) * t62 + Ifges(5,2) * t60) + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + t14 * t147 + Ifges(3,1) + Ifges(2,3) - (2 * pkin(1) * mrSges(3,2)) + (mrSges(6,1) * t146 - 0.2e1 * t12 * mrSges(6,3) + Ifges(6,2) * t109 + t127) * t109 + (mrSges(6,2) * t146 + mrSges(6,3) * t147 + Ifges(6,1) * t35 - t88 * t6 + t92 * t7 + (-Ifges(7,6) * t88 - (2 * Ifges(6,4))) * t109) * t35 + m(4) * (t124 * t96 ^ 2 + t97) + m(3) * ((pkin(1) ^ 2) + t97) + m(5) * (t37 ^ 2 + t38 ^ 2 + t71 ^ 2) + m(6) * (t12 ^ 2 + t42 ^ 2 + t149) + m(7) * (t2 ^ 2 + t3 ^ 2 + t149) + (-0.2e1 * Ifges(4,4) * t95 + Ifges(4,2) * t91) * t91 + 0.2e1 * (mrSges(4,1) * t91 + mrSges(4,2) * t95 + mrSges(3,3)) * qJ(2) - 0.2e1 * t110 * mrSges(5,3) - 0.2e1 * t96 * t115; -m(3) * pkin(1) + mrSges(3,2) - (t35 * mrSges(6,3) + t14) * t35 - t123 * mrSges(5,3) - t115 + (-mrSges(6,3) * t109 + t111) * t109 + m(7) * (t109 * t113 - t135) + m(6) * (t109 * t12 - t135) + m(5) * t110 + t96 * t120; m(3) + m(6) * (t30 + t148) + m(7) * (t125 * t30 + t148) + m(5) * t123 + t120; t111 * t46 + m(6) * (-t10 * t47 + t12 * t48) + (mrSges(4,1) * t96 + Ifges(4,5)) * t95 + (-mrSges(4,2) * t96 - Ifges(4,6)) * t91 - t157 * mrSges(6,3) + (m(5) * (t37 * t94 + t38 * t90) + t108 * mrSges(5,3)) * pkin(3) + m(7) * (t10 * t45 + t113 * t46) + t101 + t45 * t14; t95 * mrSges(4,1) - t91 * mrSges(4,2) + m(6) * t157 + m(7) * (t109 * t117 - t35 * t45) - m(5) * t108 * pkin(3) + t104; -0.2e1 * t132 + Ifges(4,3) + 0.2e1 * t39 + 0.2e1 * t40 + 0.2e1 * t41 + 0.2e1 * t43 + 0.2e1 * t106 + m(7) * (t125 * t46 ^ 2 + t45 ^ 2) + m(6) * (t47 ^ 2 + t48 ^ 2) + m(5) * (t90 ^ 2 + t94 ^ 2) * pkin(3) ^ 2 + t114; (m(6) * (-t10 * t93 + t12 * t89) - t158 * mrSges(6,3)) * pkin(4) + t101 + t152 * t75 + t150 * t74; m(7) * (t109 * t116 - t35 * t75) + m(6) * t158 * pkin(4) + t104; m(7) * (t116 * t46 + t45 * t75) + (m(6) * (t47 * t93 + t48 * t89) - t129) * pkin(4) + t106 + t105 + Ifges(5,3) + t151; -0.2e1 * t122 + 0.2e1 * t51 + 0.2e1 * t63 + 0.2e1 * t64 + 0.2e1 * t79 + m(7) * (t125 * t74 ^ 2 + t75 ^ 2) + m(6) * (t89 ^ 2 + t93 ^ 2) * pkin(4) ^ 2 + t114; -t152 * pkin(5) + t150 * pkin(10) + t102; m(7) * (pkin(10) * t109 * t125 + pkin(5) * t35) + t107; m(7) * (-pkin(5) * t45 + pkin(10) * t117) + t105 + t153; m(7) * (-pkin(5) * t75 + pkin(10) * t116) - t122 + t121 + t151 + t153; -0.2e1 * t142 + m(7) * (pkin(10) ^ 2 * t125 + pkin(5) ^ 2) + 0.2e1 * t78 + 0.2e1 * t77 + t121; mrSges(7,1) * t2 - mrSges(7,2) * t3 - Ifges(7,6) * t134 + t127; -t112 * t109; -t112 * t46 + t126; -t112 * t74 + t126; -pkin(10) * t112 + t126; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
