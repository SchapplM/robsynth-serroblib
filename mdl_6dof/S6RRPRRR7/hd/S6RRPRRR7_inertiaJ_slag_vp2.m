% Calculate joint inertia matrix for
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:55:51
% EndTime: 2019-03-09 13:55:54
% DurationCPUTime: 1.51s
% Computational Cost: add. (1963->297), mult. (3519->405), div. (0->0), fcn. (3595->8), ass. (0->115)
t92 = sin(qJ(5));
t96 = cos(qJ(5));
t67 = -t96 * mrSges(6,1) + mrSges(6,2) * t92;
t129 = t67 - mrSges(5,1);
t127 = t92 ^ 2 + t96 ^ 2;
t154 = t127 * mrSges(6,3);
t91 = sin(qJ(6));
t95 = cos(qJ(6));
t56 = t91 * t92 - t95 * t96;
t59 = t91 * t96 + t92 * t95;
t29 = mrSges(7,1) * t56 + mrSges(7,2) * t59;
t93 = sin(qJ(4));
t45 = t59 * t93;
t46 = t56 * t93;
t97 = cos(qJ(4));
t170 = -(t29 + t129) * t97 - (mrSges(5,2) - t154) * t93 - (-t45 * t59 - t46 * t56) * mrSges(7,3);
t33 = Ifges(7,1) * t59 - Ifges(7,4) * t56;
t167 = -t33 / 0.2e1;
t31 = Ifges(7,4) * t59 - Ifges(7,2) * t56;
t168 = t31 / 0.2e1;
t94 = sin(qJ(2));
t98 = cos(qJ(2));
t58 = -t93 * t98 + t94 * t97;
t132 = t58 * t96;
t57 = t93 * t94 + t97 * t98;
t66 = -t98 * pkin(2) - t94 * qJ(3) - pkin(1);
t49 = t98 * pkin(3) - t66;
t20 = pkin(4) * t57 - pkin(9) * t58 + t49;
t142 = pkin(7) - pkin(8);
t123 = t142 * t94;
t73 = t142 * t98;
t41 = t93 * t123 + t97 * t73;
t9 = t96 * t20 - t41 * t92;
t4 = pkin(5) * t57 - pkin(10) * t132 + t9;
t10 = t92 * t20 + t96 * t41;
t133 = t58 * t92;
t7 = -pkin(10) * t133 + t10;
t2 = t4 * t95 - t7 * t91;
t38 = -t97 * t123 + t73 * t93;
t21 = pkin(5) * t133 + t38;
t24 = t59 * t58;
t25 = t56 * t58;
t3 = t4 * t91 + t7 * t95;
t5 = -Ifges(7,4) * t25 - Ifges(7,2) * t24 + Ifges(7,6) * t57;
t6 = -Ifges(7,1) * t25 - Ifges(7,4) * t24 + Ifges(7,5) * t57;
t169 = t25 * t167 - t24 * t168 + t21 * t29 - (t2 * t59 + t3 * t56) * mrSges(7,3) + t129 * t38 - t41 * mrSges(5,2) - t56 * t5 / 0.2e1 + t59 * t6 / 0.2e1;
t140 = pkin(5) * t96;
t99 = -pkin(2) - pkin(3);
t64 = -t93 * qJ(3) + t97 * t99;
t62 = pkin(4) - t64;
t48 = t62 + t140;
t166 = t48 * t29;
t75 = -pkin(4) - t140;
t165 = t75 * t29;
t164 = t94 ^ 2 + t98 ^ 2;
t138 = Ifges(6,4) * t92;
t69 = Ifges(6,2) * t96 + t138;
t137 = Ifges(6,4) * t96;
t70 = Ifges(6,1) * t92 + t137;
t107 = t96 * t69 + t92 * t70 + Ifges(5,3);
t160 = -t56 * t31 + t59 * t33 + t107;
t157 = Ifges(6,5) * t132 + Ifges(6,3) * t57;
t156 = (t56 * t91 + t59 * t95) * mrSges(7,3);
t135 = Ifges(7,6) * t56;
t136 = Ifges(7,5) * t59;
t153 = t135 - t136;
t152 = -m(4) * pkin(2) - mrSges(4,1);
t68 = Ifges(6,5) * t92 + Ifges(6,6) * t96;
t151 = t136 / 0.2e1 - t135 / 0.2e1 - Ifges(5,6) + t68 / 0.2e1;
t149 = t38 ^ 2;
t148 = 0.2e1 * t38;
t147 = 0.2e1 * t49;
t146 = -0.2e1 * t66;
t145 = -0.2e1 * t67;
t141 = -pkin(10) - pkin(9);
t65 = t97 * qJ(3) + t93 * t99;
t63 = -pkin(9) + t65;
t139 = pkin(10) - t63;
t134 = t38 * t97;
t131 = t64 * mrSges(5,1);
t130 = t65 * mrSges(5,2);
t128 = t164 * pkin(7) ^ 2;
t125 = 0.2e1 * mrSges(7,3);
t124 = -Ifges(7,5) * t25 - Ifges(7,6) * t24 + Ifges(7,3) * t57;
t120 = t127 * t63;
t119 = t127 * t93;
t118 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t15 = Ifges(6,5) * t57 + (Ifges(6,1) * t96 - t138) * t58;
t117 = t9 * mrSges(6,3) - t15 / 0.2e1;
t14 = Ifges(6,6) * t57 + (-Ifges(6,2) * t92 + t137) * t58;
t116 = -t14 / 0.2e1 - t10 * mrSges(6,3);
t113 = t10 * t96 - t9 * t92;
t112 = mrSges(6,1) * t92 + mrSges(6,2) * t96;
t42 = t139 * t92;
t43 = t139 * t96;
t16 = t42 * t95 + t43 * t91;
t17 = t42 * t91 - t43 * t95;
t109 = t16 * mrSges(7,1) - t17 * mrSges(7,2) + t153;
t71 = t141 * t92;
t72 = t141 * t96;
t37 = t71 * t95 + t72 * t91;
t40 = t71 * t91 - t72 * t95;
t108 = t37 * mrSges(7,1) - t40 * mrSges(7,2) - t153;
t106 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t124;
t105 = (mrSges(7,1) * t95 - mrSges(7,2) * t91) * pkin(5);
t104 = Ifges(5,5) + t96 * t70 / 0.2e1 - t92 * t69 / 0.2e1;
t89 = t97 ^ 2;
t86 = t93 ^ 2;
t28 = mrSges(6,1) * t57 - mrSges(6,3) * t132;
t27 = -mrSges(6,2) * t57 - mrSges(6,3) * t133;
t26 = t112 * t58;
t12 = mrSges(7,1) * t57 + mrSges(7,3) * t25;
t11 = -mrSges(7,2) * t57 - mrSges(7,3) * t24;
t8 = mrSges(7,1) * t24 - mrSges(7,2) * t25;
t1 = [0.2e1 * t10 * t27 + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t12 + 0.2e1 * t21 * t8 - t24 * t5 - t25 * t6 + t26 * t148 + 0.2e1 * t9 * t28 + Ifges(2,3) + (0.2e1 * pkin(1) * mrSges(3,1) + mrSges(4,1) * t146 + (Ifges(4,3) + Ifges(3,2)) * t98) * t98 + (-0.2e1 * pkin(1) * mrSges(3,2) + mrSges(4,3) * t146 + (Ifges(3,1) + Ifges(4,1)) * t94 + 0.2e1 * (Ifges(3,4) - Ifges(4,5)) * t98) * t94 + (mrSges(5,2) * t147 + mrSges(5,3) * t148 + Ifges(5,1) * t58 - t92 * t14 + t96 * t15) * t58 + (mrSges(5,1) * t147 - 0.2e1 * t41 * mrSges(5,3) + Ifges(5,2) * t57 + (-Ifges(6,6) * t92 - (2 * Ifges(5,4))) * t58 + t124 + t157) * t57 + m(4) * (t66 ^ 2 + t128) + m(3) * (pkin(1) ^ 2 + t128) + m(5) * (t41 ^ 2 + t49 ^ 2 + t149) + m(6) * (t10 ^ 2 + t9 ^ 2 + t149) + m(7) * (t2 ^ 2 + t21 ^ 2 + t3 ^ 2) + 0.2e1 * (mrSges(4,2) + mrSges(3,3)) * pkin(7) * t164; (qJ(3) * mrSges(4,2) + Ifges(3,6) - Ifges(4,6)) * t98 + (-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t94 + m(5) * (-t38 * t64 + t41 * t65) + m(7) * (t16 * t2 + t17 * t3 + t21 * t48) + t62 * t26 + (-t65 * mrSges(5,3) - t151) * t57 + m(6) * (t113 * t63 + t38 * t62) + (-t64 * mrSges(5,3) - t104) * t58 + ((m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * t98 + (-mrSges(3,1) + t152) * t94) * pkin(7) + t48 * t8 + t16 * t12 + t17 * t11 + (t63 * t27 + t116) * t96 + (-t63 * t28 + t117) * t92 - t169; 0.2e1 * pkin(2) * mrSges(4,1) - 0.2e1 * t131 + 0.2e1 * t130 + 0.2e1 * qJ(3) * mrSges(4,3) - 0.2e1 * t166 + t62 * t145 + Ifges(4,2) + Ifges(3,3) + (t16 * t59 + t17 * t56) * t125 - 0.2e1 * t63 * t154 + m(7) * (t16 ^ 2 + t17 ^ 2 + t48 ^ 2) + m(6) * (t127 * t63 ^ 2 + t62 ^ 2) + m(5) * (t64 ^ 2 + t65 ^ 2) + m(4) * (pkin(2) ^ 2 + qJ(3) ^ 2) + t160; -t46 * t11 - t45 * t12 + (m(4) * pkin(7) + mrSges(4,2)) * t94 + (-t58 * mrSges(5,3) - t26 - t8) * t97 + (-t57 * mrSges(5,3) + t96 * t27 - t92 * t28) * t93 + m(7) * (-t2 * t45 - t21 * t97 - t3 * t46) + m(6) * (t113 * t93 - t134) + m(5) * (t41 * t93 - t134); m(7) * (-t16 * t45 - t17 * t46 - t48 * t97) + m(6) * (t119 * t63 - t62 * t97) + m(5) * (t64 * t97 + t65 * t93) + t152 - t170; m(4) + m(5) * (t86 + t89) + m(6) * (t127 * t86 + t89) + m(7) * (t45 ^ 2 + t46 ^ 2 + t89); t75 * t8 + t37 * t12 + t40 * t11 - pkin(4) * t26 + t151 * t57 + (pkin(9) * t27 - t116) * t96 + (-pkin(9) * t28 - t117) * t92 + m(6) * (-pkin(4) * t38 + pkin(9) * t113) + m(7) * (t2 * t37 + t21 * t75 + t3 * t40) + t104 * t58 + t169; t131 - t130 + t166 - t165 + (pkin(4) + t62) * t67 + 0.2e1 * t167 * t59 + 0.2e1 * t168 * t56 + m(6) * (-pkin(4) * t62 + pkin(9) * t120) + m(7) * (t16 * t37 + t17 * t40 + t48 * t75) + ((-t16 + t37) * t59 + (-t17 + t40) * t56) * mrSges(7,3) + (-pkin(9) * t127 + t120) * mrSges(6,3) - t107; m(6) * (pkin(4) * t97 + pkin(9) * t119) + m(7) * (-t37 * t45 - t40 * t46 - t75 * t97) + t170; pkin(4) * t145 + 0.2e1 * t165 + (-t37 * t59 - t40 * t56) * t125 + 0.2e1 * pkin(9) * t154 + m(7) * (t37 ^ 2 + t40 ^ 2 + t75 ^ 2) + m(6) * (pkin(9) ^ 2 * t127 + pkin(4) ^ 2) + t160; -Ifges(6,6) * t133 + t9 * mrSges(6,1) - t10 * mrSges(6,2) + (m(7) * (t2 * t95 + t3 * t91) + t91 * t11 + t95 * t12) * pkin(5) + t106 + t157; -t112 * t63 + (m(7) * (t16 * t95 + t17 * t91) + t156) * pkin(5) + t109 - t68; -t112 * t93 + m(7) * (-t45 * t95 - t46 * t91) * pkin(5) + t118; -t112 * pkin(9) + (m(7) * (t37 * t95 + t40 * t91) - t156) * pkin(5) + t108 + t68; Ifges(6,3) + Ifges(7,3) + m(7) * (t91 ^ 2 + t95 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t105; t106; t109; t118; t108; Ifges(7,3) + t105; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
