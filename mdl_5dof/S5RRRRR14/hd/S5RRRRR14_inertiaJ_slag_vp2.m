% Calculate joint inertia matrix for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR14_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR14_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR14_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR14_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:42:21
% EndTime: 2024-09-27 18:42:22
% DurationCPUTime: 0.19s
% Computational Cost: add. (1673->196), mult. (3842->262), div. (0->0), fcn. (3889->10), ass. (0->102)
t90 = sin(pkin(5));
t93 = sin(qJ(4));
t94 = sin(qJ(3));
t97 = cos(qJ(4));
t98 = cos(qJ(3));
t59 = (-t93 * t94 + t97 * t98) * t90;
t60 = (t93 * t98 + t94 * t97) * t90;
t144 = Ifges(5,5) * t60 + Ifges(5,6) * t59;
t92 = sin(qJ(5));
t96 = cos(qJ(5));
t32 = t59 * t96 - t60 * t92;
t33 = t59 * t92 + t60 * t96;
t143 = Ifges(6,5) * t33 + Ifges(6,6) * t32;
t121 = t90 * t98;
t122 = t90 * t94;
t142 = Ifges(4,5) * t122 + Ifges(4,6) * t121;
t89 = t90 ^ 2;
t125 = pkin(4) * t96;
t91 = cos(pkin(5));
t27 = -mrSges(6,2) * t91 + mrSges(6,3) * t32;
t28 = mrSges(6,1) * t91 - mrSges(6,3) * t33;
t141 = t92 * pkin(4) * t27 + t28 * t125;
t112 = Ifges(4,3) * t91 + t142;
t128 = pkin(3) * t97;
t129 = pkin(3) * t93;
t49 = -mrSges(5,2) * t91 + mrSges(5,3) * t59;
t50 = mrSges(5,1) * t91 - mrSges(5,3) * t60;
t81 = pkin(4) + t128;
t63 = -t92 * t129 + t81 * t96;
t64 = t96 * t129 + t81 * t92;
t140 = t50 * t128 + t49 * t129 + t64 * t27 + t63 * t28 + t112;
t12 = -mrSges(6,1) * t32 + mrSges(6,2) * t33;
t139 = 0.2e1 * t12;
t138 = 0.2e1 * t27;
t137 = 0.2e1 * t28;
t34 = -mrSges(5,1) * t59 + mrSges(5,2) * t60;
t136 = 0.2e1 * t34;
t135 = 0.2e1 * t49;
t134 = 0.2e1 * t50;
t68 = mrSges(4,1) * t91 - mrSges(4,3) * t122;
t133 = 0.2e1 * t68;
t69 = -mrSges(4,2) * t91 + mrSges(4,3) * t121;
t132 = 0.2e1 * t69;
t131 = m(5) * pkin(3);
t130 = m(6) * pkin(4);
t127 = pkin(3) * t98;
t126 = pkin(4) * t59;
t124 = t64 * mrSges(6,2);
t123 = t89 * (-mrSges(4,1) * t98 + mrSges(4,2) * t94);
t120 = t91 * t94;
t119 = t91 * t98;
t118 = t92 * mrSges(6,2);
t117 = Ifges(5,3) + Ifges(6,3);
t99 = cos(qJ(2));
t82 = pkin(1) * t99 + pkin(2);
t73 = t82 * t119;
t95 = sin(qJ(2));
t74 = pkin(1) * t95 + pkin(8) * t90;
t88 = t91 * pkin(3);
t38 = t73 + t88 + (-pkin(9) * t90 - t74) * t94;
t52 = t82 * t120 + t98 * t74;
t77 = pkin(9) * t121;
t40 = t77 + t52;
t16 = t93 * t38 + t97 * t40;
t80 = pkin(2) * t119;
t46 = t80 + t88 + (-pkin(8) - pkin(9)) * t122;
t67 = pkin(2) * t120 + pkin(8) * t121;
t55 = t77 + t67;
t24 = t93 * t46 + t97 * t55;
t116 = -0.2e1 * t123;
t115 = pkin(4) * t118;
t114 = Ifges(6,3) * t91 + t143;
t113 = Ifges(5,3) * t91 + t144;
t111 = t91 * pkin(4) - pkin(10) * t60;
t15 = t97 * t38 - t40 * t93;
t23 = t97 * t46 - t55 * t93;
t61 = t63 * mrSges(6,1);
t110 = Ifges(6,3) + t61 - t124;
t71 = (-pkin(2) - t127) * t90;
t62 = (-t82 - t127) * t90;
t58 = t59 * pkin(10);
t11 = t58 + t16;
t9 = t111 + t15;
t2 = -t11 * t92 + t9 * t96;
t3 = t11 * t96 + t9 * t92;
t109 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + t114;
t13 = t111 + t23;
t17 = t58 + t24;
t5 = t13 * t96 - t17 * t92;
t6 = t13 * t92 + t17 * t96;
t108 = t5 * mrSges(6,1) - t6 * mrSges(6,2) + t114;
t107 = (mrSges(3,1) * t99 - mrSges(3,2) * t95) * pkin(1);
t106 = (mrSges(5,1) * t97 - mrSges(5,2) * t93) * pkin(3);
t105 = Ifges(5,1) * t60 ^ 2 + Ifges(6,1) * t33 ^ 2 + Ifges(3,3) + (0.2e1 * Ifges(5,4) * t60 + Ifges(5,2) * t59) * t59 + (0.2e1 * Ifges(6,4) * t33 + Ifges(6,2) * t32) * t32 + ((Ifges(4,1) * t94 + Ifges(4,4) * t98) * t122 + (Ifges(4,4) * t94 + Ifges(4,2) * t98) * t121) * t90 + (t112 + t113 + t114 + t142 + t143 + t144) * t91;
t104 = t15 * mrSges(5,1) - t16 * mrSges(5,2) + t109 + t113;
t103 = t23 * mrSges(5,1) - t24 * mrSges(5,2) + t108 + t113;
t83 = mrSges(6,1) * t125;
t66 = -pkin(8) * t122 + t80;
t51 = -t74 * t94 + t73;
t41 = t71 - t126;
t37 = t62 - t126;
t1 = [0.2e1 * t107 + Ifges(2,3) + m(3) * (t95 ^ 2 + t99 ^ 2) * pkin(1) ^ 2 + m(4) * (t82 ^ 2 * t89 + t51 ^ 2 + t52 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2 + t62 ^ 2) + m(6) * (t2 ^ 2 + t3 ^ 2 + t37 ^ 2) + t82 * t116 + t62 * t136 + t51 * t133 + t52 * t132 + t16 * t135 + t15 * t134 + t37 * t139 + t3 * t138 + t2 * t137 + t105; (-pkin(2) - t82) * t123 + m(4) * (pkin(2) * t82 * t89 + t51 * t66 + t52 * t67) + m(5) * (t15 * t23 + t16 * t24 + t62 * t71) + m(6) * (t2 * t5 + t3 * t6 + t37 * t41) + (t67 + t52) * t69 + (t66 + t51) * t68 + (t23 + t15) * t50 + (t24 + t16) * t49 + (t62 + t71) * t34 + (t5 + t2) * t28 + (t6 + t3) * t27 + (t37 + t41) * t12 + t107 + t105; m(6) * (t41 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t23 ^ 2 + t24 ^ 2 + t71 ^ 2) + m(4) * (pkin(2) ^ 2 * t89 + t66 ^ 2 + t67 ^ 2) + pkin(2) * t116 + t66 * t133 + t67 * t132 + t71 * t136 + t24 * t135 + t23 * t134 + t41 * t139 + t6 * t138 + t5 * t137 + t105; m(6) * (t2 * t63 + t3 * t64) + (t15 * t97 + t16 * t93) * t131 + t51 * mrSges(4,1) - t52 * mrSges(4,2) + t104 + t140; m(6) * (t5 * t63 + t6 * t64) + (t23 * t97 + t24 * t93) * t131 + t66 * mrSges(4,1) - t67 * mrSges(4,2) + t103 + t140; -0.2e1 * t124 + Ifges(4,3) + 0.2e1 * t61 + 0.2e1 * t106 + m(6) * (t63 ^ 2 + t64 ^ 2) + m(5) * (t93 ^ 2 + t97 ^ 2) * pkin(3) ^ 2 + t117; (t2 * t96 + t3 * t92) * t130 + t104 + t141; (t5 * t96 + t6 * t92) * t130 + t103 + t141; Ifges(5,3) + t83 + t106 + (m(6) * (t63 * t96 + t64 * t92) - t118) * pkin(4) + t110; -0.2e1 * t115 + 0.2e1 * t83 + m(6) * (t92 ^ 2 + t96 ^ 2) * pkin(4) ^ 2 + t117; t109; t108; t110; Ifges(6,3) + t83 - t115; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
