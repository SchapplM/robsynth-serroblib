% Calculate joint inertia matrix for
% S6RPRRRR1
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
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:54:22
% EndTime: 2019-03-09 06:54:23
% DurationCPUTime: 0.79s
% Computational Cost: add. (1935->214), mult. (3528->320), div. (0->0), fcn. (3783->10), ass. (0->106)
t91 = cos(qJ(3));
t147 = t91 ^ 2;
t84 = sin(qJ(6));
t122 = t84 * mrSges(7,3);
t86 = sin(qJ(4));
t87 = sin(qJ(3));
t90 = cos(qJ(4));
t56 = -t86 * t87 + t90 * t91;
t57 = t86 * t91 + t87 * t90;
t85 = sin(qJ(5));
t89 = cos(qJ(5));
t34 = -t89 * t56 + t57 * t85;
t36 = t56 * t85 + t57 * t89;
t15 = -mrSges(7,2) * t34 - t122 * t36;
t88 = cos(qJ(6));
t124 = t36 * t88;
t16 = mrSges(7,1) * t34 - mrSges(7,3) * t124;
t103 = t88 * t15 - t84 * t16;
t61 = -mrSges(7,1) * t88 + mrSges(7,2) * t84;
t134 = pkin(5) * t61;
t78 = t84 ^ 2;
t130 = mrSges(7,3) * t78;
t73 = pkin(10) * t130;
t80 = t88 ^ 2;
t129 = mrSges(7,3) * t80;
t74 = pkin(10) * t129;
t146 = t73 + t74 - t134;
t82 = sin(pkin(11));
t68 = pkin(1) * t82 + pkin(7);
t131 = pkin(8) + t68;
t111 = t131 * t91;
t112 = t131 * t87;
t28 = t90 * t111 - t86 * t112;
t18 = pkin(9) * t56 + t28;
t27 = -t111 * t86 - t112 * t90;
t98 = -t57 * pkin(9) + t27;
t10 = t18 * t85 - t89 * t98;
t104 = mrSges(7,1) * t84 + mrSges(7,2) * t88;
t14 = t104 * t36;
t145 = m(7) * t10 + t14;
t136 = pkin(4) * t89;
t71 = -pkin(5) - t136;
t50 = t71 * t61;
t70 = pkin(4) * t85 + pkin(10);
t59 = t70 * t130;
t60 = t70 * t129;
t75 = mrSges(6,1) * t136;
t144 = t50 + t59 + t60 + t75;
t12 = t89 * t18 + t85 * t98;
t133 = pkin(10) * t36;
t135 = pkin(5) * t34;
t83 = cos(pkin(11));
t69 = -pkin(1) * t83 - pkin(2);
t58 = -pkin(3) * t91 + t69;
t38 = -pkin(4) * t56 + t58;
t13 = -t133 + t38 + t135;
t3 = t12 * t88 + t13 * t84;
t132 = t3 * t88;
t2 = -t12 * t84 + t13 * t88;
t106 = -t2 * t84 + t132;
t143 = m(7) * t106 + t103;
t142 = t10 ^ 2;
t141 = t34 ^ 2;
t140 = t57 ^ 2;
t139 = 0.2e1 * t10;
t138 = 0.2e1 * t38;
t137 = pkin(3) * t86;
t128 = Ifges(7,4) * t84;
t127 = Ifges(7,4) * t88;
t126 = t10 * t34;
t125 = t36 * t84;
t72 = pkin(3) * t90 + pkin(4);
t48 = t89 * t137 + t85 * t72;
t123 = t48 * mrSges(6,2);
t120 = t85 * mrSges(6,2);
t118 = Ifges(7,5) * t124 + Ifges(7,3) * t34;
t117 = Ifges(7,5) * t84 + Ifges(7,6) * t88;
t116 = t78 + t80;
t115 = t87 ^ 2 + t147;
t114 = pkin(4) * t120;
t62 = Ifges(7,2) * t88 + t128;
t63 = Ifges(7,1) * t84 + t127;
t113 = t88 * t62 + t84 * t63 + Ifges(6,3);
t46 = pkin(10) + t48;
t110 = t116 * t46;
t109 = t116 * t70;
t108 = -t56 * mrSges(5,1) + t57 * mrSges(5,2);
t107 = Ifges(5,3) + t113;
t105 = -t91 * mrSges(4,1) + t87 * mrSges(4,2);
t47 = -t137 * t85 + t72 * t89;
t29 = t34 * mrSges(6,1);
t102 = t34 * t61 - t29 + (-mrSges(6,2) + t129 + t130) * t36;
t101 = (mrSges(5,1) * t90 - mrSges(5,2) * t86) * pkin(3);
t45 = -pkin(5) - t47;
t37 = t45 * t61;
t39 = t46 * t130;
t40 = t46 * t129;
t41 = t47 * mrSges(6,1);
t100 = t113 + t37 + t39 + t40 + t41 - t123;
t99 = t102 - t108;
t8 = Ifges(7,6) * t34 + (-Ifges(7,2) * t84 + t127) * t36;
t9 = Ifges(7,5) * t34 + (Ifges(7,1) * t88 - t128) * t36;
t97 = -t12 * mrSges(6,2) + mrSges(7,3) * t132 - t122 * t2 - t62 * t125 / 0.2e1 + t63 * t124 / 0.2e1 + Ifges(6,5) * t36 + t84 * t9 / 0.2e1 + t88 * t8 / 0.2e1 + (t117 / 0.2e1 - Ifges(6,6)) * t34 + (t61 - mrSges(6,1)) * t10;
t96 = t27 * mrSges(5,1) - t28 * mrSges(5,2) + Ifges(5,5) * t57 + Ifges(5,6) * t56 + t97;
t33 = t36 ^ 2;
t1 = [-0.2e1 * t27 * t57 * mrSges(5,3) + 0.2e1 * t69 * t105 + Ifges(4,2) * t147 + Ifges(5,1) * t140 + 0.2e1 * t58 * t108 + t29 * t138 + t14 * t139 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + Ifges(2,3) + Ifges(3,3) + (0.2e1 * t28 * mrSges(5,3) + 0.2e1 * Ifges(5,4) * t57 + Ifges(5,2) * t56) * t56 + (mrSges(6,2) * t138 + mrSges(6,3) * t139 + Ifges(6,1) * t36 - t84 * t8 + t88 * t9) * t36 + (-0.2e1 * t12 * mrSges(6,3) + Ifges(6,2) * t34 + (-Ifges(7,6) * t84 - (2 * Ifges(6,4))) * t36 + t118) * t34 + m(4) * (t115 * t68 ^ 2 + t69 ^ 2) + m(5) * (t27 ^ 2 + t28 ^ 2 + t58 ^ 2) + m(6) * (t12 ^ 2 + t38 ^ 2 + t142) + m(7) * (t2 ^ 2 + t3 ^ 2 + t142) + m(3) * (t82 ^ 2 + t83 ^ 2) * pkin(1) ^ 2 + (Ifges(4,1) * t87 + 0.2e1 * Ifges(4,4) * t91) * t87 + 0.2e1 * (t83 * mrSges(3,1) - t82 * mrSges(3,2)) * pkin(1) + 0.2e1 * t115 * t68 * mrSges(4,3); t34 * t14 + t103 * t36 + m(7) * (t106 * t36 + t126) + m(6) * (t12 * t36 + t126) + m(5) * (t27 * t56 + t28 * t57); m(3) + m(7) * (t116 * t33 + t141) + m(6) * (t33 + t141) + m(5) * (t56 ^ 2 + t140) + m(4) * t115; t103 * t46 + m(6) * (-t10 * t47 + t12 * t48) + (-mrSges(4,1) * t87 - mrSges(4,2) * t91) * t68 + (-t34 * t48 - t36 * t47) * mrSges(6,3) + (m(5) * (t27 * t90 + t28 * t86) + (t56 * t86 - t57 * t90) * mrSges(5,3)) * pkin(3) + m(7) * (t10 * t45 + t106 * t46) + t96 + Ifges(4,6) * t91 + Ifges(4,5) * t87 + t45 * t14; m(7) * (t110 * t36 + t34 * t45) + m(6) * (-t34 * t47 + t36 * t48) + m(5) * (t56 * t90 + t57 * t86) * pkin(3) + t99 - t105; -0.2e1 * t123 + Ifges(4,3) + 0.2e1 * t37 + 0.2e1 * t39 + 0.2e1 * t40 + 0.2e1 * t41 + 0.2e1 * t101 + m(7) * (t116 * t46 ^ 2 + t45 ^ 2) + m(6) * (t47 ^ 2 + t48 ^ 2) + m(5) * (t86 ^ 2 + t90 ^ 2) * pkin(3) ^ 2 + t107; (m(6) * (-t10 * t89 + t12 * t85) + (-t34 * t85 - t36 * t89) * mrSges(6,3)) * pkin(4) + t96 + t145 * t71 + t143 * t70; m(7) * (t109 * t36 + t71 * t34) + m(6) * (-t34 * t89 + t36 * t85) * pkin(4) + t99; t101 + m(7) * (t109 * t46 + t45 * t71) + (m(6) * (t47 * t89 + t48 * t85) - t120) * pkin(4) + t100 + Ifges(5,3) + t144; -0.2e1 * t114 + 0.2e1 * t50 + 0.2e1 * t59 + 0.2e1 * t60 + 0.2e1 * t75 + m(7) * (t116 * t70 ^ 2 + t71 ^ 2) + m(6) * (t85 ^ 2 + t89 ^ 2) * pkin(4) ^ 2 + t107; -t145 * pkin(5) + t143 * pkin(10) + t97; m(7) * (t116 * t133 - t135) + t102; m(7) * (-pkin(5) * t45 + pkin(10) * t110) + t100 + t146; -t114 + m(7) * (-pkin(5) * t71 + pkin(10) * t109) + t113 + t144 + t146; 0.2e1 * t73 + 0.2e1 * t74 - 0.2e1 * t134 + m(7) * (pkin(10) ^ 2 * t116 + pkin(5) ^ 2) + t113; mrSges(7,1) * t2 - mrSges(7,2) * t3 - Ifges(7,6) * t125 + t118; -t14; -t104 * t46 + t117; -t104 * t70 + t117; -pkin(10) * t104 + t117; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
