% Calculate joint inertia matrix for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:05:26
% EndTime: 2019-03-08 23:05:29
% DurationCPUTime: 1.04s
% Computational Cost: add. (1747->273), mult. (3641->396), div. (0->0), fcn. (4000->12), ass. (0->109)
t110 = sin(pkin(12));
t112 = cos(pkin(12));
t90 = -t112 * mrSges(6,1) + mrSges(6,2) * t110;
t162 = t90 - mrSges(5,1);
t111 = sin(pkin(6));
t121 = cos(qJ(2));
t140 = t111 * t121;
t115 = sin(qJ(4));
t119 = cos(qJ(4));
t113 = cos(pkin(6));
t116 = sin(qJ(3));
t120 = cos(qJ(3));
t117 = sin(qJ(2));
t141 = t111 * t117;
t71 = t113 * t120 - t116 * t141;
t72 = t113 * t116 + t120 * t141;
t42 = t115 * t71 + t119 * t72;
t26 = -t110 * t42 - t112 * t140;
t27 = -t110 * t140 + t112 * t42;
t131 = -t110 * t26 + t112 * t27;
t114 = sin(qJ(6));
t118 = cos(qJ(6));
t84 = t110 * t118 + t112 * t114;
t86 = t115 * t120 + t116 * t119;
t43 = t84 * t86;
t83 = -t110 * t114 + t112 * t118;
t44 = t83 * t86;
t13 = t43 * mrSges(7,1) + mrSges(7,2) * t44;
t143 = t112 * t86;
t146 = t110 * t86;
t50 = mrSges(6,1) * t146 + mrSges(6,2) * t143;
t161 = t13 + t50;
t40 = t115 * t72 - t119 * t71;
t39 = t40 ^ 2;
t156 = -pkin(9) - pkin(8);
t135 = t156 * t116;
t95 = t156 * t120;
t64 = -t115 * t95 - t119 * t135;
t160 = t64 ^ 2;
t54 = -t83 * mrSges(7,1) + mrSges(7,2) * t84;
t159 = 0.2e1 * t54;
t158 = 0.2e1 * t64;
t109 = t120 ^ 2;
t155 = pkin(3) * t119;
t154 = t40 * t64;
t153 = t83 * mrSges(7,3);
t152 = t84 * mrSges(7,3);
t102 = -pkin(3) * t120 - pkin(2);
t85 = t115 * t116 - t119 * t120;
t51 = pkin(4) * t85 - qJ(5) * t86 + t102;
t66 = t115 * t135 - t119 * t95;
t20 = t110 * t51 + t112 * t66;
t151 = Ifges(7,5) * t84 + Ifges(7,6) * t83;
t150 = Ifges(6,4) * t110;
t149 = Ifges(6,4) * t112;
t19 = -t110 * t66 + t112 * t51;
t148 = t110 * t19;
t145 = t112 * t20;
t98 = pkin(3) * t115 + qJ(5);
t142 = t112 * t98;
t139 = t110 ^ 2 + t112 ^ 2;
t138 = t116 ^ 2 + t109;
t137 = 0.2e1 * mrSges(7,3);
t136 = Ifges(7,5) * t44 - Ifges(7,6) * t43 + Ifges(7,3) * t85;
t99 = -pkin(5) * t112 - pkin(4);
t134 = t139 * qJ(5);
t55 = Ifges(7,4) * t84 + Ifges(7,2) * t83;
t56 = Ifges(7,1) * t84 + Ifges(7,4) * t83;
t92 = Ifges(6,2) * t112 + t150;
t93 = Ifges(6,1) * t110 + t149;
t133 = t110 * t93 + t112 * t92 + t55 * t83 + t56 * t84 + Ifges(5,3);
t132 = t145 - t148;
t130 = -t116 * t71 + t120 * t72;
t129 = 0.2e1 * t139 * mrSges(6,3);
t128 = (mrSges(5,1) * t119 - mrSges(5,2) * t115) * pkin(3);
t127 = t54 + t90;
t5 = -t114 * t27 + t118 * t26;
t6 = t114 * t26 + t118 * t27;
t126 = -t42 * mrSges(5,2) - t152 * t5 + t6 * t153 + t131 * mrSges(6,3) + (t54 + t162) * t40;
t10 = Ifges(7,4) * t44 - Ifges(7,2) * t43 + Ifges(7,6) * t85;
t11 = Ifges(7,1) * t44 - Ifges(7,4) * t43 + Ifges(7,5) * t85;
t12 = -pkin(10) * t146 + t20;
t9 = pkin(5) * t85 - pkin(10) * t143 + t19;
t2 = -t114 * t12 + t118 * t9;
t29 = Ifges(6,6) * t85 + (-Ifges(6,2) * t110 + t149) * t86;
t3 = t114 * t9 + t118 * t12;
t30 = Ifges(6,5) * t85 + (Ifges(6,1) * t112 - t150) * t86;
t34 = pkin(5) * t146 + t64;
t125 = -t66 * mrSges(5,2) - t152 * t2 + t3 * t153 + t34 * t54 + mrSges(6,3) * t145 - t43 * t55 / 0.2e1 + t44 * t56 / 0.2e1 + t110 * t30 / 0.2e1 + t112 * t29 / 0.2e1 - t92 * t146 / 0.2e1 + t93 * t143 / 0.2e1 + t83 * t10 / 0.2e1 - Ifges(5,6) * t85 + t84 * t11 / 0.2e1 + Ifges(5,5) * t86 + t162 * t64 + (Ifges(6,5) * t110 + Ifges(6,6) * t112 + t151) * t85 / 0.2e1;
t106 = t111 ^ 2;
t104 = t112 * pkin(10);
t101 = -pkin(4) - t155;
t97 = t106 * t121 ^ 2;
t94 = -mrSges(4,1) * t120 + mrSges(4,2) * t116;
t91 = qJ(5) * t112 + t104;
t89 = (-pkin(10) - qJ(5)) * t110;
t88 = t99 - t155;
t74 = t104 + t142;
t73 = (-pkin(10) - t98) * t110;
t62 = t114 * t89 + t118 * t91;
t61 = -t114 * t91 + t118 * t89;
t57 = mrSges(5,1) * t85 + mrSges(5,2) * t86;
t53 = mrSges(6,1) * t85 - mrSges(6,3) * t143;
t52 = -mrSges(6,2) * t85 - mrSges(6,3) * t146;
t48 = t114 * t73 + t118 * t74;
t47 = -t114 * t74 + t118 * t73;
t23 = mrSges(7,1) * t85 - mrSges(7,3) * t44;
t22 = -mrSges(7,2) * t85 - mrSges(7,3) * t43;
t1 = [m(2) + m(6) * (t26 ^ 2 + t27 ^ 2 + t39) + m(7) * (t5 ^ 2 + t6 ^ 2 + t39) + m(5) * (t42 ^ 2 + t39 + t97) + m(4) * (t71 ^ 2 + t72 ^ 2 + t97) + m(3) * (t106 * t117 ^ 2 + t113 ^ 2 + t97); -t42 * t85 * mrSges(5,3) + t6 * t22 + t5 * t23 + t26 * t53 + t27 * t52 + t130 * mrSges(4,3) + (t86 * mrSges(5,3) + t161) * t40 + (-t117 * mrSges(3,2) + (mrSges(3,1) - t57 - t94) * t121) * t111 + m(6) * (t19 * t26 + t20 * t27 + t154) + m(7) * (t2 * t5 + t3 * t6 + t34 * t40) + m(5) * (-t102 * t140 + t42 * t66 + t154) + m(4) * (pkin(2) * t140 + pkin(8) * t130); Ifges(4,2) * t109 - 0.2e1 * pkin(2) * t94 - t43 * t10 + 0.2e1 * t102 * t57 + t44 * t11 + 0.2e1 * t34 * t13 + 0.2e1 * t19 * t53 + 0.2e1 * t2 * t23 + 0.2e1 * t20 * t52 + 0.2e1 * t3 * t22 + t50 * t158 + Ifges(3,3) + (Ifges(4,1) * t116 + 0.2e1 * Ifges(4,4) * t120) * t116 + 0.2e1 * t138 * pkin(8) * mrSges(4,3) + (mrSges(5,3) * t158 + Ifges(5,1) * t86 - t110 * t29 + t112 * t30) * t86 + m(4) * (pkin(8) ^ 2 * t138 + pkin(2) ^ 2) + m(5) * (t102 ^ 2 + t66 ^ 2 + t160) + m(6) * (t19 ^ 2 + t20 ^ 2 + t160) + m(7) * (t2 ^ 2 + t3 ^ 2 + t34 ^ 2) + (-0.2e1 * t66 * mrSges(5,3) + (Ifges(6,3) + Ifges(5,2)) * t85 + (Ifges(6,5) * t112 - Ifges(6,6) * t110 - (2 * Ifges(5,4))) * t86 + t136) * t85; t71 * mrSges(4,1) - t72 * mrSges(4,2) + m(6) * (t101 * t40 + t131 * t98) + m(7) * (t40 * t88 + t47 * t5 + t48 * t6) + m(5) * (t115 * t42 - t119 * t40) * pkin(3) + t126; t125 + (m(5) * (t115 * t66 - t119 * t64) + (-t115 * t85 - t119 * t86) * mrSges(5,3)) * pkin(3) + m(6) * (t101 * t64 + t132 * t98) + (-t19 * mrSges(6,3) - t98 * t53) * t110 + t52 * t142 + Ifges(4,6) * t120 + Ifges(4,5) * t116 + t88 * t13 + t101 * t50 + m(7) * (t2 * t47 + t3 * t48 + t34 * t88) + t47 * t23 + t48 * t22 + (-mrSges(4,1) * t116 - mrSges(4,2) * t120) * pkin(8); 0.2e1 * t101 * t90 + t88 * t159 + Ifges(4,3) + 0.2e1 * t128 + (-t47 * t84 + t48 * t83) * t137 + t98 * t129 + m(7) * (t47 ^ 2 + t48 ^ 2 + t88 ^ 2) + m(6) * (t139 * t98 ^ 2 + t101 ^ 2) + m(5) * (t115 ^ 2 + t119 ^ 2) * pkin(3) ^ 2 + t133; m(6) * (-pkin(4) * t40 + qJ(5) * t131) + m(7) * (t40 * t99 + t5 * t61 + t6 * t62) + t126; t125 + m(6) * (-pkin(4) * t64 + qJ(5) * t132) + t99 * t13 + t61 * t23 + t62 * t22 - pkin(4) * t50 - mrSges(6,3) * t148 + m(7) * (t2 * t61 + t3 * t62 + t34 * t99) + (-t110 * t53 + t112 * t52) * qJ(5); (t101 - pkin(4)) * t90 + (t88 + t99) * t54 + t128 + m(7) * (t47 * t61 + t48 * t62 + t88 * t99) + m(6) * (-pkin(4) * t101 + t134 * t98) + ((-t47 - t61) * t84 + (t48 + t62) * t83) * mrSges(7,3) + (t139 * t98 + t134) * mrSges(6,3) + t133; -0.2e1 * pkin(4) * t90 + t99 * t159 + (-t61 * t84 + t62 * t83) * t137 + qJ(5) * t129 + m(7) * (t61 ^ 2 + t62 ^ 2 + t99 ^ 2) + m(6) * (qJ(5) ^ 2 * t139 + pkin(4) ^ 2) + t133; 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * t40; m(6) * t64 + m(7) * t34 + t161; m(6) * t101 + m(7) * t88 + t127; -m(6) * pkin(4) + m(7) * t99 + t127; m(6) + m(7); mrSges(7,1) * t5 - mrSges(7,2) * t6; mrSges(7,1) * t2 - mrSges(7,2) * t3 + t136; mrSges(7,1) * t47 - mrSges(7,2) * t48 + t151; mrSges(7,1) * t61 - mrSges(7,2) * t62 + t151; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
