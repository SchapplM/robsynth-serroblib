% Calculate joint inertia matrix for
% S6RPRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:10:52
% EndTime: 2019-03-09 06:10:55
% DurationCPUTime: 1.05s
% Computational Cost: add. (2027->218), mult. (3831->296), div. (0->0), fcn. (4328->8), ass. (0->97)
t162 = Ifges(6,1) + Ifges(7,1);
t161 = Ifges(7,4) + Ifges(6,5);
t160 = Ifges(7,2) + Ifges(6,3);
t96 = sin(qJ(5));
t142 = Ifges(7,5) * t96;
t144 = Ifges(6,4) * t96;
t100 = cos(qJ(4));
t101 = cos(qJ(3));
t94 = sin(pkin(10));
t95 = cos(pkin(10));
t98 = sin(qJ(3));
t59 = t101 * t95 - t94 * t98;
t60 = t101 * t94 + t95 * t98;
t97 = sin(qJ(4));
t49 = -t100 * t59 + t60 * t97;
t50 = t100 * t60 + t59 * t97;
t99 = cos(qJ(5));
t159 = (t162 * t99 + t142 - t144) * t50 + t161 * t49;
t141 = Ifges(7,5) * t99;
t143 = Ifges(6,4) * t99;
t158 = t162 * t96 - t141 + t143;
t127 = t96 ^ 2 + t99 ^ 2;
t157 = (Ifges(6,6) - Ifges(7,6)) * t99 + t161 * t96;
t156 = (mrSges(6,3) + mrSges(7,2)) * t127;
t135 = pkin(7) + qJ(2);
t65 = t135 * t94;
t66 = t135 * t95;
t51 = -t101 * t65 - t66 * t98;
t111 = -pkin(8) * t60 + t51;
t52 = t101 * t66 - t98 * t65;
t39 = pkin(8) * t59 + t52;
t22 = -t100 * t111 + t39 * t97;
t155 = t22 ^ 2;
t91 = t95 ^ 2;
t154 = 0.2e1 * t22;
t79 = -pkin(2) * t95 - pkin(1);
t53 = -pkin(3) * t59 + t79;
t153 = 0.2e1 * t53;
t152 = 0.2e1 * t59;
t69 = -t99 * mrSges(7,1) - t96 * mrSges(7,3);
t151 = 0.2e1 * t69;
t148 = m(7) * t96;
t24 = t100 * t39 + t111 * t97;
t25 = pkin(4) * t49 - pkin(9) * t50 + t53;
t7 = t99 * t24 + t96 * t25;
t3 = qJ(6) * t49 + t7;
t146 = t3 * t99;
t145 = t7 * t99;
t140 = Ifges(6,6) * t49;
t139 = pkin(3) * t100;
t138 = t50 * t96;
t137 = t50 * t99;
t84 = t96 * mrSges(7,2);
t136 = t96 * mrSges(6,3);
t28 = -mrSges(6,2) * t49 - t136 * t50;
t31 = mrSges(7,3) * t49 - t50 * t84;
t134 = t28 + t31;
t29 = mrSges(6,1) * t49 - mrSges(6,3) * t137;
t30 = -t49 * mrSges(7,1) + mrSges(7,2) * t137;
t133 = t29 - t30;
t80 = pkin(3) * t97 + pkin(9);
t132 = t127 * pkin(9) * t80;
t131 = t127 * t80 ^ 2;
t129 = t127 * pkin(9) ^ 2;
t128 = t94 ^ 2 + t91;
t126 = qJ(6) * t99;
t124 = -t95 * mrSges(3,1) + t94 * mrSges(3,2);
t123 = -t59 * mrSges(4,1) + t60 * mrSges(4,2);
t120 = Ifges(7,6) * t138 + t137 * t161 + t160 * t49;
t6 = -t24 * t96 + t25 * t99;
t4 = -pkin(5) * t49 - t6;
t119 = t4 * t96 + t146;
t118 = -t6 * t96 + t145;
t70 = -t99 * mrSges(6,1) + t96 * mrSges(6,2);
t117 = t96 * mrSges(6,1) + t99 * mrSges(6,2);
t116 = t96 * mrSges(7,1) - t99 * mrSges(7,3);
t115 = pkin(5) * t99 + qJ(6) * t96;
t114 = pkin(5) * t96 - t126;
t64 = -pkin(4) - t115;
t71 = -Ifges(7,3) * t99 + t142;
t72 = Ifges(6,2) * t99 + t144;
t113 = Ifges(5,3) + (-t71 + t72) * t99 + t158 * t96;
t112 = (mrSges(5,1) * t100 - mrSges(5,2) * t97) * pkin(3);
t109 = -t133 * t96 + t134 * t99;
t108 = mrSges(7,2) * t126 - pkin(5) * t84 + t157;
t107 = 0.2e1 * t156;
t106 = -m(7) * t114 - t116 - t117;
t14 = Ifges(7,6) * t49 + (Ifges(7,3) * t96 + t141) * t50;
t15 = t140 + (-Ifges(6,2) * t96 + t143) * t50;
t9 = t114 * t50 + t22;
t105 = -t24 * mrSges(5,2) + mrSges(7,2) * t146 + mrSges(6,3) * t145 + Ifges(5,5) * t50 - t136 * t6 + t4 * t84 + t9 * t69 + (-t14 / 0.2e1 + t15 / 0.2e1) * t99 + (t70 - mrSges(5,1)) * t22 + t159 * t96 / 0.2e1 + (t71 / 0.2e1 - t72 / 0.2e1) * t138 + t158 * t137 / 0.2e1 + (-Ifges(5,6) + t157 / 0.2e1) * t49;
t81 = -pkin(4) - t139;
t56 = t64 - t139;
t44 = t50 * mrSges(5,2);
t27 = t117 * t50;
t26 = t116 * t50;
t1 = [t52 * mrSges(4,3) * t152 - 0.2e1 * pkin(1) * t124 + Ifges(3,2) * t91 + 0.2e1 * t79 * t123 + t44 * t153 + Ifges(4,2) * t59 ^ 2 + t27 * t154 + 0.2e1 * t7 * t28 + 0.2e1 * t6 * t29 + 0.2e1 * t4 * t30 + 0.2e1 * t3 * t31 + 0.2e1 * t9 * t26 + Ifges(2,3) + (Ifges(3,1) * t94 + 0.2e1 * Ifges(3,4) * t95) * t94 + 0.2e1 * t128 * qJ(2) * mrSges(3,3) + (-0.2e1 * t51 * mrSges(4,3) + Ifges(4,1) * t60 + Ifges(4,4) * t152) * t60 + (mrSges(5,1) * t153 - 0.2e1 * t24 * mrSges(5,3) + Ifges(5,2) * t49 + t120) * t49 + (mrSges(5,3) * t154 + Ifges(5,1) * t50 - 0.2e1 * Ifges(5,4) * t49 + t159 * t99 + (t14 - t15 - t140) * t96) * t50 + m(3) * (qJ(2) ^ 2 * t128 + pkin(1) ^ 2) + m(4) * (t51 ^ 2 + t52 ^ 2 + t79 ^ 2) + m(5) * (t24 ^ 2 + t53 ^ 2 + t155) + m(7) * (t3 ^ 2 + t4 ^ 2 + t9 ^ 2) + m(6) * (t6 ^ 2 + t7 ^ 2 + t155); -m(3) * pkin(1) + t49 * mrSges(5,1) + t44 + t133 * t99 + t134 * t96 + m(7) * (t3 * t96 - t4 * t99) + m(6) * (t6 * t99 + t7 * t96) + m(5) * t53 + m(4) * t79 + t123 + t124; m(3) + m(4) + m(5) + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * t127; t105 + (m(5) * (-t100 * t22 + t24 * t97) + (-t100 * t50 - t97 * t49) * mrSges(5,3)) * pkin(3) + m(6) * (t118 * t80 + t22 * t81) + m(7) * (t119 * t80 + t56 * t9) + t109 * t80 + t81 * t27 + t56 * t26 + Ifges(4,6) * t59 + Ifges(4,5) * t60 + t51 * mrSges(4,1) - t52 * mrSges(4,2); 0; t56 * t151 + 0.2e1 * t81 * t70 + Ifges(4,3) + 0.2e1 * t112 + m(7) * (t56 ^ 2 + t131) + m(6) * (t81 ^ 2 + t131) + m(5) * (t100 ^ 2 + t97 ^ 2) * pkin(3) ^ 2 + t107 * t80 + t113; t109 * pkin(9) + t105 + m(7) * (pkin(9) * t119 + t64 * t9) + m(6) * (-pkin(4) * t22 + pkin(9) * t118) + t64 * t26 - pkin(4) * t27; 0; (t81 - pkin(4)) * t70 + (t56 + t64) * t69 + t112 + m(7) * (t56 * t64 + t132) + m(6) * (-pkin(4) * t81 + t132) + t113 + (pkin(9) + t80) * t156; -0.2e1 * pkin(4) * t70 + t64 * t151 + m(7) * (t64 ^ 2 + t129) + m(6) * (pkin(4) ^ 2 + t129) + t107 * pkin(9) + t113; -Ifges(6,6) * t138 - pkin(5) * t30 + m(7) * (-pkin(5) * t4 + qJ(6) * t3) + qJ(6) * t31 + t3 * mrSges(7,3) + t6 * mrSges(6,1) - t4 * mrSges(7,1) - t7 * mrSges(6,2) + t120; m(7) * t115 - t69 - t70; t106 * t80 + t108; pkin(9) * t106 + t108; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t160; m(7) * t4 + t30; -m(7) * t99; t148 * t80 + t84; pkin(9) * t148 + t84; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
