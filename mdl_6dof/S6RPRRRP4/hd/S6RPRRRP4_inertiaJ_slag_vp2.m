% Calculate joint inertia matrix for
% S6RPRRRP4
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
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:06:59
% EndTime: 2019-03-09 06:07:02
% DurationCPUTime: 1.28s
% Computational Cost: add. (2011->232), mult. (3796->316), div. (0->0), fcn. (4301->8), ass. (0->93)
t161 = Ifges(6,4) + Ifges(7,4);
t156 = Ifges(7,6) + Ifges(6,6);
t100 = cos(qJ(4));
t101 = cos(qJ(3));
t94 = sin(pkin(10));
t95 = cos(pkin(10));
t98 = sin(qJ(3));
t60 = t101 * t95 - t94 * t98;
t61 = t101 * t94 + t95 * t98;
t97 = sin(qJ(4));
t49 = -t100 * t60 + t61 * t97;
t160 = t156 * t49;
t159 = Ifges(6,1) + Ifges(7,1);
t158 = Ifges(6,5) + Ifges(7,5);
t157 = Ifges(6,2) + Ifges(7,2);
t99 = cos(qJ(5));
t155 = t161 * t99;
t96 = sin(qJ(5));
t154 = t161 * t96;
t153 = Ifges(6,3) + Ifges(7,3);
t50 = t100 * t61 + t60 * t97;
t152 = (-t157 * t96 + t155) * t50 + t160;
t151 = (t159 * t99 - t154) * t50 + t158 * t49;
t150 = t157 * t99 + t154;
t149 = t159 * t96 + t155;
t113 = t156 * t99 + t158 * t96;
t121 = t96 ^ 2 + t99 ^ 2;
t148 = 0.2e1 * t121;
t125 = pkin(7) + qJ(2);
t66 = t125 * t94;
t67 = t125 * t95;
t51 = -t101 * t66 - t67 * t98;
t106 = -pkin(8) * t61 + t51;
t52 = t101 * t67 - t98 * t66;
t39 = pkin(8) * t60 + t52;
t22 = -t100 * t106 + t39 * t97;
t147 = t22 ^ 2;
t91 = t95 ^ 2;
t146 = 0.2e1 * t22;
t79 = -pkin(2) * t95 - pkin(1);
t53 = -pkin(3) * t60 + t79;
t145 = 0.2e1 * t53;
t144 = 0.2e1 * t60;
t85 = t96 * mrSges(7,2);
t70 = -t99 * mrSges(7,1) + t85;
t143 = 0.2e1 * t70;
t142 = m(7) * pkin(5);
t138 = pkin(5) * t99;
t24 = t100 * t39 + t106 * t97;
t25 = pkin(4) * t49 - pkin(9) * t50 + t53;
t6 = t99 * t24 + t96 * t25;
t137 = t6 * t99;
t136 = mrSges(6,2) * t96;
t131 = t50 * t96;
t130 = t50 * t99;
t129 = t96 * mrSges(7,3);
t128 = t99 * mrSges(7,3);
t29 = -mrSges(6,2) * t49 - mrSges(6,3) * t131;
t127 = t99 * t29;
t26 = mrSges(7,1) * t131 + mrSges(7,2) * t130;
t122 = t94 ^ 2 + t91;
t84 = t99 * qJ(6);
t120 = 0.2e1 * mrSges(7,3);
t81 = -pkin(3) * t100 - pkin(4);
t80 = pkin(3) * t97 + pkin(9);
t117 = t121 * t80;
t116 = -t95 * mrSges(3,1) + t94 * mrSges(3,2);
t115 = -t60 * mrSges(4,1) + t61 * mrSges(4,2);
t5 = -t24 * t96 + t99 * t25;
t114 = t158 * t130 + t153 * t49;
t112 = t149 * t96 + t150 * t99 + Ifges(5,3);
t1 = pkin(5) * t49 - t50 * t84 + t5;
t111 = -t5 * mrSges(6,3) - t1 * mrSges(7,3);
t110 = -t5 * t96 + t137;
t109 = mrSges(6,1) * t96 + mrSges(6,2) * t99;
t108 = mrSges(6,3) * t148;
t107 = (mrSges(5,1) * t100 - mrSges(5,2) * t97) * pkin(3);
t3 = -qJ(6) * t131 + t6;
t71 = -mrSges(6,1) * t99 + t136;
t8 = pkin(5) * t131 + t22;
t105 = -t24 * mrSges(5,2) + mrSges(6,3) * t137 + Ifges(5,5) * t50 + t3 * t128 + t8 * t70 + (-mrSges(5,1) + t71) * t22 + t151 * t96 / 0.2e1 + t152 * t99 / 0.2e1 - t150 * t131 / 0.2e1 + t149 * t130 / 0.2e1 + (-Ifges(5,6) + t113 / 0.2e1) * t49;
t82 = -pkin(4) - t138;
t72 = pkin(9) * t99 + t84;
t69 = (-qJ(6) - pkin(9)) * t96;
t68 = t81 - t138;
t57 = t80 * t99 + t84;
t56 = (-qJ(6) - t80) * t96;
t44 = t50 * mrSges(5,2);
t31 = mrSges(6,1) * t49 - mrSges(6,3) * t130;
t30 = mrSges(7,1) * t49 - t128 * t50;
t28 = -mrSges(7,2) * t49 - t129 * t50;
t27 = t109 * t50;
t2 = [t52 * mrSges(4,3) * t144 + Ifges(3,2) * t91 - 0.2e1 * pkin(1) * t116 + 0.2e1 * t79 * t115 + Ifges(4,2) * t60 ^ 2 + t44 * t145 + 0.2e1 * t8 * t26 + t27 * t146 + 0.2e1 * t3 * t28 + 0.2e1 * t6 * t29 + 0.2e1 * t1 * t30 + 0.2e1 * t5 * t31 + Ifges(2,3) + (Ifges(3,1) * t94 + 0.2e1 * Ifges(3,4) * t95) * t94 + 0.2e1 * t122 * qJ(2) * mrSges(3,3) + (-0.2e1 * t51 * mrSges(4,3) + Ifges(4,1) * t61 + Ifges(4,4) * t144) * t61 + (mrSges(5,1) * t145 - 0.2e1 * t24 * mrSges(5,3) + Ifges(5,2) * t49 + t114) * t49 + (mrSges(5,3) * t146 + Ifges(5,1) * t50 - 0.2e1 * Ifges(5,4) * t49 + t151 * t99 + (-t152 - t160) * t96) * t50 + m(3) * (qJ(2) ^ 2 * t122 + pkin(1) ^ 2) + m(4) * (t51 ^ 2 + t52 ^ 2 + t79 ^ 2) + m(5) * (t24 ^ 2 + t53 ^ 2 + t147) + m(7) * (t1 ^ 2 + t3 ^ 2 + t8 ^ 2) + m(6) * (t5 ^ 2 + t6 ^ 2 + t147); -m(3) * pkin(1) + t49 * mrSges(5,1) + t44 + (t30 + t31) * t99 + (t28 + t29) * t96 + m(6) * (t5 * t99 + t6 * t96) + m(7) * (t1 * t99 + t3 * t96) + m(5) * t53 + m(4) * t79 + t115 + t116; m(3) + m(4) + m(5) + (m(6) / 0.2e1 + m(7) / 0.2e1) * t148; (m(5) * (-t100 * t22 + t24 * t97) + (-t100 * t50 - t97 * t49) * mrSges(5,3)) * pkin(3) + t105 + m(6) * (t110 * t80 + t22 * t81) + (-t80 * t31 + t111) * t96 + m(7) * (t1 * t56 + t3 * t57 + t68 * t8) + t80 * t127 + t81 * t27 + Ifges(4,5) * t61 + t68 * t26 - t52 * mrSges(4,2) + t56 * t30 + t57 * t28 + Ifges(4,6) * t60 + t51 * mrSges(4,1); m(7) * (t56 * t99 + t57 * t96); t68 * t143 + 0.2e1 * t81 * t71 + Ifges(4,3) + 0.2e1 * t107 + (-t56 * t96 + t57 * t99) * t120 + t80 * t108 + m(7) * (t56 ^ 2 + t57 ^ 2 + t68 ^ 2) + m(6) * (t121 * t80 ^ 2 + t81 ^ 2) + m(5) * (t100 ^ 2 + t97 ^ 2) * pkin(3) ^ 2 + t112; (-pkin(9) * t31 + t111) * t96 + m(7) * (t1 * t69 + t3 * t72 + t8 * t82) + t105 + m(6) * (-pkin(4) * t22 + pkin(9) * t110) + pkin(9) * t127 + t82 * t26 + t69 * t30 + t72 * t28 - pkin(4) * t27; m(7) * (t69 * t99 + t72 * t96); (t81 - pkin(4)) * t71 + (t68 + t82) * t70 + t107 + m(7) * (t56 * t69 + t57 * t72 + t68 * t82) + m(6) * (-pkin(4) * t81 + pkin(9) * t117) + ((t57 + t72) * t99 + (-t56 - t69) * t96) * mrSges(7,3) + (pkin(9) * t121 + t117) * mrSges(6,3) + t112; -0.2e1 * pkin(4) * t71 + t82 * t143 + (-t69 * t96 + t72 * t99) * t120 + pkin(9) * t108 + m(7) * (t69 ^ 2 + t72 ^ 2 + t82 ^ 2) + m(6) * (pkin(9) ^ 2 * t121 + pkin(4) ^ 2) + t112; mrSges(6,1) * t5 + mrSges(7,1) * t1 - mrSges(6,2) * t6 - mrSges(7,2) * t3 - t156 * t131 + (m(7) * t1 + t30) * pkin(5) + t114; -t136 - t85 + (mrSges(6,1) + mrSges(7,1) + t142) * t99; mrSges(7,1) * t56 - mrSges(7,2) * t57 - t109 * t80 + (m(7) * t56 - t129) * pkin(5) + t113; mrSges(7,1) * t69 - mrSges(7,2) * t72 - t109 * pkin(9) + (m(7) * t69 - t129) * pkin(5) + t113; (0.2e1 * mrSges(7,1) + t142) * pkin(5) + t153; m(7) * t8 + t26; 0; m(7) * t68 + t70; m(7) * t82 + t70; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
