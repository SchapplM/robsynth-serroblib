% Calculate joint inertia matrix for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:43:37
% EndTime: 2019-03-09 16:43:40
% DurationCPUTime: 1.13s
% Computational Cost: add. (1348->235), mult. (2402->301), div. (0->0), fcn. (2285->6), ass. (0->91)
t161 = Ifges(6,1) + Ifges(7,1);
t160 = Ifges(7,4) + Ifges(6,5);
t159 = Ifges(7,2) + Ifges(6,3);
t96 = cos(qJ(5));
t140 = Ifges(7,5) * t96;
t142 = Ifges(6,4) * t96;
t94 = sin(qJ(3));
t95 = sin(qJ(2));
t97 = cos(qJ(3));
t98 = cos(qJ(2));
t56 = t94 * t95 - t97 * t98;
t57 = t94 * t98 + t95 * t97;
t93 = sin(qJ(5));
t158 = t160 * t57 + (t161 * t93 - t140 + t142) * t56;
t141 = Ifges(7,5) * t93;
t143 = Ifges(6,4) * t93;
t157 = t161 * t96 + t141 - t143;
t127 = -t93 ^ 2 - t96 ^ 2;
t156 = t160 * t96 + (-Ifges(6,6) + Ifges(7,6)) * t93;
t66 = mrSges(6,1) * t93 + mrSges(6,2) * t96;
t155 = mrSges(5,3) + t66;
t154 = (mrSges(6,3) + mrSges(7,2)) * t127;
t144 = pkin(2) * t94;
t76 = qJ(4) + t144;
t153 = t76 ^ 2;
t80 = -pkin(2) * t98 - pkin(1);
t113 = -qJ(4) * t57 + t80;
t24 = pkin(3) * t56 + t113;
t152 = -0.2e1 * t24;
t65 = mrSges(7,1) * t93 - mrSges(7,3) * t96;
t151 = 0.2e1 * t65;
t150 = 0.2e1 * t80;
t147 = pkin(3) + pkin(9);
t146 = -pkin(8) - pkin(7);
t145 = m(7) * t96;
t139 = t56 * t93;
t138 = t56 * t96;
t136 = t93 * mrSges(7,2);
t83 = t96 * mrSges(7,2);
t18 = t147 * t56 + t113;
t71 = t146 * t95;
t72 = t146 * t98;
t35 = -t97 * t71 - t72 * t94;
t20 = pkin(4) * t57 + t35;
t5 = t96 * t18 + t93 * t20;
t25 = mrSges(6,1) * t57 - mrSges(6,3) * t139;
t26 = -t57 * mrSges(7,1) + t56 * t136;
t134 = t25 - t26;
t27 = -mrSges(6,2) * t57 + mrSges(6,3) * t138;
t28 = mrSges(7,3) * t57 + t56 * t83;
t133 = t27 + t28;
t79 = -pkin(2) * t97 - pkin(3);
t75 = -pkin(9) + t79;
t132 = t127 * t75 * t147;
t131 = t127 * t75 ^ 2;
t130 = t127 * t147 ^ 2;
t128 = t95 ^ 2 + t98 ^ 2;
t126 = qJ(4) * t76;
t37 = t71 * t94 - t72 * t97;
t123 = t35 ^ 2 + t37 ^ 2;
t121 = 0.2e1 * t155;
t119 = Ifges(6,6) * t138 + t160 * t139 + t159 * t57;
t2 = qJ(6) * t57 + t5;
t4 = -t18 * t93 + t20 * t96;
t3 = -pkin(5) * t57 - t4;
t118 = t2 * t93 - t3 * t96;
t117 = t4 * t96 + t5 * t93;
t62 = t93 * pkin(5) - qJ(6) * t96 + qJ(4);
t116 = -t96 * mrSges(6,1) + t93 * mrSges(6,2);
t115 = -t96 * mrSges(7,1) - t93 * mrSges(7,3);
t114 = pkin(5) * t96 + qJ(6) * t93;
t112 = (mrSges(4,1) * t97 - mrSges(4,2) * t94) * pkin(2);
t111 = -0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * t127;
t67 = Ifges(7,3) * t93 + t140;
t68 = -Ifges(6,2) * t93 + t142;
t110 = Ifges(5,1) + Ifges(4,3) + t157 * t96 + (t67 - t68) * t93;
t109 = t133 * t93 + t134 * t96;
t107 = mrSges(5,2) + t154;
t106 = 0.2e1 * t154;
t105 = -t114 * mrSges(7,2) + t156;
t104 = m(7) * t114 - t115 - t116;
t14 = t57 * Ifges(7,6) + (-Ifges(7,3) * t96 + t141) * t56;
t15 = Ifges(6,6) * t57 + (Ifges(6,2) * t96 + t143) * t56;
t21 = -pkin(4) * t56 + t37;
t7 = (-pkin(4) - t114) * t56 + t37;
t103 = t3 * t83 + t21 * t66 + t7 * t65 + (-mrSges(4,2) + mrSges(5,3)) * t37 + (t14 / 0.2e1 - t15 / 0.2e1) * t93 + (Ifges(5,5) - Ifges(4,6)) * t56 + (mrSges(5,2) - mrSges(4,1)) * t35 + t158 * t96 / 0.2e1 + (-t67 / 0.2e1 + t68 / 0.2e1) * t138 + t157 * t139 / 0.2e1 + (Ifges(4,5) - Ifges(5,4) + t156 / 0.2e1) * t57;
t100 = qJ(4) ^ 2;
t46 = t62 + t144;
t23 = t116 * t56;
t22 = t115 * t56;
t1 = [-0.2e1 * pkin(1) * (-t98 * mrSges(3,1) + t95 * mrSges(3,2)) + t95 * (Ifges(3,1) * t95 + Ifges(3,4) * t98) + t98 * (Ifges(3,4) * t95 + Ifges(3,2) * t98) + 0.2e1 * t21 * t23 + 0.2e1 * t4 * t25 + 0.2e1 * t3 * t26 + 0.2e1 * t5 * t27 + 0.2e1 * t2 * t28 + 0.2e1 * t7 * t22 + Ifges(2,3) + m(3) * (t128 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t80 ^ 2 + t123) + m(5) * (t24 ^ 2 + t123) + m(6) * (t21 ^ 2 + t4 ^ 2 + t5 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t7 ^ 2) + (mrSges(4,2) * t150 + mrSges(5,3) * t152 + (Ifges(4,1) + Ifges(5,2)) * t57 + t119) * t57 + (mrSges(4,1) * t150 + mrSges(5,2) * t152 + (Ifges(4,2) + Ifges(5,3)) * t56 + (-t14 + t15) * t96 + t158 * t93 + (-Ifges(7,6) * t96 - (2 * Ifges(4,4)) - (2 * Ifges(5,6))) * t57) * t56 + 0.2e1 * t128 * pkin(7) * mrSges(3,3) + 0.2e1 * (t35 * t57 - t37 * t56) * (mrSges(5,1) + mrSges(4,3)); t103 + m(6) * (t117 * t75 + t21 * t76) + m(7) * (t118 * t75 + t46 * t7) + t109 * t75 + m(5) * (t35 * t79 + t37 * t76) + (-t95 * mrSges(3,1) - t98 * mrSges(3,2)) * pkin(7) + Ifges(3,5) * t95 + Ifges(3,6) * t98 + t76 * t23 - t2 * t136 + t46 * t22 + (-t76 * t56 + t79 * t57) * mrSges(5,1) - t117 * mrSges(6,3) + ((-t56 * t94 - t57 * t97) * mrSges(4,3) + m(4) * (-t35 * t97 + t37 * t94)) * pkin(2); 0.2e1 * t79 * mrSges(5,2) + t46 * t151 + Ifges(3,3) + t76 * t121 + 0.2e1 * t112 + t106 * t75 + m(7) * (t46 ^ 2 - t131) + m(6) * (-t131 + t153) + m(5) * (t79 ^ 2 + t153) + m(4) * (t94 ^ 2 + t97 ^ 2) * pkin(2) ^ 2 + t110; t103 + m(7) * (-t118 * t147 + t62 * t7) + m(6) * (qJ(4) * t21 - t117 * t147) + t62 * t22 + qJ(4) * t23 + (-t4 * mrSges(6,3) - t134 * t147) * t96 + (-t2 * mrSges(7,2) - t5 * mrSges(6,3) - t133 * t147) * t93 + (-pkin(3) * t57 - qJ(4) * t56) * mrSges(5,1) + m(5) * (-pkin(3) * t35 + qJ(4) * t37); (t62 + t46) * t65 + t112 + (t79 - pkin(3)) * mrSges(5,2) + m(5) * (-pkin(3) * t79 + t126) + m(6) * (t126 + t132) + m(7) * (t46 * t62 + t132) + t110 + t155 * (t76 + qJ(4)) + (t75 - t147) * t154; -0.2e1 * pkin(3) * mrSges(5,2) + t62 * t151 + qJ(4) * t121 + m(7) * (t62 ^ 2 - t130) + m(6) * (t100 - t130) + m(5) * (pkin(3) ^ 2 + t100) - t106 * t147 + t110; m(5) * t35 + m(6) * t117 + m(7) * t118 + t57 * mrSges(5,1) + t109; m(5) * t79 + t75 * t111 + t107; -m(5) * pkin(3) - t111 * t147 + t107; m(5) + t111; -Ifges(7,6) * t138 + m(7) * (-pkin(5) * t3 + qJ(6) * t2) - t5 * mrSges(6,2) - t3 * mrSges(7,1) - pkin(5) * t26 + t2 * mrSges(7,3) + qJ(6) * t28 + t4 * mrSges(6,1) + t119; t104 * t75 + t105; -t104 * t147 + t105; t104; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t159; m(7) * t3 + t26; -t75 * t145 + t83; t145 * t147 + t83; -t145; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
