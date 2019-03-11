% Calculate joint inertia matrix for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:22:58
% EndTime: 2019-03-09 06:23:00
% DurationCPUTime: 1.10s
% Computational Cost: add. (1247->222), mult. (2227->297), div. (0->0), fcn. (2098->6), ass. (0->95)
t174 = Ifges(6,1) + Ifges(7,1);
t173 = Ifges(7,4) + Ifges(6,5);
t101 = cos(qJ(5));
t98 = sin(qJ(5));
t135 = t101 ^ 2 + t98 ^ 2;
t172 = mrSges(6,3) + mrSges(7,2);
t171 = Ifges(7,2) + Ifges(6,3);
t154 = Ifges(7,5) * t98;
t155 = Ifges(6,4) * t98;
t100 = sin(qJ(3));
t102 = cos(qJ(4));
t103 = cos(qJ(3));
t99 = sin(qJ(4));
t65 = t100 * t99 - t102 * t103;
t66 = t100 * t102 + t103 * t99;
t170 = t173 * t66 + (-t174 * t101 - t154 + t155) * t65;
t73 = -mrSges(6,1) * t101 + mrSges(6,2) * t98;
t169 = t73 - mrSges(5,1);
t132 = Ifges(7,5) * t101;
t133 = Ifges(6,4) * t101;
t168 = t174 * t98 - t132 + t133;
t167 = t173 * t98 + (Ifges(6,6) - Ifges(7,6)) * t101;
t83 = pkin(3) * t99 + pkin(9);
t166 = t135 * t83;
t165 = t135 * t66;
t164 = t103 ^ 2;
t163 = t172 * t135;
t104 = -pkin(1) - pkin(7);
t142 = -pkin(8) + t104;
t126 = t142 * t103;
t68 = t142 * t100;
t30 = -t102 * t126 + t68 * t99;
t162 = t30 ^ 2;
t61 = t65 ^ 2;
t161 = -2 * mrSges(5,3);
t72 = -mrSges(7,1) * t101 - mrSges(7,3) * t98;
t160 = 0.2e1 * t72;
t156 = m(7) * t98;
t153 = Ifges(6,6) * t66;
t152 = pkin(3) * t102;
t81 = t100 * pkin(3) + qJ(2);
t23 = pkin(4) * t66 + pkin(9) * t65 + t81;
t32 = t102 * t68 + t126 * t99;
t8 = t101 * t32 + t98 * t23;
t3 = qJ(6) * t66 + t8;
t151 = t101 * t3;
t150 = t101 * t8;
t149 = t30 * t65;
t148 = t65 * t98;
t86 = t98 * mrSges(7,2);
t143 = t98 * mrSges(6,3);
t141 = t166 * t66;
t140 = t165 * pkin(9);
t139 = t166 * pkin(9);
t138 = t135 * t83 ^ 2;
t136 = t135 * pkin(9) ^ 2;
t134 = t100 ^ 2 + t164;
t131 = t101 * t65;
t130 = qJ(6) * t101;
t129 = m(4) * t134;
t127 = t134 * mrSges(4,3);
t26 = -t66 * mrSges(7,1) - mrSges(7,2) * t131;
t123 = -Ifges(7,6) * t148 - t173 * t131 + t171 * t66;
t7 = t101 * t23 - t32 * t98;
t4 = -pkin(5) * t66 - t7;
t122 = t4 * t98 + t151;
t121 = -t7 * t98 + t150;
t120 = -t98 * mrSges(6,1) - t101 * mrSges(6,2);
t119 = -t98 * mrSges(7,1) + t101 * mrSges(7,3);
t118 = -pkin(5) * t98 + t130;
t117 = t102 * t65 - t66 * t99;
t69 = -pkin(5) * t101 - qJ(6) * t98 - pkin(4);
t74 = -Ifges(7,3) * t101 + t154;
t75 = Ifges(6,2) * t101 + t155;
t116 = Ifges(5,3) + t168 * t98 + (-t74 + t75) * t101;
t115 = (mrSges(5,1) * t102 - mrSges(5,2) * t99) * pkin(3);
t24 = -t66 * mrSges(6,2) + t143 * t65;
t25 = mrSges(6,1) * t66 + mrSges(6,3) * t131;
t27 = t66 * mrSges(7,3) + t65 * t86;
t113 = (-t25 + t26) * t98 + (t24 + t27) * t101;
t112 = -t66 * mrSges(5,2) + (t72 + t169) * t65 + t172 * t165;
t111 = mrSges(7,2) * t130 - pkin(5) * t86 + t167;
t110 = 0.2e1 * t163;
t109 = m(7) * t118 + t119 + t120;
t14 = Ifges(7,6) * t66 + (-Ifges(7,3) * t98 - t132) * t65;
t15 = t153 + (Ifges(6,2) * t98 - t133) * t65;
t9 = t118 * t65 + t30;
t108 = -t32 * mrSges(5,2) + mrSges(7,2) * t151 + mrSges(6,3) * t150 - Ifges(5,5) * t65 - t143 * t7 + t4 * t86 + t9 * t72 + t169 * t30 + t170 * t98 / 0.2e1 + (-t74 / 0.2e1 + t75 / 0.2e1) * t148 - t168 * t131 / 0.2e1 + (-t14 / 0.2e1 + t15 / 0.2e1) * t101 + (-Ifges(5,6) + t167 / 0.2e1) * t66;
t105 = qJ(2) ^ 2;
t84 = -pkin(4) - t152;
t60 = t66 ^ 2;
t53 = t69 - t152;
t22 = t120 * t65;
t21 = t119 * t65;
t1 = [Ifges(4,1) * t164 + 0.2e1 * t30 * t22 + 0.2e1 * t9 * t21 + 0.2e1 * t8 * t24 + 0.2e1 * t7 * t25 + 0.2e1 * t4 * t26 + 0.2e1 * t3 * t27 + Ifges(3,1) + Ifges(2,3) - (2 * pkin(1) * mrSges(3,2)) - 0.2e1 * t104 * t127 + (0.2e1 * mrSges(5,1) * t81 + Ifges(5,2) * t66 + t161 * t32 + t123) * t66 + (-0.2e1 * t81 * mrSges(5,2) + t30 * t161 + Ifges(5,1) * t65 + 0.2e1 * Ifges(5,4) * t66 - t170 * t101 + (-t14 + t15 + t153) * t98) * t65 + m(4) * (t104 ^ 2 * t134 + t105) + m(3) * ((pkin(1) ^ 2) + t105) + m(5) * (t32 ^ 2 + t81 ^ 2 + t162) + m(6) * (t7 ^ 2 + t8 ^ 2 + t162) + m(7) * (t3 ^ 2 + t4 ^ 2 + t9 ^ 2) + (-0.2e1 * Ifges(4,4) * t103 + Ifges(4,2) * t100) * t100 + 0.2e1 * (mrSges(4,1) * t100 + mrSges(4,2) * t103 + mrSges(3,3)) * qJ(2); -m(3) * pkin(1) - t61 * mrSges(5,3) + mrSges(3,2) + (t21 + t22) * t65 - t127 + (-mrSges(5,3) * t66 + t113) * t66 + m(7) * (t122 * t66 + t65 * t9) + m(6) * (t121 * t66 + t149) + m(5) * (t32 * t66 + t149) + t104 * t129; m(3) + m(5) * (t60 + t61) + t129 + 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * (t135 * t60 + t61); t113 * t83 + (mrSges(4,1) * t104 + Ifges(4,5)) * t103 + (-mrSges(4,2) * t104 - Ifges(4,6)) * t100 + t108 + t84 * t22 + t53 * t21 + m(7) * (t122 * t83 + t53 * t9) + m(6) * (t121 * t83 + t30 * t84) + (m(5) * (-t102 * t30 + t32 * t99) + t117 * mrSges(5,3)) * pkin(3); t103 * mrSges(4,1) - t100 * mrSges(4,2) + m(7) * (t53 * t65 + t141) + m(6) * (t65 * t84 + t141) - m(5) * t117 * pkin(3) + t112; t53 * t160 + 0.2e1 * t84 * t73 + Ifges(4,3) + 0.2e1 * t115 + m(7) * (t53 ^ 2 + t138) + m(6) * (t84 ^ 2 + t138) + m(5) * (t102 ^ 2 + t99 ^ 2) * pkin(3) ^ 2 + t110 * t83 + t116; t113 * pkin(9) + t108 + t69 * t21 - pkin(4) * t22 + m(6) * (-pkin(4) * t30 + pkin(9) * t121) + m(7) * (pkin(9) * t122 + t69 * t9); m(7) * (t65 * t69 + t140) + m(6) * (-pkin(4) * t65 + t140) + t112; (t84 - pkin(4)) * t73 + (t53 + t69) * t72 + t115 + m(7) * (t53 * t69 + t139) + m(6) * (-pkin(4) * t84 + t139) + t116 + (pkin(9) + t83) * t163; -0.2e1 * pkin(4) * t73 + t69 * t160 + m(7) * (t69 ^ 2 + t136) + m(6) * (pkin(4) ^ 2 + t136) + t110 * pkin(9) + t116; Ifges(6,6) * t148 - pkin(5) * t26 + m(7) * (-pkin(5) * t4 + qJ(6) * t3) + qJ(6) * t27 + t3 * mrSges(7,3) - t8 * mrSges(6,2) + t7 * mrSges(6,1) - t4 * mrSges(7,1) + t123; t109 * t66; t109 * t83 + t111; pkin(9) * t109 + t111; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t171; m(7) * t4 + t26; t66 * t156; t83 * t156 + t86; pkin(9) * t156 + t86; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
