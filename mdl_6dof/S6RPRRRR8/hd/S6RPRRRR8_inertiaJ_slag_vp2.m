% Calculate joint inertia matrix for
% S6RPRRRR8
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:19:47
% EndTime: 2019-03-09 07:19:50
% DurationCPUTime: 1.19s
% Computational Cost: add. (1905->274), mult. (3418->386), div. (0->0), fcn. (3591->8), ass. (0->109)
t109 = sin(qJ(5));
t113 = cos(qJ(5));
t145 = t109 ^ 2 + t113 ^ 2;
t172 = mrSges(6,3) * t145;
t88 = -mrSges(6,1) * t113 + mrSges(6,2) * t109;
t171 = t88 - mrSges(5,1);
t115 = cos(qJ(3));
t170 = t115 ^ 2;
t110 = sin(qJ(4));
t111 = sin(qJ(3));
t114 = cos(qJ(4));
t82 = t110 * t111 - t114 * t115;
t148 = t113 * t82;
t83 = t110 * t115 + t111 * t114;
t169 = -Ifges(6,5) * t148 + Ifges(6,3) * t83;
t108 = sin(qJ(6));
t146 = Ifges(6,5) * t109 + Ifges(6,6) * t113;
t112 = cos(qJ(6));
t126 = t108 * t109 - t112 * t113;
t157 = t126 * mrSges(7,3);
t168 = -t108 * pkin(5) * t157 + t146;
t116 = -pkin(1) - pkin(7);
t155 = -pkin(8) + t116;
t138 = t155 * t115;
t86 = t155 * t111;
t52 = t110 * t86 - t114 * t138;
t167 = t52 ^ 2;
t166 = t82 ^ 2;
t165 = -2 * mrSges(5,3);
t84 = t108 * t113 + t109 * t112;
t49 = mrSges(7,1) * t126 + mrSges(7,2) * t84;
t164 = 0.2e1 * t49;
t160 = Ifges(6,6) * t83;
t159 = pkin(3) * t114;
t158 = t52 * t82;
t156 = t84 * mrSges(7,3);
t94 = t111 * pkin(3) + qJ(2);
t44 = pkin(4) * t83 + pkin(9) * t82 + t94;
t54 = t110 * t138 + t114 * t86;
t15 = t109 * t44 + t113 * t54;
t154 = Ifges(7,5) * t84 - Ifges(7,6) * t126;
t153 = Ifges(6,4) * t109;
t152 = Ifges(6,4) * t113;
t14 = -t109 * t54 + t113 * t44;
t151 = t109 * t14;
t150 = t109 * t82;
t149 = t113 * t15;
t96 = pkin(3) * t110 + pkin(9);
t147 = t113 * t96;
t144 = t111 ^ 2 + t170;
t143 = 0.2e1 * mrSges(7,3);
t35 = t84 * t82;
t37 = t126 * t82;
t142 = Ifges(7,5) * t37 + Ifges(7,6) * t35 + Ifges(7,3) * t83;
t141 = t112 * t156;
t98 = -pkin(5) * t113 - pkin(4);
t34 = t84 * t83;
t36 = t126 * t83;
t140 = -t34 * mrSges(7,1) + t36 * mrSges(7,2);
t139 = m(4) * t144;
t137 = t145 * pkin(9);
t136 = t145 * t96;
t135 = t144 * mrSges(4,3);
t50 = Ifges(7,4) * t84 - Ifges(7,2) * t126;
t51 = Ifges(7,1) * t84 - Ifges(7,4) * t126;
t89 = Ifges(6,2) * t113 + t153;
t90 = Ifges(6,1) * t109 + t152;
t134 = t109 * t90 + t113 * t89 - t126 * t50 + t84 * t51 + Ifges(5,3);
t133 = -mrSges(6,1) * t109 - mrSges(6,2) * t113;
t132 = t149 - t151;
t45 = -mrSges(6,2) * t83 + mrSges(6,3) * t150;
t46 = mrSges(6,1) * t83 + mrSges(6,3) * t148;
t131 = -t109 * t46 + t113 * t45;
t130 = t110 * t83 - t114 * t82;
t73 = (-pkin(10) - t96) * t109;
t102 = t113 * pkin(10);
t74 = t102 + t147;
t41 = -t108 * t74 + t112 * t73;
t42 = t108 * t73 + t112 * t74;
t129 = t41 * mrSges(7,1) - t42 * mrSges(7,2) + t154;
t91 = (-pkin(10) - pkin(9)) * t109;
t92 = pkin(9) * t113 + t102;
t58 = -t108 * t92 + t112 * t91;
t59 = t108 * t91 + t112 * t92;
t128 = t58 * mrSges(7,1) - t59 * mrSges(7,2) + t154;
t127 = 0.2e1 * t172;
t10 = pkin(10) * t150 + t15;
t7 = pkin(5) * t83 + pkin(10) * t148 + t14;
t3 = -t10 * t108 + t112 * t7;
t4 = t10 * t112 + t108 * t7;
t125 = t3 * mrSges(7,1) - t4 * mrSges(7,2) + t142;
t124 = (mrSges(5,1) * t114 - mrSges(5,2) * t110) * pkin(3);
t123 = (mrSges(7,1) * t112 - mrSges(7,2) * t108) * pkin(5);
t122 = t34 * t156 + t36 * t157 + (t49 + t171) * t82 + (-mrSges(5,2) + t172) * t83;
t23 = t160 + (Ifges(6,2) * t109 - t152) * t82;
t24 = Ifges(6,5) * t83 + (-Ifges(6,1) * t113 + t153) * t82;
t25 = -pkin(5) * t150 + t52;
t8 = Ifges(7,4) * t37 + Ifges(7,2) * t35 + Ifges(7,6) * t83;
t9 = Ifges(7,1) * t37 + Ifges(7,4) * t35 + Ifges(7,5) * t83;
t121 = -t54 * mrSges(5,2) - t3 * t156 - t4 * t157 + t25 * t49 + mrSges(6,3) * t149 + t35 * t50 / 0.2e1 + t37 * t51 / 0.2e1 + t109 * t24 / 0.2e1 + t113 * t23 / 0.2e1 - t126 * t8 / 0.2e1 + t89 * t150 / 0.2e1 - t90 * t148 / 0.2e1 + t84 * t9 / 0.2e1 - Ifges(5,6) * t83 - Ifges(5,5) * t82 + t171 * t52 + (t154 + t146) * t83 / 0.2e1;
t117 = qJ(2) ^ 2;
t97 = -pkin(4) - t159;
t87 = t98 - t159;
t77 = t83 ^ 2;
t43 = t133 * t82;
t19 = mrSges(7,1) * t83 - mrSges(7,3) * t37;
t18 = -mrSges(7,2) * t83 + mrSges(7,3) * t35;
t11 = -mrSges(7,1) * t35 + mrSges(7,2) * t37;
t1 = [Ifges(4,1) * t170 + 0.2e1 * t15 * t45 + 0.2e1 * t14 * t46 + 0.2e1 * t52 * t43 + t35 * t8 + t37 * t9 + 0.2e1 * t4 * t18 + 0.2e1 * t3 * t19 + 0.2e1 * t25 * t11 + Ifges(3,1) + Ifges(2,3) - (2 * pkin(1) * mrSges(3,2)) - 0.2e1 * t116 * t135 + (0.2e1 * t94 * mrSges(5,1) + Ifges(5,2) * t83 + t54 * t165 + t142 + t169) * t83 + (-0.2e1 * t94 * mrSges(5,2) + t52 * t165 + Ifges(5,1) * t82 + 0.2e1 * Ifges(5,4) * t83 - t113 * t24 + (t23 + t160) * t109) * t82 + m(4) * (t144 * t116 ^ 2 + t117) + m(3) * ((pkin(1) ^ 2) + t117) + m(5) * (t54 ^ 2 + t94 ^ 2 + t167) + m(6) * (t14 ^ 2 + t15 ^ 2 + t167) + m(7) * (t25 ^ 2 + t3 ^ 2 + t4 ^ 2) + (-0.2e1 * Ifges(4,4) * t115 + Ifges(4,2) * t111) * t111 + 0.2e1 * (mrSges(4,1) * t111 + mrSges(4,2) * t115 + mrSges(3,3)) * qJ(2); -m(3) * pkin(1) - t166 * mrSges(5,3) - t36 * t18 - t34 * t19 + mrSges(3,2) + (t11 + t43) * t82 - t135 + (-mrSges(5,3) * t83 + t131) * t83 + m(7) * (t25 * t82 - t3 * t34 - t36 * t4) + m(6) * (t132 * t83 + t158) + m(5) * (t54 * t83 + t158) + t116 * t139; m(3) + m(7) * (t34 ^ 2 + t36 ^ 2 + t166) + m(6) * (t145 * t77 + t166) + m(5) * (t77 + t166) + t139; t121 + (m(5) * (t110 * t54 - t114 * t52) - t130 * mrSges(5,3)) * pkin(3) + m(6) * (t132 * t96 + t52 * t97) + (-t14 * mrSges(6,3) - t96 * t46) * t109 + m(7) * (t25 * t87 + t3 * t41 + t4 * t42) + (t116 * mrSges(4,1) + Ifges(4,5)) * t115 + t45 * t147 + t97 * t43 + t87 * t11 + t42 * t18 + (-t116 * mrSges(4,2) - Ifges(4,6)) * t111 + t41 * t19; t115 * mrSges(4,1) - t111 * mrSges(4,2) + m(7) * (-t34 * t41 - t36 * t42 + t82 * t87) + m(6) * (t83 * t136 + t82 * t97) + m(5) * t130 * pkin(3) + t122; t87 * t164 + 0.2e1 * t97 * t88 + Ifges(4,3) + 0.2e1 * t124 + (-t126 * t42 - t41 * t84) * t143 + t96 * t127 + m(7) * (t41 ^ 2 + t42 ^ 2 + t87 ^ 2) + m(6) * (t145 * t96 ^ 2 + t97 ^ 2) + m(5) * (t110 ^ 2 + t114 ^ 2) * pkin(3) ^ 2 + t134; t131 * pkin(9) + m(7) * (t25 * t98 + t3 * t58 + t4 * t59) + t121 + m(6) * (-pkin(4) * t52 + t132 * pkin(9)) - mrSges(6,3) * t151 + t98 * t11 + t58 * t19 + t59 * t18 - pkin(4) * t43; m(7) * (-t34 * t58 - t36 * t59 + t82 * t98) + m(6) * (-pkin(4) * t82 + t83 * t137) + t122; (t97 - pkin(4)) * t88 + (t87 + t98) * t49 + t124 + m(7) * (t41 * t58 + t42 * t59 + t87 * t98) + m(6) * (-pkin(4) * t97 + pkin(9) * t136) + ((-t41 - t58) * t84 - (t42 + t59) * t126) * mrSges(7,3) + (t136 + t137) * mrSges(6,3) + t134; -0.2e1 * pkin(4) * t88 + t98 * t164 + (-t126 * t59 - t58 * t84) * t143 + pkin(9) * t127 + m(7) * (t58 ^ 2 + t59 ^ 2 + t98 ^ 2) + m(6) * (t145 * pkin(9) ^ 2 + pkin(4) ^ 2) + t134; Ifges(6,6) * t150 + t14 * mrSges(6,1) - t15 * mrSges(6,2) + (m(7) * (t108 * t4 + t112 * t3) + t108 * t18 + t112 * t19) * pkin(5) + t125 + t169; t133 * t83 + m(7) * (-t108 * t36 - t112 * t34) * pkin(5) + t140; t133 * t96 + (m(7) * (t108 * t42 + t112 * t41) - t141) * pkin(5) + t129 + t168; t133 * pkin(9) + (m(7) * (t108 * t59 + t112 * t58) - t141) * pkin(5) + t128 + t168; Ifges(6,3) + Ifges(7,3) + m(7) * (t108 ^ 2 + t112 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t123; t125; t140; t129; t128; Ifges(7,3) + t123; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
