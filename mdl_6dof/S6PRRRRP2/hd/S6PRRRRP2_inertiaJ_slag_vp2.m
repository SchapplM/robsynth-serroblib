% Calculate joint inertia matrix for
% S6PRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:01:14
% EndTime: 2019-03-09 00:01:16
% DurationCPUTime: 1.16s
% Computational Cost: add. (1292->243), mult. (2720->330), div. (0->0), fcn. (2779->10), ass. (0->107)
t184 = Ifges(6,1) + Ifges(7,1);
t183 = Ifges(7,4) + Ifges(6,5);
t182 = mrSges(6,3) + mrSges(7,2);
t181 = Ifges(7,2) + Ifges(6,3);
t110 = cos(qJ(5));
t106 = sin(qJ(5));
t151 = Ifges(7,5) * t106;
t153 = Ifges(6,4) * t106;
t107 = sin(qJ(4));
t108 = sin(qJ(3));
t111 = cos(qJ(4));
t112 = cos(qJ(3));
t70 = t107 * t108 - t111 * t112;
t71 = t107 * t112 + t108 * t111;
t180 = (t184 * t110 + t151 - t153) * t71 + t183 * t70;
t77 = -mrSges(6,1) * t110 + mrSges(6,2) * t106;
t179 = t77 - mrSges(5,1);
t150 = Ifges(7,5) * t110;
t152 = Ifges(6,4) * t110;
t178 = t184 * t106 - t150 + t152;
t177 = t106 ^ 2 + t110 ^ 2;
t104 = sin(pkin(6));
t113 = cos(qJ(2));
t142 = t104 * t113;
t105 = cos(pkin(6));
t109 = sin(qJ(2));
t143 = t104 * t109;
t60 = t105 * t112 - t108 * t143;
t61 = t105 * t108 + t112 * t143;
t33 = t107 * t60 + t111 * t61;
t22 = -t106 * t142 + t110 * t33;
t146 = t110 * t22;
t20 = t106 * t33 + t110 * t142;
t176 = t106 * t20 + t146;
t175 = (Ifges(6,6) - Ifges(7,6)) * t110 + t183 * t106;
t91 = -pkin(3) * t112 - pkin(2);
t38 = pkin(4) * t70 - pkin(10) * t71 + t91;
t168 = -pkin(9) - pkin(8);
t139 = t168 * t108;
t83 = t168 * t112;
t52 = t107 * t139 - t111 * t83;
t7 = -t106 * t52 + t110 * t38;
t8 = t106 * t38 + t110 * t52;
t130 = -t106 * t7 + t110 * t8;
t174 = -m(7) * pkin(5) - mrSges(7,1);
t173 = t182 * t177;
t31 = t107 * t61 - t111 * t60;
t30 = t31 ^ 2;
t50 = -t107 * t83 - t111 * t139;
t172 = t50 ^ 2;
t171 = 0.2e1 * t50;
t76 = -mrSges(7,1) * t110 - mrSges(7,3) * t106;
t170 = 0.2e1 * t76;
t103 = t112 ^ 2;
t165 = Ifges(6,6) * t70;
t164 = pkin(3) * t111;
t3 = qJ(6) * t70 + t8;
t162 = t110 * t3;
t160 = t31 * t50;
t148 = t106 * t71;
t39 = -mrSges(6,2) * t70 - mrSges(6,3) * t148;
t42 = -mrSges(7,2) * t148 + mrSges(7,3) * t70;
t159 = t39 + t42;
t145 = t110 * t71;
t40 = mrSges(6,1) * t70 - mrSges(6,3) * t145;
t41 = -t70 * mrSges(7,1) + mrSges(7,2) * t145;
t158 = -t40 + t41;
t89 = pkin(3) * t107 + pkin(10);
t157 = t177 * pkin(10) * t89;
t156 = t177 * t89 ^ 2;
t154 = t177 * pkin(10) ^ 2;
t93 = t106 * mrSges(7,2);
t147 = t106 * t89;
t144 = qJ(6) * t110;
t140 = t108 ^ 2 + t103;
t135 = t176 * pkin(10);
t134 = t89 * t146 + t147 * t20;
t132 = Ifges(7,6) * t148 + t183 * t145 + t181 * t70;
t4 = -pkin(5) * t70 - t7;
t131 = t106 * t4 + t162;
t129 = t106 * mrSges(6,1) + t110 * mrSges(6,2);
t128 = t106 * mrSges(7,1) - t110 * mrSges(7,3);
t127 = pkin(5) * t106 - t144;
t126 = -t108 * t60 + t112 * t61;
t73 = -pkin(5) * t110 - qJ(6) * t106 - pkin(4);
t79 = -Ifges(7,3) * t110 + t151;
t80 = Ifges(6,2) * t110 + t153;
t125 = Ifges(5,3) + (-t79 + t80) * t110 + t178 * t106;
t124 = (mrSges(5,1) * t111 - mrSges(5,2) * t107) * pkin(3);
t123 = t106 * t158 + t110 * t159;
t121 = mrSges(7,2) * t144 - pkin(5) * t93 + t175;
t120 = -t33 * mrSges(5,2) + (t76 + t179) * t31 + t182 * t176;
t119 = 0.2e1 * t173;
t118 = -m(7) * t127 - t128 - t129;
t10 = t127 * t71 + t50;
t25 = Ifges(7,6) * t70 + (Ifges(7,3) * t106 + t150) * t71;
t26 = t165 + (-Ifges(6,2) * t106 + t152) * t71;
t117 = -t52 * mrSges(5,2) + mrSges(7,2) * t162 + Ifges(5,5) * t71 + t10 * t76 + t4 * t93 + t179 * t50 + t180 * t106 / 0.2e1 + (t79 / 0.2e1 - t80 / 0.2e1) * t148 + t178 * t145 / 0.2e1 + (-t25 / 0.2e1 + t26 / 0.2e1) * t110 + t130 * mrSges(6,3) + (-Ifges(5,6) + t175 / 0.2e1) * t70;
t99 = t104 ^ 2;
t90 = -pkin(4) - t164;
t87 = t99 * t113 ^ 2;
t78 = -mrSges(4,1) * t112 + mrSges(4,2) * t108;
t62 = t73 - t164;
t43 = mrSges(5,1) * t70 + mrSges(5,2) * t71;
t37 = t129 * t71;
t36 = t128 * t71;
t1 = [m(2) + m(5) * (t33 ^ 2 + t30 + t87) + m(4) * (t60 ^ 2 + t61 ^ 2 + t87) + m(3) * (t109 ^ 2 * t99 + t105 ^ 2 + t87) + (m(6) + m(7)) * (t20 ^ 2 + t22 ^ 2 + t30); -t33 * t70 * mrSges(5,3) + t159 * t22 + t158 * t20 + t126 * mrSges(4,3) + (t71 * mrSges(5,3) + t36 + t37) * t31 + (-t109 * mrSges(3,2) + (mrSges(3,1) - t43 - t78) * t113) * t104 + m(6) * (-t20 * t7 + t22 * t8 + t160) + m(7) * (t10 * t31 + t20 * t4 + t22 * t3) + m(5) * (-t142 * t91 + t33 * t52 + t160) + m(4) * (pkin(2) * t142 + pkin(8) * t126); Ifges(4,2) * t103 - 0.2e1 * pkin(2) * t78 + 0.2e1 * t10 * t36 + 0.2e1 * t3 * t42 + t37 * t171 + 0.2e1 * t8 * t39 + 0.2e1 * t4 * t41 + 0.2e1 * t7 * t40 + 0.2e1 * t91 * t43 + Ifges(3,3) + (Ifges(4,1) * t108 + 0.2e1 * Ifges(4,4) * t112) * t108 + 0.2e1 * t140 * pkin(8) * mrSges(4,3) + (-0.2e1 * mrSges(5,3) * t52 + Ifges(5,2) * t70 + t132) * t70 + m(4) * (pkin(8) ^ 2 * t140 + pkin(2) ^ 2) + m(5) * (t52 ^ 2 + t91 ^ 2 + t172) + m(6) * (t7 ^ 2 + t8 ^ 2 + t172) + m(7) * (t10 ^ 2 + t3 ^ 2 + t4 ^ 2) + (mrSges(5,3) * t171 + Ifges(5,1) * t71 - 0.2e1 * Ifges(5,4) * t70 + t180 * t110 + (t25 - t26 - t165) * t106) * t71; t60 * mrSges(4,1) - t61 * mrSges(4,2) + m(6) * (t31 * t90 + t134) + m(7) * (t31 * t62 + t134) + m(5) * (t107 * t33 - t111 * t31) * pkin(3) + t120; (-t108 * mrSges(4,1) - t112 * mrSges(4,2)) * pkin(8) + t117 + Ifges(4,6) * t112 + Ifges(4,5) * t108 + t90 * t37 + t62 * t36 + t123 * t89 + m(7) * (t10 * t62 + t131 * t89) + m(6) * (t130 * t89 + t50 * t90) + (m(5) * (t107 * t52 - t111 * t50) + (-t107 * t70 - t111 * t71) * mrSges(5,3)) * pkin(3); t62 * t170 + 0.2e1 * t90 * t77 + Ifges(4,3) + 0.2e1 * t124 + m(7) * (t62 ^ 2 + t156) + m(6) * (t90 ^ 2 + t156) + m(5) * (t107 ^ 2 + t111 ^ 2) * pkin(3) ^ 2 + t119 * t89 + t125; m(6) * (-pkin(4) * t31 + t135) + m(7) * (t31 * t73 + t135) + t120; t123 * pkin(10) + t117 + t73 * t36 - pkin(4) * t37 + m(6) * (-pkin(4) * t50 + pkin(10) * t130) + m(7) * (pkin(10) * t131 + t10 * t73); (t90 - pkin(4)) * t77 + (t62 + t73) * t76 + t124 + m(7) * (t62 * t73 + t157) + m(6) * (-pkin(4) * t90 + t157) + t125 + (pkin(10) + t89) * t173; -0.2e1 * pkin(4) * t77 + t73 * t170 + m(7) * (t73 ^ 2 + t154) + m(6) * (pkin(4) ^ 2 + t154) + t119 * pkin(10) + t125; (m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3)) * t22 + (-mrSges(6,1) + t174) * t20; m(7) * (-pkin(5) * t4 + qJ(6) * t3) + t3 * mrSges(7,3) - t4 * mrSges(7,1) - pkin(5) * t41 + qJ(6) * t42 - Ifges(6,6) * t148 - t8 * mrSges(6,2) + t7 * mrSges(6,1) + t132; t118 * t89 + t121; pkin(10) * t118 + t121; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t181; m(7) * t20; m(7) * t4 + t41; m(7) * t147 + t93; m(7) * t106 * pkin(10) + t93; t174; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
