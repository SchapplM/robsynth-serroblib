% Calculate joint inertia matrix for
% S6RPRRRR2
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:33
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:32:57
% EndTime: 2018-11-23 16:32:58
% DurationCPUTime: 1.13s
% Computational Cost: add. (1952->269), mult. (3631->384), div. (0->0), fcn. (3809->10), ass. (0->108)
t109 = sin(qJ(5));
t113 = cos(qJ(5));
t140 = t109 ^ 2 + t113 ^ 2;
t166 = mrSges(6,3) * t140;
t115 = cos(qJ(3));
t165 = t115 ^ 2;
t108 = sin(qJ(6));
t112 = cos(qJ(6));
t82 = t108 * t113 + t109 * t112;
t110 = sin(qJ(4));
t111 = sin(qJ(3));
t114 = cos(qJ(4));
t83 = t110 * t115 + t111 * t114;
t38 = t82 * t83;
t79 = -t108 * t109 + t112 * t113;
t39 = t79 * t83;
t15 = t38 * mrSges(7,1) + t39 * mrSges(7,2);
t130 = mrSges(6,1) * t109 + mrSges(6,2) * t113;
t47 = t130 * t83;
t164 = t15 + t47;
t143 = t113 * t83;
t80 = t110 * t111 - t114 * t115;
t163 = Ifges(6,5) * t143 + Ifges(6,3) * t80;
t141 = Ifges(6,5) * t109 + Ifges(6,6) * t113;
t151 = t79 * mrSges(7,3);
t162 = t108 * pkin(5) * t151 + t141;
t106 = sin(pkin(11));
t93 = pkin(1) * t106 + pkin(7);
t155 = pkin(8) + t93;
t135 = t155 * t111;
t65 = t155 * t115;
t42 = t110 * t65 + t114 * t135;
t161 = t42 ^ 2;
t75 = t80 ^ 2;
t160 = 0.2e1 * t42;
t50 = -mrSges(7,1) * t79 + mrSges(7,2) * t82;
t159 = 0.2e1 * t50;
t107 = cos(pkin(11));
t94 = -pkin(1) * t107 - pkin(2);
t84 = -pkin(3) * t115 + t94;
t158 = 0.2e1 * t84;
t156 = pkin(4) * t80;
t153 = pkin(3) * t114;
t152 = t42 * t80;
t150 = t82 * mrSges(7,3);
t37 = -pkin(9) * t83 + t156 + t84;
t44 = -t110 * t135 + t114 * t65;
t14 = t109 * t37 + t113 * t44;
t149 = Ifges(7,5) * t82 + Ifges(7,6) * t79;
t148 = Ifges(6,4) * t109;
t147 = Ifges(6,4) * t113;
t13 = -t109 * t44 + t113 * t37;
t146 = t109 * t13;
t145 = t109 * t83;
t144 = t113 * t14;
t96 = pkin(3) * t110 + pkin(9);
t142 = t113 * t96;
t139 = t111 ^ 2 + t165;
t138 = 0.2e1 * mrSges(7,3);
t137 = Ifges(7,5) * t39 - Ifges(7,6) * t38 + Ifges(7,3) * t80;
t136 = t112 * t150;
t98 = -pkin(5) * t113 - pkin(4);
t134 = t140 * pkin(9);
t133 = t140 * t96;
t51 = Ifges(7,4) * t82 + Ifges(7,2) * t79;
t52 = Ifges(7,1) * t82 + Ifges(7,4) * t79;
t87 = Ifges(6,2) * t113 + t148;
t88 = Ifges(6,1) * t109 + t147;
t132 = t109 * t88 + t113 * t87 + t79 * t51 + t82 * t52 + Ifges(5,3);
t131 = -t115 * mrSges(4,1) + t111 * mrSges(4,2);
t129 = t144 - t146;
t48 = -mrSges(6,2) * t80 - mrSges(6,3) * t145;
t49 = mrSges(6,1) * t80 - mrSges(6,3) * t143;
t128 = -t109 * t49 + t113 * t48;
t73 = (-pkin(10) - t96) * t109;
t101 = t113 * pkin(10);
t74 = t101 + t142;
t45 = -t108 * t74 + t112 * t73;
t46 = t108 * t73 + t112 * t74;
t127 = t45 * mrSges(7,1) - t46 * mrSges(7,2) + t149;
t89 = (-pkin(10) - pkin(9)) * t109;
t90 = pkin(9) * t113 + t101;
t56 = -t108 * t90 + t112 * t89;
t57 = t108 * t89 + t112 * t90;
t126 = t56 * mrSges(7,1) - t57 * mrSges(7,2) + t149;
t125 = 0.2e1 * t166;
t5 = pkin(5) * t80 - pkin(10) * t143 + t13;
t8 = -pkin(10) * t145 + t14;
t3 = -t108 * t8 + t112 * t5;
t4 = t108 * t5 + t112 * t8;
t124 = t3 * mrSges(7,1) - t4 * mrSges(7,2) + t137;
t123 = (mrSges(5,1) * t114 - mrSges(5,2) * t110) * pkin(3);
t122 = (mrSges(7,1) * t112 - mrSges(7,2) * t108) * pkin(5);
t66 = t80 * mrSges(5,1);
t86 = -mrSges(6,1) * t113 + mrSges(6,2) * t109;
t121 = t38 * t150 + t39 * t151 - t66 + (t50 + t86) * t80 + (-mrSges(5,2) + t166) * t83;
t10 = Ifges(7,1) * t39 - Ifges(7,4) * t38 + Ifges(7,5) * t80;
t22 = pkin(5) * t145 + t42;
t24 = Ifges(6,6) * t80 + (-Ifges(6,2) * t109 + t147) * t83;
t25 = Ifges(6,5) * t80 + (Ifges(6,1) * t113 - t148) * t83;
t9 = Ifges(7,4) * t39 - Ifges(7,2) * t38 + Ifges(7,6) * t80;
t120 = -t44 * mrSges(5,2) - t3 * t150 + t4 * t151 + mrSges(6,3) * t144 + t22 * t50 - t38 * t51 / 0.2e1 + t39 * t52 / 0.2e1 + t109 * t25 / 0.2e1 + t113 * t24 / 0.2e1 - t87 * t145 / 0.2e1 + t88 * t143 / 0.2e1 + t79 * t9 / 0.2e1 + t82 * t10 / 0.2e1 - Ifges(5,6) * t80 + Ifges(5,5) * t83 + (t86 - mrSges(5,1)) * t42 + (t149 + t141) * t80 / 0.2e1;
t97 = -pkin(4) - t153;
t85 = t98 - t153;
t76 = t83 ^ 2;
t19 = mrSges(7,1) * t80 - mrSges(7,3) * t39;
t18 = -mrSges(7,2) * t80 - mrSges(7,3) * t38;
t1 = [0.2e1 * t94 * t131 + Ifges(4,2) * t165 + t66 * t158 + t47 * t160 + 0.2e1 * t14 * t48 + 0.2e1 * t13 * t49 - t38 * t9 + t39 * t10 + 0.2e1 * t4 * t18 + 0.2e1 * t3 * t19 + 0.2e1 * t22 * t15 + Ifges(2,3) + Ifges(3,3) + (mrSges(5,2) * t158 + mrSges(5,3) * t160 + Ifges(5,1) * t83 - t109 * t24 + t113 * t25) * t83 + (-0.2e1 * t44 * mrSges(5,3) + Ifges(5,2) * t80 + (-Ifges(6,6) * t109 - (2 * Ifges(5,4))) * t83 + t137 + t163) * t80 + m(4) * (t139 * t93 ^ 2 + t94 ^ 2) + m(5) * (t44 ^ 2 + t84 ^ 2 + t161) + m(6) * (t13 ^ 2 + t14 ^ 2 + t161) + m(7) * (t22 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(3) * (t106 ^ 2 + t107 ^ 2) * pkin(1) ^ 2 + (Ifges(4,1) * t111 + 0.2e1 * Ifges(4,4) * t115) * t111 + 0.2e1 * (t107 * mrSges(3,1) - t106 * mrSges(3,2)) * pkin(1) + 0.2e1 * t139 * t93 * mrSges(4,3); t39 * t18 - t38 * t19 + t128 * t83 + t164 * t80 + m(7) * (t22 * t80 - t3 * t38 + t39 * t4) + m(6) * (t129 * t83 + t152) + m(5) * (t44 * t83 + t152); m(3) + m(7) * (t38 ^ 2 + t39 ^ 2 + t75) + m(6) * (t140 * t76 + t75) + m(5) * (t76 + t75) + m(4) * t139; t120 + (m(5) * (t110 * t44 - t114 * t42) + (-t110 * t80 - t114 * t83) * mrSges(5,3)) * pkin(3) + m(6) * (t129 * t96 + t42 * t97) + (-t13 * mrSges(6,3) - t96 * t49) * t109 + m(7) * (t22 * t85 + t3 * t45 + t4 * t46) + t48 * t142 + (-t93 * mrSges(4,2) + Ifges(4,6)) * t115 + t97 * t47 + t85 * t15 + (-t93 * mrSges(4,1) + Ifges(4,5)) * t111 + t45 * t19 + t46 * t18; m(7) * (-t38 * t45 + t39 * t46 + t80 * t85) + m(6) * (t83 * t133 + t80 * t97) + m(5) * (t110 * t83 - t114 * t80) * pkin(3) + t121 - t131; t85 * t159 + 0.2e1 * t97 * t86 + Ifges(4,3) + 0.2e1 * t123 + (-t45 * t82 + t46 * t79) * t138 + t96 * t125 + m(7) * (t45 ^ 2 + t46 ^ 2 + t85 ^ 2) + m(6) * (t140 * t96 ^ 2 + t97 ^ 2) + m(5) * (t110 ^ 2 + t114 ^ 2) * pkin(3) ^ 2 + t132; t128 * pkin(9) + m(7) * (t22 * t98 + t3 * t56 + t4 * t57) + t120 + m(6) * (-pkin(4) * t42 + t129 * pkin(9)) + t98 * t15 - pkin(4) * t47 + t56 * t19 + t57 * t18 - mrSges(6,3) * t146; m(7) * (-t38 * t56 + t39 * t57 + t80 * t98) + m(6) * (t83 * t134 - t156) + t121; (t97 - pkin(4)) * t86 + (t85 + t98) * t50 + t123 + m(7) * (t45 * t56 + t46 * t57 + t85 * t98) + m(6) * (-pkin(4) * t97 + pkin(9) * t133) + ((-t45 - t56) * t82 + (t46 + t57) * t79) * mrSges(7,3) + (t133 + t134) * mrSges(6,3) + t132; -0.2e1 * pkin(4) * t86 + t98 * t159 + (-t56 * t82 + t57 * t79) * t138 + pkin(9) * t125 + m(7) * (t56 ^ 2 + t57 ^ 2 + t98 ^ 2) + m(6) * (t140 * pkin(9) ^ 2 + pkin(4) ^ 2) + t132; -Ifges(6,6) * t145 + t13 * mrSges(6,1) - t14 * mrSges(6,2) + (t112 * t19 + m(7) * (t108 * t4 + t112 * t3) + t108 * t18) * pkin(5) + t124 + t163; m(7) * (t108 * t39 - t112 * t38) * pkin(5) - t164; -t130 * t96 + (m(7) * (t108 * t46 + t112 * t45) - t136) * pkin(5) + t127 + t162; -t130 * pkin(9) + (m(7) * (t108 * t57 + t112 * t56) - t136) * pkin(5) + t126 + t162; Ifges(6,3) + Ifges(7,3) + m(7) * (t108 ^ 2 + t112 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t122; t124; -t15; t127; t126; Ifges(7,3) + t122; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
