% Calculate Gravitation load on the joints for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:14:12
% EndTime: 2019-03-09 01:14:18
% DurationCPUTime: 2.11s
% Computational Cost: add. (2401->196), mult. (6932->312), div. (0->0), fcn. (9003->18), ass. (0->109)
t187 = m(6) + m(7);
t94 = sin(qJ(6));
t98 = cos(qJ(6));
t189 = m(7) * pkin(5) + t98 * mrSges(7,1) - t94 * mrSges(7,2) + mrSges(6,1);
t177 = -m(7) * pkin(13) + mrSges(6,2) - mrSges(7,3);
t156 = sin(pkin(6));
t168 = sin(qJ(2));
t141 = t156 * t168;
t93 = sin(pkin(7));
t134 = t93 * t141;
t158 = cos(pkin(8));
t100 = cos(qJ(3));
t159 = cos(pkin(7));
t137 = t159 * t156;
t125 = t168 * t137;
t170 = cos(qJ(2));
t142 = t170 * t156;
t97 = sin(qJ(3));
t80 = -t100 * t125 - t142 * t97;
t92 = sin(pkin(8));
t70 = t158 * t134 - t80 * t92;
t174 = -t94 * mrSges(7,1) - t98 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t95 = sin(qJ(5));
t99 = cos(qJ(5));
t188 = pkin(4) * t187 - t177 * t95 + t189 * t99 + mrSges(5,1);
t185 = -t93 * mrSges(4,3) + mrSges(3,2);
t184 = -t92 * mrSges(5,3) + mrSges(4,2);
t155 = sin(pkin(14));
t160 = cos(pkin(6));
t138 = t160 * t155;
t157 = cos(pkin(14));
t120 = t138 * t170 + t157 * t168;
t135 = t156 * t155;
t108 = -t120 * t159 + t135 * t93;
t85 = -t138 * t168 + t157 * t170;
t103 = -t100 * t108 + t85 * t97;
t109 = t120 * t93 + t135 * t159;
t183 = t103 * t158 - t109 * t92;
t139 = t160 * t157;
t119 = -t139 * t170 + t155 * t168;
t136 = t157 * t156;
t180 = -t119 * t159 - t93 * t136;
t84 = t139 * t168 + t155 * t170;
t104 = -t100 * t180 + t84 * t97;
t107 = t119 * t93 - t136 * t159;
t182 = t104 * t158 - t107 * t92;
t122 = t137 * t170 + t160 * t93;
t111 = -t100 * t122 + t141 * t97;
t121 = -t142 * t93 + t159 * t160;
t181 = t111 * t158 - t121 * t92;
t175 = -t187 * pkin(12) + t174;
t172 = pkin(10) * t93;
t171 = pkin(11) * t92;
t169 = cos(qJ(4));
t165 = t92 * t93;
t164 = t92 * t95;
t163 = t92 * t99;
t161 = pkin(2) * t142 + pkin(10) * t134;
t153 = t93 * t158;
t96 = sin(qJ(4));
t152 = t96 * t158;
t151 = t97 * t159;
t150 = t100 * t159;
t149 = t169 * t165;
t148 = -t119 * pkin(2) + t172 * t84;
t147 = -t120 * pkin(2) + t172 * t85;
t59 = t84 * t100 + t180 * t97;
t146 = -t104 * pkin(3) + t171 * t59;
t60 = t85 * t100 + t108 * t97;
t145 = -t103 * pkin(3) + t171 * t60;
t76 = t100 * t141 + t122 * t97;
t144 = -t111 * pkin(3) + t171 * t76;
t143 = t158 * t169;
t132 = t92 * t134;
t81 = t100 * t142 - t125 * t97;
t131 = t81 * pkin(3) + t70 * pkin(11) + t161;
t65 = t119 * t97 - t150 * t84;
t47 = t153 * t84 - t65 * t92;
t67 = t120 * t97 - t150 * t85;
t48 = t153 * t85 - t67 * t92;
t66 = -t100 * t119 - t151 * t84;
t116 = t66 * pkin(3) + pkin(11) * t47 + t148;
t68 = -t100 * t120 - t151 * t85;
t115 = t68 * pkin(3) + pkin(11) * t48 + t147;
t58 = t111 * t92 + t121 * t158;
t51 = t81 * t169 + (t158 * t80 + t132) * t96;
t50 = -t132 * t169 - t143 * t80 + t81 * t96;
t46 = -t111 * t169 - t152 * t76;
t45 = -t111 * t96 + t143 * t76;
t41 = t103 * t92 + t109 * t158;
t40 = t104 * t92 + t107 * t158;
t39 = t76 * t169 - t181 * t96;
t38 = t169 * t181 + t76 * t96;
t32 = -t103 * t169 - t152 * t60;
t31 = -t103 * t96 + t143 * t60;
t30 = -t104 * t169 - t152 * t59;
t29 = -t104 * t96 + t143 * t59;
t28 = t68 * t169 + (t158 * t67 + t165 * t85) * t96;
t27 = -t143 * t67 - t149 * t85 + t68 * t96;
t26 = t66 * t169 + (t158 * t65 + t165 * t84) * t96;
t25 = -t143 * t65 - t149 * t84 + t66 * t96;
t20 = t60 * t169 - t183 * t96;
t19 = t169 * t183 + t60 * t96;
t18 = t59 * t169 - t182 * t96;
t17 = t169 * t182 + t59 * t96;
t16 = t39 * t99 + t58 * t95;
t4 = t20 * t99 + t41 * t95;
t2 = t18 * t99 + t40 * t95;
t1 = [(-m(2) - m(3) - m(4) - m(5) - t187) * g(3) (-m(4) * t161 - m(5) * t131 - mrSges(3,1) * t142 - t81 * mrSges(4,1) - t51 * mrSges(5,1) + mrSges(3,2) * t141 - t80 * mrSges(4,2) - mrSges(4,3) * t134 - t70 * mrSges(5,3) - t187 * (t51 * pkin(4) + pkin(12) * t50 + t131) - t189 * (t51 * t99 + t70 * t95) + t174 * t50 + t177 * (t51 * t95 - t70 * t99)) * g(3) + (-m(4) * t148 - m(5) * t116 + mrSges(3,1) * t119 - t66 * mrSges(4,1) - t26 * mrSges(5,1) - t65 * mrSges(4,2) - t47 * mrSges(5,3) + t185 * t84 - t187 * (t26 * pkin(4) + t25 * pkin(12) + t116) + t177 * (t26 * t95 - t47 * t99) - t189 * (t26 * t99 + t47 * t95) + t174 * t25) * g(2) + (-m(4) * t147 - m(5) * t115 + mrSges(3,1) * t120 - t68 * mrSges(4,1) - t28 * mrSges(5,1) - t67 * mrSges(4,2) - t48 * mrSges(5,3) + t185 * t85 - t187 * (t28 * pkin(4) + t27 * pkin(12) + t115) + t177 * (t28 * t95 - t48 * t99) - t189 * (t28 * t99 + t48 * t95) + t174 * t27) * g(1) (-m(5) * t144 + mrSges(4,1) * t111 - t46 * mrSges(5,1) + t184 * t76 - t187 * (t46 * pkin(4) + pkin(12) * t45 + t144) - t189 * (t164 * t76 + t46 * t99) + t174 * t45 + t177 * (-t163 * t76 + t46 * t95)) * g(3) + (-m(5) * t146 + mrSges(4,1) * t104 - t30 * mrSges(5,1) - t187 * (t30 * pkin(4) + pkin(12) * t29 + t146) + t177 * (-t163 * t59 + t30 * t95) + t184 * t59 - t189 * (t164 * t59 + t30 * t99) + t174 * t29) * g(2) + (-m(5) * t145 + mrSges(4,1) * t103 - t32 * mrSges(5,1) + t184 * t60 - t187 * (t32 * pkin(4) + pkin(12) * t31 + t145) - t189 * (t164 * t60 + t32 * t99) + t174 * t31 + t177 * (-t163 * t60 + t32 * t95)) * g(1) (t175 * t39 + t188 * t38) * g(3) + (t188 * t17 + t175 * t18) * g(2) + (t175 * t20 + t188 * t19) * g(1) (t177 * t16 - t189 * (-t39 * t95 + t58 * t99)) * g(3) + (t177 * t2 - t189 * (-t18 * t95 + t40 * t99)) * g(2) + (t177 * t4 - t189 * (-t20 * t95 + t41 * t99)) * g(1), -g(1) * ((t19 * t98 - t4 * t94) * mrSges(7,1) + (-t19 * t94 - t4 * t98) * mrSges(7,2)) - g(2) * ((t17 * t98 - t2 * t94) * mrSges(7,1) + (-t17 * t94 - t2 * t98) * mrSges(7,2)) - g(3) * ((-t16 * t94 + t38 * t98) * mrSges(7,1) + (-t16 * t98 - t38 * t94) * mrSges(7,2))];
taug  = t1(:);
