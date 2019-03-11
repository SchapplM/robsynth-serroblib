% Calculate Gravitation load on the joints for
% S6RRRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 23:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR12_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:29:16
% EndTime: 2019-03-09 23:29:22
% DurationCPUTime: 2.14s
% Computational Cost: add. (1445->195), mult. (3469->286), div. (0->0), fcn. (4347->16), ass. (0->89)
t94 = sin(qJ(6));
t97 = cos(qJ(6));
t180 = m(7) * pkin(5) + t97 * mrSges(7,1) - t94 * mrSges(7,2) + mrSges(6,1);
t179 = -m(7) * pkin(12) + mrSges(6,2) - mrSges(7,3);
t143 = cos(pkin(7));
t92 = sin(pkin(6));
t135 = t92 * t143;
t144 = cos(pkin(6));
t168 = cos(qJ(2));
t121 = t144 * t168;
t165 = sin(qJ(2));
t166 = sin(qJ(1));
t99 = cos(qJ(1));
t69 = -t121 * t99 + t165 * t166;
t91 = sin(pkin(7));
t186 = t99 * t135 - t69 * t91;
t96 = sin(qJ(3));
t134 = t96 * t143;
t150 = t92 * t99;
t142 = t91 * t150;
t167 = cos(qJ(3));
t120 = t144 * t165;
t70 = t120 * t99 + t166 * t168;
t26 = -t134 * t69 - t96 * t142 + t167 * t70;
t95 = sin(qJ(4));
t98 = cos(qJ(4));
t189 = t186 * t98 + t26 * t95;
t161 = t186 * t95;
t188 = -t26 * t98 + t161;
t172 = -m(5) * pkin(3) - t98 * mrSges(5,1) + t95 * mrSges(5,2) - mrSges(4,1);
t106 = t121 * t166 + t165 * t99;
t53 = t106 * t91 + t166 * t135;
t136 = t91 * t144;
t50 = t96 * t136 + (t134 * t168 + t165 * t167) * t92;
t139 = t92 * t168;
t68 = -t139 * t91 + t143 * t144;
t185 = -t50 * t95 + t68 * t98;
t138 = t92 * t166;
t178 = t106 * t143 - t91 * t138;
t71 = -t120 * t166 + t168 * t99;
t30 = t71 * t167 - t178 * t96;
t9 = -t30 * t95 + t53 * t98;
t90 = qJ(4) + pkin(13);
t87 = sin(t90);
t88 = cos(t90);
t4 = -t186 * t87 + t26 * t88;
t183 = -t186 * t88 - t26 * t87;
t182 = -m(4) - m(5);
t181 = m(6) + m(7);
t115 = -m(5) * pkin(11) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t177 = -t94 * mrSges(7,1) - t97 * mrSges(7,2) + t115;
t174 = -t179 * t87 + t180 * t88 - t172;
t173 = mrSges(3,2) + (-t95 * mrSges(5,1) - t98 * mrSges(5,2) - mrSges(4,3)) * t91;
t169 = pkin(10) * t91;
t159 = t53 * t95;
t155 = t87 * t91;
t154 = t88 * t91;
t137 = t92 * t165;
t132 = t91 * t137;
t146 = pkin(2) * t139 + pkin(10) * t132;
t145 = t99 * pkin(1) + pkin(9) * t138;
t124 = (pkin(4) * t95 + pkin(10)) * t91;
t123 = -pkin(1) * t166 + pkin(9) * t150;
t119 = t143 * t167;
t118 = t143 * t165;
t114 = t95 * t132;
t107 = -t70 * pkin(2) + t186 * pkin(10) + t123;
t25 = t119 * t69 + t142 * t167 + t70 * t96;
t101 = t71 * pkin(2) + t53 * pkin(10) + t145;
t29 = t178 * t167 + t71 * t96;
t86 = pkin(4) * t98 + pkin(3);
t93 = -qJ(5) - pkin(11);
t100 = pkin(4) * t159 - t29 * t93 + t30 * t86 + t101;
t66 = t106 * pkin(2);
t64 = t69 * pkin(2);
t62 = (-t118 * t96 + t167 * t168) * t92;
t61 = (t118 * t167 + t168 * t96) * t92;
t49 = -t119 * t139 - t136 * t167 + t137 * t96;
t40 = -t106 * t167 - t134 * t71;
t39 = -t106 * t96 + t119 * t71;
t38 = -t134 * t70 - t167 * t69;
t37 = t119 * t70 - t69 * t96;
t20 = t50 * t88 + t68 * t87;
t10 = t30 * t98 + t159;
t8 = t30 * t88 + t53 * t87;
t7 = t30 * t87 - t53 * t88;
t2 = t29 * t94 + t8 * t97;
t1 = t29 * t97 - t8 * t94;
t3 = [(-t99 * mrSges(2,1) + t166 * mrSges(2,2) - m(3) * t145 - t71 * mrSges(3,1) + t106 * mrSges(3,2) - mrSges(3,3) * t138 - m(4) * t101 - t30 * mrSges(4,1) - t53 * mrSges(4,3) - m(5) * (t30 * pkin(3) + t101) - t10 * mrSges(5,1) - t9 * mrSges(5,2) - m(6) * t100 - t8 * mrSges(6,1) - m(7) * (t8 * pkin(5) + t100) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t179 * t7 + t115 * t29) * g(2) + (t166 * mrSges(2,1) - m(3) * t123 + t70 * mrSges(3,1) - t69 * mrSges(3,2) - m(4) * t107 + t26 * mrSges(4,1) - t186 * mrSges(4,3) - m(5) * (-pkin(3) * t26 + t107) - t188 * mrSges(5,1) - t189 * mrSges(5,2) + (-mrSges(3,3) * t92 + mrSges(2,2)) * t99 + t179 * t183 + t180 * t4 - t177 * t25 - t181 * (pkin(4) * t161 + t25 * t93 - t26 * t86 + t107)) * g(1) (-(mrSges(3,1) * t168 - mrSges(3,2) * t165) * t92 - m(4) * t146 - t62 * mrSges(4,1) - mrSges(4,3) * t132 - m(5) * (pkin(3) * t62 + t146) - (t62 * t98 + t114) * mrSges(5,1) - (t132 * t98 - t62 * t95) * mrSges(5,2) - t181 * (pkin(4) * t114 - t61 * t93 + t62 * t86 + t146) + t179 * (-t132 * t88 + t62 * t87) - t180 * (t132 * t87 + t62 * t88) + t177 * t61) * g(3) + (mrSges(3,1) * t69 + t182 * (t169 * t70 - t64) + t172 * t38 + t173 * t70 - t181 * (t124 * t70 - t37 * t93 + t38 * t86 - t64) + t179 * (-t154 * t70 + t38 * t87) - t180 * (t155 * t70 + t38 * t88) + t177 * t37) * g(2) + (t106 * mrSges(3,1) + t182 * (t169 * t71 - t66) + t172 * t40 + t173 * t71 - t181 * (t124 * t71 - t39 * t93 + t40 * t86 - t66) + t179 * (-t154 * t71 + t40 * t87) - t180 * (t155 * t71 + t40 * t88) + t177 * t39) * g(1) (-t181 * (-t49 * t86 - t50 * t93) + t177 * t50 + t174 * t49) * g(3) + (-t181 * (-t25 * t86 - t26 * t93) + t177 * t26 + t174 * t25) * g(2) + (-t181 * (-t29 * t86 - t30 * t93) + t177 * t30 + t174 * t29) * g(1) (-t185 * mrSges(5,1) - (-t50 * t98 - t68 * t95) * mrSges(5,2) + t179 * t20 - t180 * (-t50 * t87 + t68 * t88)) * g(3) + (mrSges(5,1) * t189 - mrSges(5,2) * t188 + t179 * t4 - t180 * t183) * g(2) + (-t9 * mrSges(5,1) + t10 * mrSges(5,2) + t179 * t8 + t180 * t7) * g(1) + (-g(1) * t9 + g(2) * t189 - g(3) * t185) * t181 * pkin(4), t181 * (-g(1) * t29 - g(2) * t25 - g(3) * t49) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t25 * t97 - t4 * t94) * mrSges(7,1) + (-t25 * t94 - t4 * t97) * mrSges(7,2)) - g(3) * ((-t20 * t94 + t49 * t97) * mrSges(7,1) + (-t20 * t97 - t49 * t94) * mrSges(7,2))];
taug  = t3(:);
