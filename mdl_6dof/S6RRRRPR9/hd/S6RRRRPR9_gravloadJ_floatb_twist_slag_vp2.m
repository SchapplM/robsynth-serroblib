% Calculate Gravitation load on the joints for
% S6RRRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:47:49
% EndTime: 2019-03-09 22:47:53
% DurationCPUTime: 1.77s
% Computational Cost: add. (955->155), mult. (1538->195), div. (0->0), fcn. (1782->14), ass. (0->82)
t74 = sin(pkin(12));
t76 = cos(pkin(12));
t195 = -t76 * mrSges(6,2) - (m(7) * pkin(5) + mrSges(6,1)) * t74;
t185 = -mrSges(6,3) - mrSges(7,3);
t72 = pkin(12) + qJ(6);
t67 = sin(t72);
t68 = cos(t72);
t193 = -t68 * mrSges(7,1) + t67 * mrSges(7,2);
t191 = t76 * mrSges(6,1) - t74 * mrSges(6,2) + mrSges(5,1);
t65 = pkin(5) * t76 + pkin(4);
t85 = -m(6) * pkin(4) - m(7) * t65 - t191;
t183 = t85 + t193;
t73 = qJ(3) + qJ(4);
t70 = cos(t73);
t192 = -m(4) * pkin(2) + t183 * t70 - mrSges(3,1);
t77 = -pkin(11) - qJ(5);
t190 = -m(6) * qJ(5) + m(7) * t77 + mrSges(5,2);
t188 = t191 - t193;
t88 = t190 + t185;
t163 = m(6) + m(7);
t186 = m(5) + t163;
t182 = -m(4) * pkin(9) + mrSges(3,2) - mrSges(5,3);
t181 = mrSges(5,2) + t185;
t97 = mrSges(7,1) * t67 + mrSges(7,2) * t68;
t179 = t195 - t97;
t122 = cos(pkin(6));
t75 = sin(pkin(6));
t79 = sin(qJ(2));
t138 = t75 * t79;
t78 = sin(qJ(3));
t81 = cos(qJ(3));
t171 = t122 * t81 - t78 * t138;
t136 = t75 * t81;
t108 = t79 * t122;
t143 = cos(qJ(2));
t144 = cos(qJ(1));
t80 = sin(qJ(1));
t51 = -t108 * t80 + t143 * t144;
t27 = t80 * t136 - t51 * t78;
t114 = t75 * t144;
t49 = t108 * t144 + t143 * t80;
t69 = sin(t73);
t21 = t114 * t70 + t49 * t69;
t22 = -t69 * t114 + t49 * t70;
t170 = t181 * t22 + t188 * t21;
t137 = t75 * t80;
t25 = -t137 * t70 + t51 * t69;
t26 = t137 * t69 + t51 * t70;
t169 = t181 * t26 + t188 * t25;
t42 = -t122 * t70 + t138 * t69;
t43 = t122 * t69 + t138 * t70;
t168 = t181 * t43 + t188 * t42;
t166 = -mrSges(4,3) + t182;
t159 = t81 * mrSges(4,1) - t78 * mrSges(4,2);
t165 = -t88 * t69 + t159 - t192;
t164 = t166 + t179;
t153 = t166 + t195;
t152 = -t21 * t65 - t22 * t77;
t145 = -t25 * t65 - t26 * t77;
t129 = -t42 * t65 - t43 * t77;
t123 = t144 * pkin(1) + pkin(8) * t137;
t121 = t78 * t137;
t113 = t75 * t143;
t111 = -pkin(1) * t80 + pkin(8) * t114;
t59 = t78 * t114;
t109 = -t49 * t81 + t59;
t106 = -t21 * pkin(4) + t22 * qJ(5);
t105 = -t25 * pkin(4) + qJ(5) * t26;
t104 = -t42 * pkin(4) + qJ(5) * t43;
t102 = t27 * pkin(3);
t100 = t122 * t143;
t95 = t171 * pkin(3);
t90 = t114 * t81 + t49 * t78;
t86 = t90 * pkin(3);
t82 = -pkin(10) - pkin(9);
t66 = pkin(3) * t81 + pkin(2);
t50 = t100 * t80 + t144 * t79;
t48 = -t100 * t144 + t79 * t80;
t28 = t51 * t81 + t121;
t2 = t26 * t68 + t50 * t67;
t1 = -t26 * t67 + t50 * t68;
t3 = [(-t144 * mrSges(2,1) - m(3) * t123 - t51 * mrSges(3,1) - m(4) * (pkin(2) * t51 + t123) - t28 * mrSges(4,1) - t27 * mrSges(4,2) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + (-mrSges(3,3) * t75 + mrSges(2,2)) * t80 + t85 * t26 + t153 * t50 + t88 * t25 - t186 * (pkin(3) * t121 - t50 * t82 + t51 * t66 + t123)) * g(2) + (t80 * mrSges(2,1) + t144 * mrSges(2,2) - m(3) * t111 + t49 * mrSges(3,1) - mrSges(3,3) * t114 - m(4) * (-pkin(2) * t49 + t111) - t109 * mrSges(4,1) - t90 * mrSges(4,2) - t183 * t22 + (-t153 + t97) * t48 - t88 * t21 + t186 * (-pkin(3) * t59 - t48 * t82 + t49 * t66 - t111)) * g(1) (-t186 * (-t48 * t66 - t49 * t82) + t164 * t49 + t165 * t48) * g(2) + (-t186 * (-t50 * t66 - t51 * t82) + t164 * t51 + t165 * t50) * g(1) + (-mrSges(4,3) * t138 + ((t186 * t82 + t179 + t182) * t79 + (t190 * t69 + t192) * t143) * t75 + (t185 * t69 - t186 * t66 - t159) * t113) * g(3) (-t171 * mrSges(4,1) - (-t122 * t78 - t136 * t79) * mrSges(4,2) - m(5) * t95 - m(6) * (t104 + t95) - m(7) * (t95 + t129) + t168) * g(3) + (t90 * mrSges(4,1) - t109 * mrSges(4,2) + m(5) * t86 - m(6) * (t106 - t86) - m(7) * (-t86 + t152) + t170) * g(2) + (-mrSges(4,1) * t27 + mrSges(4,2) * t28 - m(5) * t102 - m(6) * (t102 + t105) - m(7) * (t102 + t145) + t169) * g(1) (-m(6) * t104 - m(7) * t129 + t168) * g(3) + (-m(6) * t106 - m(7) * t152 + t170) * g(2) + (-m(6) * t105 - m(7) * t145 + t169) * g(1), t163 * (-g(1) * t25 - g(2) * t21 - g(3) * t42) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t22 * t67 + t48 * t68) * mrSges(7,1) + (-t22 * t68 - t48 * t67) * mrSges(7,2)) - g(3) * ((-t113 * t68 - t43 * t67) * mrSges(7,1) + (t113 * t67 - t43 * t68) * mrSges(7,2))];
taug  = t3(:);
