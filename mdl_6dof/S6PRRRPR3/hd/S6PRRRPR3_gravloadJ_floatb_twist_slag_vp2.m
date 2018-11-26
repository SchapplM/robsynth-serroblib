% Calculate Gravitation load on the joints for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:23
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:23:29
% EndTime: 2018-11-23 15:23:30
% DurationCPUTime: 0.93s
% Computational Cost: add. (1231->125), mult. (1271->160), div. (0->0), fcn. (1230->16), ass. (0->71)
t69 = sin(qJ(6));
t72 = cos(qJ(6));
t146 = mrSges(7,1) * t69 + mrSges(7,2) * t72 - mrSges(5,2) + mrSges(6,3);
t144 = mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t66 = sin(pkin(11));
t67 = sin(pkin(6));
t118 = t66 * t67;
t70 = sin(qJ(3));
t73 = cos(qJ(3));
t108 = cos(pkin(11));
t105 = pkin(6) + qJ(2);
t88 = sin(t105) / 0.2e1;
t106 = pkin(6) - qJ(2);
t92 = sin(t106);
t52 = t88 - t92 / 0.2e1;
t74 = cos(qJ(2));
t84 = t108 * t74 - t66 * t52;
t142 = t73 * t118 - t84 * t70;
t89 = cos(t105) / 0.2e1;
t93 = cos(t106);
t53 = t89 - t93 / 0.2e1;
t68 = cos(pkin(6));
t141 = t53 * t70 + t68 * t73;
t139 = m(6) + m(7);
t65 = qJ(3) + qJ(4);
t63 = sin(t65);
t64 = cos(t65);
t35 = -t53 * t63 - t68 * t64;
t36 = -t53 * t64 + t63 * t68;
t138 = t144 * t35 - t146 * t36;
t19 = -t118 * t64 + t63 * t84;
t20 = t118 * t63 + t64 * t84;
t137 = t144 * t19 - t146 * t20;
t100 = t67 * t108;
t85 = t108 * t52 + t66 * t74;
t17 = t100 * t64 + t63 * t85;
t18 = -t100 * t63 + t64 * t85;
t136 = t144 * t17 - t146 * t18;
t135 = m(4) * pkin(2) + t73 * mrSges(4,1) - t70 * mrSges(4,2) + t144 * t64 + t146 * t63 + mrSges(3,1);
t134 = m(4) * pkin(8) + m(7) * pkin(5) + t72 * mrSges(7,1) - t69 * mrSges(7,2) + mrSges(6,1) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t131 = pkin(10) * t19;
t130 = pkin(10) * t35;
t129 = t17 * pkin(10);
t71 = sin(qJ(2));
t77 = t93 / 0.2e1 + t89;
t39 = -t108 * t77 + t66 * t71;
t126 = t39 * t64;
t42 = t108 * t71 + t66 * t77;
t125 = t42 * t64;
t51 = t88 + t92 / 0.2e1;
t123 = t51 * t64;
t62 = pkin(3) * t73 + pkin(2);
t75 = -pkin(9) - pkin(8);
t114 = -t39 * t62 - t75 * t85;
t113 = -t42 * t62 - t75 * t84;
t110 = t51 * t62 + t53 * t75;
t109 = qJ(5) * t63;
t99 = -t17 * pkin(4) + t18 * qJ(5);
t98 = -t19 * pkin(4) + qJ(5) * t20;
t97 = -t35 * pkin(4) + qJ(5) * t36;
t96 = -pkin(4) * t126 - t39 * t109 + t114;
t95 = -pkin(4) * t125 - t42 * t109 + t113;
t94 = pkin(4) * t123 + t51 * t109 + t110;
t91 = t142 * pkin(3);
t90 = t141 * pkin(3);
t81 = -t100 * t73 - t70 * t85;
t80 = t91 + t98;
t79 = t90 + t97;
t78 = t81 * pkin(3);
t76 = t78 + t99;
t1 = [(-m(2) - m(3) - m(4) - m(5) - t139) * g(3) (-m(5) * t110 - m(6) * t94 - m(7) * (pkin(10) * t123 + t94) + t134 * t53 - t135 * t51) * g(3) + (-m(5) * t114 - m(6) * t96 - m(7) * (-pkin(10) * t126 + t96) - t134 * t85 + t135 * t39) * g(2) + (-m(5) * t113 - m(6) * t95 - m(7) * (-pkin(10) * t125 + t95) - t134 * t84 + t135 * t42) * g(1) (-t141 * mrSges(4,1) - (t53 * t73 - t68 * t70) * mrSges(4,2) - m(5) * t90 - m(6) * t79 - m(7) * (t79 - t130) + t138) * g(3) + (-t81 * mrSges(4,1) - (t100 * t70 - t73 * t85) * mrSges(4,2) - m(5) * t78 - m(6) * t76 - m(7) * (t76 - t129) + t136) * g(2) + (-t142 * mrSges(4,1) - (-t118 * t70 - t73 * t84) * mrSges(4,2) - m(5) * t91 - m(6) * t80 - m(7) * (t80 - t131) + t137) * g(1) (-m(6) * t97 - m(7) * (t97 - t130) + t138) * g(3) + (-m(6) * t99 - m(7) * (t99 - t129) + t136) * g(2) + (-m(6) * t98 - m(7) * (t98 - t131) + t137) * g(1), t139 * (-g(1) * t19 - g(2) * t17 - g(3) * t35) -g(1) * ((t19 * t72 - t42 * t69) * mrSges(7,1) + (-t19 * t69 - t42 * t72) * mrSges(7,2)) - g(2) * ((t17 * t72 - t39 * t69) * mrSges(7,1) + (-t17 * t69 - t39 * t72) * mrSges(7,2)) - g(3) * ((t35 * t72 + t51 * t69) * mrSges(7,1) + (-t35 * t69 + t51 * t72) * mrSges(7,2))];
taug  = t1(:);
