% Calculate Gravitation load on the joints for
% S6RRRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:53:54
% EndTime: 2019-03-09 16:53:58
% DurationCPUTime: 1.50s
% Computational Cost: add. (795->141), mult. (1412->196), div. (0->0), fcn. (1627->12), ass. (0->66)
t123 = -m(4) * pkin(9) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t63 = cos(qJ(5));
t50 = pkin(5) * t63 + pkin(4);
t136 = -m(6) * pkin(4) - m(7) * t50 - mrSges(5,1);
t71 = -m(6) * pkin(10) + m(7) * (-qJ(6) - pkin(10)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t135 = mrSges(6,1) + mrSges(7,1);
t126 = mrSges(6,2) + mrSges(7,2);
t119 = m(7) * pkin(5);
t59 = sin(qJ(5));
t69 = t119 * t59 - t123;
t55 = qJ(3) + pkin(11);
t52 = sin(t55);
t53 = cos(t55);
t60 = sin(qJ(3));
t64 = cos(qJ(3));
t134 = -m(4) * pkin(2) - t64 * mrSges(4,1) + t60 * mrSges(4,2) + t136 * t53 + t52 * t71 - mrSges(3,1);
t125 = m(5) + m(6) + m(7);
t56 = sin(pkin(6));
t61 = sin(qJ(2));
t105 = t56 * t61;
t94 = cos(pkin(6));
t133 = -t60 * t105 + t94 * t64;
t103 = t56 * t64;
t116 = cos(qJ(1));
t65 = cos(qJ(2));
t62 = sin(qJ(1));
t86 = t62 * t94;
t36 = t116 * t65 - t61 * t86;
t19 = t62 * t103 - t36 * t60;
t80 = t94 * t116;
t34 = t61 * t80 + t62 * t65;
t89 = t56 * t116;
t14 = t34 * t53 - t52 * t89;
t33 = t61 * t62 - t65 * t80;
t131 = t14 * t59 - t33 * t63;
t130 = -t14 * t63 - t33 * t59;
t73 = t34 * t60 + t64 * t89;
t124 = -t126 * t59 + t135 * t63 - t136;
t122 = -t119 - t135;
t113 = t34 * t59;
t111 = t36 * t59;
t107 = t53 * t59;
t106 = t53 * t63;
t104 = t56 * t62;
t102 = t56 * t65;
t101 = t59 * t65;
t100 = t63 * t65;
t95 = t116 * pkin(1) + pkin(8) * t104;
t93 = t60 * t104;
t88 = -pkin(1) * t62 + pkin(8) * t89;
t44 = t60 * t89;
t87 = -t34 * t64 + t44;
t35 = t116 * t61 + t65 * t86;
t51 = pkin(3) * t64 + pkin(2);
t58 = -qJ(4) - pkin(9);
t82 = pkin(3) * t93 - t35 * t58 + t36 * t51 + t95;
t18 = t52 * t104 + t36 * t53;
t5 = -t18 * t59 + t35 * t63;
t76 = pkin(3) * t44 + t33 * t58 - t34 * t51 + t88;
t13 = t34 * t52 + t53 * t89;
t28 = t53 * t105 + t94 * t52;
t27 = t52 * t105 - t94 * t53;
t20 = t36 * t64 + t93;
t17 = -t53 * t104 + t36 * t52;
t6 = t18 * t63 + t35 * t59;
t1 = [(-t116 * mrSges(2,1) - m(3) * t95 - t36 * mrSges(3,1) - m(4) * (pkin(2) * t36 + t95) - t20 * mrSges(4,1) - t19 * mrSges(4,2) - m(5) * t82 - t18 * mrSges(5,1) - m(6) * (pkin(4) * t18 + t82) - m(7) * (t18 * t50 + t82) + (-mrSges(3,3) * t56 + mrSges(2,2)) * t62 - t135 * t6 - t126 * t5 - t69 * t35 + t71 * t17) * g(2) + (t62 * mrSges(2,1) + t116 * mrSges(2,2) - m(3) * t88 + t34 * mrSges(3,1) - mrSges(3,3) * t89 - m(4) * (-pkin(2) * t34 + t88) - t87 * mrSges(4,1) - t73 * mrSges(4,2) - m(5) * t76 + t14 * mrSges(5,1) - m(6) * (-pkin(4) * t14 + t76) - m(7) * (-t14 * t50 + t76) - t135 * t130 - t126 * t131 + t69 * t33 - t71 * t13) * g(1) (-t113 * t119 - t125 * (-t33 * t51 - t34 * t58) - t135 * (-t33 * t106 + t113) - t126 * (t33 * t107 + t34 * t63) + t123 * t34 - t134 * t33) * g(2) + (-t111 * t119 - t125 * (-t35 * t51 - t36 * t58) - t126 * (t35 * t107 + t36 * t63) - t135 * (-t35 * t106 + t111) + t123 * t36 - t134 * t35) * g(1) + (-t125 * t51 * t102 + (t134 * t65 + (-t100 * t135 + t126 * t101) * t53 + (t125 * t58 - t126 * t63 - t135 * t59 - t69) * t61) * t56) * g(3) (-t133 * mrSges(4,1) - (-t61 * t103 - t94 * t60) * mrSges(4,2) + t71 * t28 + t124 * t27) * g(3) + (t73 * mrSges(4,1) - t87 * mrSges(4,2) + t124 * t13 + t71 * t14) * g(2) + (-mrSges(4,1) * t19 + mrSges(4,2) * t20 + t124 * t17 + t71 * t18) * g(1) + (-g(1) * t19 + g(2) * t73 - g(3) * t133) * t125 * pkin(3), t125 * (-g(1) * t35 - g(2) * t33 + g(3) * t102) (-t126 * (t56 * t101 - t28 * t63) + t122 * (-t56 * t100 - t28 * t59)) * g(3) + (-t122 * t131 - t126 * t130) * g(2) + (t122 * t5 + t126 * t6) * g(1) (-g(1) * t17 - g(2) * t13 - g(3) * t27) * m(7)];
taug  = t1(:);
