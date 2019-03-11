% Calculate Gravitation load on the joints for
% S6RRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:32:49
% EndTime: 2019-03-09 18:32:53
% DurationCPUTime: 1.38s
% Computational Cost: add. (899->161), mult. (1256->213), div. (0->0), fcn. (1423->14), ass. (0->74)
t144 = mrSges(6,2) - mrSges(7,3);
t75 = sin(qJ(6));
t79 = cos(qJ(6));
t143 = mrSges(7,1) * t79 - mrSges(7,2) * t75 + mrSges(6,1);
t145 = -m(7) * pkin(5) - t143;
t100 = -m(7) * pkin(11) + t144;
t74 = -qJ(4) - pkin(9);
t83 = -m(4) * pkin(9) + m(5) * t74 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t139 = t75 * mrSges(7,1) + t79 * mrSges(7,2) - t83;
t72 = qJ(3) + pkin(12);
t68 = qJ(5) + t72;
t63 = sin(t68);
t64 = cos(t68);
t80 = cos(qJ(3));
t69 = t80 * pkin(3);
t65 = t69 + pkin(2);
t66 = sin(t72);
t67 = cos(t72);
t76 = sin(qJ(3));
t141 = -m(4) * pkin(2) - m(5) * t65 - t80 * mrSges(4,1) - t67 * mrSges(5,1) + t76 * mrSges(4,2) + t66 * mrSges(5,2) + t100 * t63 + t145 * t64 - mrSges(3,1);
t140 = m(6) + m(7);
t138 = m(5) * pkin(3) + mrSges(4,1);
t112 = cos(pkin(6));
t73 = sin(pkin(6));
t77 = sin(qJ(2));
t122 = t73 * t77;
t34 = t112 * t64 - t122 * t63;
t35 = t112 * t63 + t122 * t64;
t137 = -t143 * t34 + t144 * t35;
t78 = sin(qJ(1));
t121 = t73 * t78;
t101 = t78 * t112;
t125 = cos(qJ(1));
t81 = cos(qJ(2));
t48 = -t101 * t77 + t125 * t81;
t17 = -t121 * t64 + t48 * t63;
t18 = t121 * t63 + t48 * t64;
t136 = t143 * t17 + t144 * t18;
t109 = t73 * t125;
t98 = t112 * t125;
t46 = t77 * t98 + t78 * t81;
t13 = -t64 * t109 - t46 * t63;
t14 = -t63 * t109 + t46 * t64;
t135 = -t143 * t13 + t144 * t14;
t128 = pkin(3) * t76;
t120 = t73 * t80;
t119 = t73 * t81;
t54 = pkin(4) * t66 + t128;
t55 = pkin(4) * t67 + t69;
t116 = t55 * t121 - t48 * t54;
t114 = t112 * t55 - t54 * t122;
t113 = t125 * pkin(1) + pkin(8) * t121;
t107 = -t78 * pkin(1) + pkin(8) * t109;
t106 = t13 * pkin(5) + t14 * pkin(11);
t105 = -t17 * pkin(5) + pkin(11) * t18;
t104 = t34 * pkin(5) + pkin(11) * t35;
t103 = t66 * t109 - t46 * t67;
t58 = t76 * t109;
t102 = -t46 * t80 + t58;
t47 = t101 * t81 + t125 * t77;
t53 = pkin(2) + t55;
t71 = -pkin(10) + t74;
t99 = t54 * t121 - t47 * t71 + t48 * t53 + t113;
t92 = -t109 * t55 - t46 * t54;
t21 = t120 * t78 - t48 * t76;
t88 = t109 * t67 + t46 * t66;
t87 = t109 * t80 + t46 * t76;
t45 = t77 * t78 - t81 * t98;
t22 = t121 * t76 + t48 * t80;
t20 = t121 * t66 + t48 * t67;
t19 = t121 * t67 - t48 * t66;
t2 = t18 * t79 + t47 * t75;
t1 = -t18 * t75 + t47 * t79;
t3 = [(-t125 * mrSges(2,1) - m(3) * t113 - t48 * mrSges(3,1) - m(4) * (pkin(2) * t48 + t113) - t22 * mrSges(4,1) - t21 * mrSges(4,2) - m(5) * (t48 * t65 + t113) - t20 * mrSges(5,1) - t19 * mrSges(5,2) - m(6) * t99 - t18 * mrSges(6,1) - m(7) * (pkin(5) * t18 + t99) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + (mrSges(2,2) + (-m(5) * t128 - mrSges(3,3)) * t73) * t78 + t100 * t17 + t83 * t47) * g(2) + (t78 * mrSges(2,1) + t125 * mrSges(2,2) - m(3) * t107 + t46 * mrSges(3,1) - mrSges(3,3) * t109 - m(4) * (-pkin(2) * t46 + t107) - t102 * mrSges(4,1) - t87 * mrSges(4,2) - m(5) * (pkin(3) * t58 - t46 * t65 + t107) - t103 * mrSges(5,1) - t88 * mrSges(5,2) + t100 * t13 - t145 * t14 + t139 * t45 + t140 * (-t54 * t109 - t45 * t71 + t46 * t53 - t107)) * g(1) (-t140 * (-t45 * t53 - t46 * t71) - t139 * t46 - t141 * t45) * g(2) + (-t140 * (-t47 * t53 - t48 * t71) - t139 * t48 - t141 * t47) * g(1) + (-t140 * t53 * t119 + (t141 * t81 + (t140 * t71 - t139) * t77) * t73) * g(3) (-(-t112 * t76 - t120 * t77) * mrSges(4,2) - (t112 * t67 - t122 * t66) * mrSges(5,1) - (-t112 * t66 - t122 * t67) * mrSges(5,2) - m(6) * t114 - m(7) * (t104 + t114) - t138 * (t112 * t80 - t122 * t76) + t137) * g(3) + (-t102 * mrSges(4,2) + t88 * mrSges(5,1) - t103 * mrSges(5,2) - m(6) * t92 - m(7) * (t106 + t92) + t138 * t87 + t135) * g(2) + (t22 * mrSges(4,2) - t19 * mrSges(5,1) + t20 * mrSges(5,2) - m(6) * t116 - m(7) * (t105 + t116) - t138 * t21 + t136) * g(1) (m(5) + t140) * (-g(1) * t47 - g(2) * t45 + g(3) * t119) (-m(7) * t104 + t137) * g(3) + (-m(7) * t106 + t135) * g(2) + (-m(7) * t105 + t136) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((-t14 * t75 + t45 * t79) * mrSges(7,1) + (-t14 * t79 - t45 * t75) * mrSges(7,2)) - g(3) * ((-t119 * t79 - t35 * t75) * mrSges(7,1) + (t119 * t75 - t35 * t79) * mrSges(7,2))];
taug  = t3(:);
