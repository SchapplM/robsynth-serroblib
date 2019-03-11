% Calculate Gravitation load on the joints for
% S6RRPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:36:34
% EndTime: 2019-03-09 12:36:37
% DurationCPUTime: 1.12s
% Computational Cost: add. (841->137), mult. (1472->191), div. (0->0), fcn. (1730->12), ass. (0->63)
t93 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t91 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t125 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t68 = sin(pkin(11));
t70 = cos(pkin(11));
t130 = -m(4) * pkin(2) - t70 * mrSges(4,1) + t68 * mrSges(4,2) - mrSges(3,1);
t67 = pkin(11) + qJ(4);
t64 = sin(t67);
t65 = cos(t67);
t127 = t65 * mrSges(5,1) - mrSges(5,2) * t64 - t130;
t120 = cos(qJ(1));
t69 = sin(pkin(6));
t102 = t69 * t120;
t73 = sin(qJ(2));
t74 = sin(qJ(1));
t76 = cos(qJ(2));
t106 = cos(pkin(6));
t89 = t106 * t120;
t46 = t73 * t89 + t74 * t76;
t18 = -t64 * t102 + t46 * t65;
t45 = t73 * t74 - t76 * t89;
t72 = sin(qJ(5));
t75 = cos(qJ(5));
t1 = t18 * t72 - t45 * t75;
t132 = t18 * t75 + t45 * t72;
t131 = m(6) + m(7);
t128 = mrSges(6,3) + mrSges(7,2);
t126 = t91 * t72 - t93 * t75 - mrSges(5,1);
t124 = mrSges(5,2) - t128;
t122 = pkin(4) * t65;
t119 = t45 * t64;
t96 = t74 * t106;
t47 = t120 * t73 + t76 * t96;
t116 = t47 * t64;
t115 = t65 * t72;
t114 = t65 * t75;
t113 = t69 * t73;
t112 = t69 * t74;
t111 = t69 * t76;
t110 = t75 * t76;
t63 = pkin(3) * t70 + pkin(2);
t71 = -pkin(9) - qJ(3);
t109 = -t45 * t63 - t46 * t71;
t48 = t120 * t76 - t73 * t96;
t108 = -t47 * t63 - t48 * t71;
t107 = t120 * pkin(1) + pkin(8) * t112;
t105 = t64 * t111;
t104 = t72 * t111;
t101 = -pkin(1) * t74 + pkin(8) * t102;
t17 = -t65 * t102 - t46 * t64;
t92 = t68 * t102;
t90 = t68 * pkin(3) * t112 - t47 * t71 + t48 * t63 + t107;
t82 = pkin(3) * t92 + t45 * t71 - t46 * t63 + t101;
t80 = -t131 * pkin(10) + t124;
t49 = t63 * t111;
t35 = t106 * t64 + t113 * t65;
t34 = t106 * t65 - t113 * t64;
t22 = t112 * t64 + t48 * t65;
t21 = -t112 * t65 + t48 * t64;
t15 = t110 * t69 + t35 * t72;
t6 = t22 * t75 + t47 * t72;
t5 = t22 * t72 - t47 * t75;
t2 = [(-t120 * mrSges(2,1) - m(5) * t90 - t22 * mrSges(5,1) + t130 * t48 + (mrSges(2,2) + (-mrSges(4,1) * t68 - mrSges(4,2) * t70 - mrSges(3,3)) * t69) * t74 - t93 * t6 + t91 * t5 + t125 * t47 + t80 * t21 - t131 * (t22 * pkin(4) + t90) + (-m(3) - m(4)) * t107) * g(2) + (t74 * mrSges(2,1) + t120 * mrSges(2,2) - m(3) * t101 + t46 * mrSges(3,1) - mrSges(3,3) * t102 - m(4) * (-pkin(2) * t46 + t101) - (-t46 * t70 + t92) * mrSges(4,1) - (t102 * t70 + t46 * t68) * mrSges(4,2) - m(5) * t82 + t18 * mrSges(5,1) + t93 * t132 - t91 * t1 - t125 * t45 + t80 * t17 + t131 * (pkin(4) * t18 - t82)) * g(1) (-m(5) * t109 - t131 * (-pkin(10) * t119 - t45 * t122 + t109) - t93 * (-t114 * t45 + t46 * t72) + t91 * (-t115 * t45 - t46 * t75) + t128 * t119 + t125 * t46 + t127 * t45) * g(2) + (-m(5) * t108 - t131 * (-pkin(10) * t116 - t47 * t122 + t108) + t91 * (-t115 * t47 - t48 * t75) + t128 * t116 - t93 * (-t114 * t47 + t48 * t72) + t125 * t48 + t127 * t47) * g(1) + (-m(5) * t49 - t131 * (pkin(10) * t105 + t111 * t122 - t113 * t71 + t49) + t91 * (t104 * t65 - t113 * t75) - t128 * t105 + (-t93 * t110 * t65 - t127 * t76 + (m(5) * t71 - t93 * t72 + t125) * t73) * t69) * g(3) (-g(1) * t47 - g(2) * t45 + g(3) * t111) * (m(4) + m(5) + t131) (-t131 * (t34 * pkin(4) + pkin(10) * t35) + t124 * t35 + t126 * t34) * g(3) + (-t131 * (t17 * pkin(4) + pkin(10) * t18) + t124 * t18 + t126 * t17) * g(2) + (-t131 * (-t21 * pkin(4) + pkin(10) * t22) + t124 * t22 - t126 * t21) * g(1) (t91 * (t35 * t75 - t104) + t93 * t15) * g(3) + (t93 * t1 + t91 * t132) * g(2) + (t5 * t93 + t6 * t91) * g(1) (-g(1) * t5 - g(2) * t1 - g(3) * t15) * m(7)];
taug  = t2(:);
