% Calculate Gravitation load on the joints for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
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
% Datum: 2018-11-23 16:11
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR13_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR13_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR13_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:10:36
% EndTime: 2018-11-23 16:10:38
% DurationCPUTime: 1.13s
% Computational Cost: add. (2940->127), mult. (3127->171), div. (0->0), fcn. (3073->22), ass. (0->83)
t140 = m(6) + m(7);
t136 = m(5) + t140;
t96 = -t136 * qJ(4) + mrSges(4,2) - mrSges(5,3);
t66 = sin(qJ(6));
t70 = cos(qJ(6));
t97 = -pkin(10) * t140 - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t147 = -t66 * mrSges(7,1) - t70 * mrSges(7,2) + t97;
t139 = m(7) * pkin(5) + t70 * mrSges(7,1) - t66 * mrSges(7,2) + mrSges(6,1);
t146 = m(7) * pkin(11) - mrSges(6,2) + mrSges(7,3);
t62 = sin(pkin(6));
t69 = sin(qJ(1));
t126 = t62 * t69;
t61 = sin(pkin(7));
t64 = cos(pkin(7));
t122 = sin(pkin(12));
t73 = cos(qJ(1));
t116 = pkin(6) - pkin(12);
t101 = cos(t116) / 0.2e1;
t115 = pkin(6) + pkin(12);
t107 = cos(t115);
t83 = t101 + t107 / 0.2e1;
t79 = t122 * t73 + t69 * t83;
t33 = t64 * t126 + t79 * t61;
t44 = t122 * t69 - t73 * t83;
t31 = t62 * t64 * t73 - t44 * t61;
t145 = pkin(3) * t136 - t147;
t119 = pkin(7) + qJ(3);
t103 = sin(t119) / 0.2e1;
t120 = pkin(7) - qJ(3);
t110 = sin(t120);
t137 = t103 - t110 / 0.2e1;
t52 = t101 - t107 / 0.2e1;
t72 = cos(qJ(3));
t100 = sin(t115) / 0.2e1;
t106 = sin(t116);
t82 = t100 + t106 / 0.2e1;
t144 = t137 * t82 + t52 * t72;
t51 = t100 - t106 / 0.2e1;
t63 = cos(pkin(12));
t46 = -t69 * t51 + t63 * t73;
t143 = -t137 * t79 + t46 * t72;
t45 = t51 * t73 + t69 * t63;
t142 = -t137 * t44 + t45 * t72;
t68 = sin(qJ(3));
t85 = t103 + t110 / 0.2e1;
t80 = t62 * t85;
t111 = cos(t119);
t104 = t111 / 0.2e1;
t112 = cos(t120);
t105 = t112 / 0.2e1;
t88 = t105 + t104;
t17 = t44 * t88 + t45 * t68 + t73 * t80;
t67 = sin(qJ(5));
t71 = cos(qJ(5));
t141 = t17 * t71 + t31 * t67;
t6 = t17 * t67 - t31 * t71;
t132 = -t139 * t67 + t146 * t71 + t96;
t125 = -mrSges(4,3) - mrSges(5,1);
t123 = qJ(2) * t62;
t124 = t73 * pkin(1) + t69 * t123;
t113 = -pkin(1) * t69 + t73 * t123;
t99 = t105 - t111 / 0.2e1;
t92 = t62 * t99;
t90 = -t45 * pkin(2) + t31 * pkin(9) + t113;
t18 = -t73 * t92 + t142;
t89 = -pkin(3) * t18 + t90;
t87 = t104 - t112 / 0.2e1;
t81 = t62 * t87;
t76 = t46 * pkin(2) + t33 * pkin(9) + t124;
t23 = t69 * t92 + t143;
t75 = t23 * pkin(3) + t76;
t74 = t33 * pkin(4) + t75;
t65 = cos(pkin(6));
t43 = -t61 * t82 + t65 * t64;
t27 = t65 * t99 + t144;
t26 = t52 * t68 - t65 * t85 - t82 * t88;
t22 = t46 * t68 - t69 * t80 + t79 * t88;
t10 = t26 * t67 + t43 * t71;
t8 = t22 * t67 + t33 * t71;
t7 = -t22 * t71 + t33 * t67;
t2 = t23 * t66 + t70 * t8;
t1 = t23 * t70 - t66 * t8;
t3 = [(-t73 * mrSges(2,1) + t69 * mrSges(2,2) - m(3) * t124 - t46 * mrSges(3,1) + t79 * mrSges(3,2) - mrSges(3,3) * t126 - m(4) * t76 - m(5) * t75 - m(6) * t74 - t8 * mrSges(6,1) - m(7) * (t8 * pkin(5) + t74) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t146 * t7 + t125 * t33 + t96 * t22 + t97 * t23) * g(2) + (t69 * mrSges(2,1) - m(3) * t113 + t45 * mrSges(3,1) - t44 * mrSges(3,2) - m(4) * t90 - m(5) * t89 + (-t62 * mrSges(3,3) + mrSges(2,2)) * t73 + t125 * t31 - t146 * t141 - t96 * t17 + t139 * t6 - t147 * t18 + t140 * (-t31 * pkin(4) - t89)) * g(1) (-g(3) * t65 + (-g(1) * t69 + g(2) * t73) * t62) * (m(3) + m(4) + t136) (t132 * (-t65 * t87 + t144) + t145 * t26) * g(3) + (t132 * (t73 * t81 + t142) + t145 * t17) * g(2) + (t132 * (-t69 * t81 + t143) + t145 * t22) * g(1), t136 * (-g(1) * t22 - g(2) * t17 - g(3) * t26) (-t139 * (t26 * t71 - t43 * t67) - t146 * t10) * g(3) + (-t139 * t141 - t146 * t6) * g(2) + (t139 * t7 - t146 * t8) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t18 * t70 - t6 * t66) * mrSges(7,1) + (-t18 * t66 - t6 * t70) * mrSges(7,2)) - g(3) * ((-t10 * t66 + t27 * t70) * mrSges(7,1) + (-t10 * t70 - t27 * t66) * mrSges(7,2))];
taug  = t3(:);
