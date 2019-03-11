% Calculate Gravitation load on the joints for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR13_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR13_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR13_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:45:09
% EndTime: 2019-03-09 14:45:12
% DurationCPUTime: 1.09s
% Computational Cost: add. (622->117), mult. (1365->154), div. (0->0), fcn. (1588->12), ass. (0->62)
t61 = mrSges(5,2) - m(6) * pkin(10) + m(7) * (-pkin(11) - pkin(10)) - mrSges(6,3) - mrSges(7,3);
t115 = m(7) * pkin(5);
t84 = m(5) + m(6) + m(7);
t79 = m(4) + t84;
t114 = qJ(3) * t79;
t44 = qJ(5) + qJ(6);
t41 = sin(t44);
t42 = cos(t44);
t50 = cos(qJ(5));
t111 = -t41 * mrSges(7,1) - t50 * mrSges(6,2) - t42 * mrSges(7,2);
t46 = sin(qJ(5));
t82 = mrSges(4,2) - mrSges(3,1) - mrSges(5,3);
t83 = t46 * t115;
t113 = -t46 * mrSges(6,1) + t111 + t82 - t83;
t40 = pkin(5) * t50 + pkin(4);
t104 = m(6) * pkin(4) + m(7) * t40 + t50 * mrSges(6,1) + t42 * mrSges(7,1) - t46 * mrSges(6,2) - t41 * mrSges(7,2) + mrSges(5,1);
t112 = t84 * pkin(9);
t47 = sin(qJ(4));
t51 = cos(qJ(4));
t88 = mrSges(4,3) - mrSges(3,2);
t110 = -t104 * t47 - t61 * t51 - t88;
t109 = pkin(2) * t79 + t112 - t113;
t106 = mrSges(6,1) + t115;
t103 = -t82 + t112;
t102 = t88 + t114;
t99 = t110 - t114;
t48 = sin(qJ(2));
t52 = cos(qJ(2));
t53 = cos(qJ(1));
t49 = sin(qJ(1));
t85 = cos(pkin(6));
t75 = t49 * t85;
t28 = t53 * t48 + t52 * t75;
t45 = sin(pkin(6));
t92 = t45 * t49;
t12 = t28 * t47 + t51 * t92;
t29 = -t48 * t75 + t52 * t53;
t5 = -t12 * t41 + t29 * t42;
t6 = t12 * t42 + t29 * t41;
t97 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t74 = t53 * t85;
t26 = t48 * t49 - t52 * t74;
t90 = t45 * t53;
t15 = -t26 * t47 + t51 * t90;
t27 = t48 * t74 + t49 * t52;
t96 = (t15 * t41 + t27 * t42) * mrSges(7,1) + (t15 * t42 - t27 * t41) * mrSges(7,2);
t91 = t45 * t52;
t25 = -t47 * t91 + t51 * t85;
t93 = t45 * t48;
t94 = (-t25 * t41 + t42 * t93) * mrSges(7,1) + (-t25 * t42 - t41 * t93) * mrSges(7,2);
t87 = pkin(2) * t91 + qJ(3) * t93;
t86 = t53 * pkin(1) + pkin(8) * t92;
t81 = t29 * pkin(2) + t86;
t76 = -t49 * pkin(1) + pkin(8) * t90;
t73 = pkin(3) * t92 + t81;
t72 = -t27 * pkin(2) + t76;
t71 = mrSges(2,2) + (-mrSges(4,1) - mrSges(3,3)) * t45;
t7 = -t12 * t46 + t29 * t50;
t13 = t26 * t51 + t47 * t90;
t11 = -t28 * t51 + t47 * t92;
t8 = t12 * t50 + t29 * t46;
t1 = [(-t53 * mrSges(2,1) - m(3) * t86 - m(4) * t81 - m(5) * t73 - t12 * mrSges(5,1) - m(6) * (pkin(4) * t12 + t73) - t8 * mrSges(6,1) - t7 * mrSges(6,2) - m(7) * (t12 * t40 + t73) - t6 * mrSges(7,1) - t5 * mrSges(7,2) - t102 * t28 + t71 * t49 + t61 * t11 + (-t103 - t83) * t29) * g(2) + (t49 * mrSges(2,1) - m(3) * t76 - m(4) * t72 + t102 * t26 + t71 * t53 + t61 * t13 - t104 * t15 + (t106 * t46 + t103 - t111) * t27 + t84 * (-pkin(3) * t90 - t72)) * g(1) (-m(4) * t87 - t84 * (pkin(9) * t91 + t87) + (t110 * t48 + t113 * t52) * t45) * g(3) + (t109 * t26 + t99 * t27) * g(2) + (t109 * t28 + t99 * t29) * g(1) (-g(1) * t28 - g(2) * t26 + g(3) * t91) * t79 (t61 * t25 - t104 * (-t47 * t85 - t51 * t91)) * g(3) + (-t104 * t13 - t61 * t15) * g(2) + (t104 * t11 + t61 * t12) * g(1) (-(-t25 * t50 - t46 * t93) * mrSges(6,2) - t94 - t106 * (-t25 * t46 + t50 * t93)) * g(3) + (-(t15 * t50 - t27 * t46) * mrSges(6,2) - t96 - t106 * (t15 * t46 + t27 * t50)) * g(2) + (t8 * mrSges(6,2) - t106 * t7 - t97) * g(1), -g(1) * t97 - g(2) * t96 - g(3) * t94];
taug  = t1(:);
