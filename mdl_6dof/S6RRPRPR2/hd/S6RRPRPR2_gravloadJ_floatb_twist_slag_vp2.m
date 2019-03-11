% Calculate Gravitation load on the joints for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:12:05
% EndTime: 2019-03-09 10:12:07
% DurationCPUTime: 0.64s
% Computational Cost: add. (462->113), mult. (420->120), div. (0->0), fcn. (347->10), ass. (0->61)
t41 = sin(qJ(6));
t44 = cos(qJ(6));
t110 = mrSges(7,1) * t41 + mrSges(7,2) * t44;
t39 = qJ(2) + pkin(10);
t34 = sin(t39);
t35 = cos(t39);
t42 = sin(qJ(2));
t45 = cos(qJ(2));
t109 = -t45 * mrSges(3,1) - t35 * mrSges(4,1) + t42 * mrSges(3,2) + t34 * mrSges(4,2);
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t101 = g(1) * t46 + g(2) * t43;
t36 = qJ(4) + t39;
t31 = sin(t36);
t32 = cos(t36);
t108 = (-mrSges(5,1) + mrSges(6,2)) * t32 + (mrSges(5,2) - mrSges(6,3)) * t31;
t61 = m(7) * (-pkin(4) - pkin(9)) - mrSges(7,3);
t107 = t110 * t32 + t31 * t61;
t96 = m(6) + m(7);
t24 = t31 * qJ(5);
t84 = t32 * t46;
t106 = pkin(4) * t84 + t46 * t24;
t86 = t31 * t46;
t103 = -mrSges(6,2) * t86 - mrSges(6,3) * t84 - t107 * t46;
t87 = t31 * t43;
t102 = -mrSges(6,2) * t87 + (-mrSges(6,3) * t32 - t107) * t43;
t99 = -t32 * mrSges(7,3) - t110 * t31 + t108;
t37 = t45 * pkin(2);
t98 = -mrSges(2,1) - m(3) * pkin(1) - m(4) * (t37 + pkin(1)) + t108 + t109;
t40 = -qJ(3) - pkin(7);
t38 = -pkin(8) + t40;
t97 = -m(3) * pkin(7) + m(4) * t40 - m(7) * (pkin(5) - t38) - mrSges(6,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t93 = g(3) * t32;
t29 = t32 * pkin(4);
t92 = t42 * pkin(2);
t83 = t43 * t41;
t82 = t43 * t44;
t81 = t46 * t41;
t80 = t46 * t44;
t77 = t29 + t24;
t76 = pkin(3) * t35 + t37;
t75 = qJ(5) * t32;
t12 = pkin(1) + t76;
t5 = t46 * t12;
t68 = -t43 * t38 + t5;
t65 = -t12 - t24;
t64 = t76 + t77;
t14 = t43 * t75;
t63 = -pkin(4) * t87 + t14;
t16 = t46 * t75;
t62 = -pkin(4) * t86 + t16;
t55 = mrSges(5,1) * t31 + mrSges(5,2) * t32;
t28 = t32 * pkin(9);
t13 = -pkin(3) * t34 - t92;
t7 = t46 * t13;
t6 = t43 * t13;
t4 = -t31 * t83 + t80;
t3 = t31 * t82 + t81;
t2 = t31 * t81 + t82;
t1 = t31 * t80 - t83;
t8 = [(-m(5) * t68 - m(6) * (t68 + t106) - m(7) * (pkin(9) * t84 + t106 + t5) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - mrSges(7,3) * t84 + t98 * t46 + t97 * t43) * g(2) + (-t4 * mrSges(7,1) + t3 * mrSges(7,2) + ((m(5) + m(6)) * t38 + t97) * t46 + (m(5) * t12 - m(6) * (t65 - t29) - m(7) * t65 - t61 * t32 - t98) * t43) * g(1) (-m(6) * (t6 + t63) - m(7) * (t14 + t6) + t102) * g(2) + (-m(6) * (t62 + t7) - m(7) * (t16 + t7) + t103) * g(1) + (-m(4) * t37 - m(5) * t76 - m(6) * t64 - m(7) * (t28 + t64) + t99 + t109) * g(3) + t101 * (m(4) * t92 - m(5) * t13 + mrSges(3,1) * t42 + mrSges(4,1) * t34 + mrSges(3,2) * t45 + mrSges(4,2) * t35 + t55) (-g(1) * t43 + g(2) * t46) * (m(4) + m(5) + t96) t101 * t55 + (-m(6) * t63 - m(7) * t14 + t102) * g(2) + (-m(6) * t62 - m(7) * t16 + t103) * g(1) + (-m(6) * t77 - m(7) * (t28 + t77) + t99) * g(3) (-t101 * t31 + t93) * t96, -g(1) * (t1 * mrSges(7,1) - t2 * mrSges(7,2)) - g(2) * (t3 * mrSges(7,1) + t4 * mrSges(7,2)) - (-mrSges(7,1) * t44 + mrSges(7,2) * t41) * t93];
taug  = t8(:);
