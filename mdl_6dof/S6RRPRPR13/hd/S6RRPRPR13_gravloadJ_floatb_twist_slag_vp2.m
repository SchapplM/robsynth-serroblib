% Calculate Gravitation load on the joints for
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR13_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR13_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR13_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR13_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:24:11
% EndTime: 2019-03-09 11:24:14
% DurationCPUTime: 1.03s
% Computational Cost: add. (561->98), mult. (1240->129), div. (0->0), fcn. (1429->12), ass. (0->55)
t36 = pkin(11) + qJ(6);
t33 = sin(t36);
t34 = cos(t36);
t37 = sin(pkin(11));
t39 = cos(pkin(11));
t49 = -m(7) * (pkin(5) * t39 + pkin(4)) - m(6) * pkin(4) - mrSges(6,1) * t39 + mrSges(6,2) * t37 - mrSges(5,1);
t114 = -mrSges(7,1) * t34 + mrSges(7,2) * t33 + t49;
t100 = m(6) + m(7);
t108 = m(5) + t100;
t73 = m(4) + t108;
t113 = qJ(3) * t73;
t112 = -t39 * mrSges(6,2) - (m(7) * pkin(5) + mrSges(6,1)) * t37 - mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t111 = pkin(9) * t108;
t53 = mrSges(5,2) - m(6) * qJ(5) + m(7) * (-pkin(10) - qJ(5)) - mrSges(6,3) - mrSges(7,3);
t61 = t33 * mrSges(7,1) + t34 * mrSges(7,2);
t105 = t112 - t61;
t41 = sin(qJ(4));
t44 = cos(qJ(4));
t80 = mrSges(4,3) - mrSges(3,2);
t104 = t114 * t41 - t44 * t53 - t80;
t103 = pkin(2) * t73 - t105 + t111;
t94 = t80 + t113;
t92 = -t111 + t112;
t90 = t104 - t113;
t38 = sin(pkin(6));
t42 = sin(qJ(2));
t85 = t38 * t42;
t43 = sin(qJ(1));
t84 = t38 * t43;
t45 = cos(qJ(2));
t83 = t38 * t45;
t46 = cos(qJ(1));
t82 = t38 * t46;
t79 = pkin(2) * t83 + qJ(3) * t85;
t78 = t46 * pkin(1) + pkin(8) * t84;
t77 = cos(pkin(6));
t68 = t43 * t77;
t21 = -t42 * t68 + t45 * t46;
t75 = t21 * pkin(2) + t78;
t70 = -pkin(1) * t43 + pkin(8) * t82;
t67 = t46 * t77;
t19 = t42 * t67 + t43 * t45;
t65 = -t19 * pkin(2) + t70;
t64 = mrSges(2,2) + (-mrSges(4,1) - mrSges(3,3)) * t38;
t18 = t42 * t43 - t45 * t67;
t7 = -t18 * t41 + t44 * t82;
t5 = t18 * t44 + t41 * t82;
t20 = t46 * t42 + t45 * t68;
t17 = -t41 * t83 + t44 * t77;
t16 = t41 * t77 + t44 * t83;
t4 = t20 * t41 + t44 * t84;
t3 = -t20 * t44 + t41 * t84;
t2 = t21 * t33 + t34 * t4;
t1 = t21 * t34 - t33 * t4;
t6 = [(-m(3) * t78 - m(4) * t75 - t46 * mrSges(2,1) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t94 * t20 + t92 * t21 + t53 * t3 + t49 * t4 + t64 * t43 - t108 * (pkin(3) * t84 + t75)) * g(2) + (t43 * mrSges(2,1) - m(3) * t70 - m(4) * t65 + t94 * t18 + t53 * t5 + t64 * t46 + t114 * t7 + (t61 - t92) * t19 + t108 * (-pkin(3) * t82 - t65)) * g(1) (-m(4) * t79 - t108 * (pkin(9) * t83 + t79) + (t104 * t42 + t105 * t45) * t38) * g(3) + (t103 * t18 + t90 * t19) * g(2) + (t103 * t20 + t90 * t21) * g(1) (-g(1) * t20 - g(2) * t18 + g(3) * t83) * t73 (-t114 * t16 + t53 * t17) * g(3) + (t114 * t5 - t53 * t7) * g(2) + (-t114 * t3 + t53 * t4) * g(1), t100 * (-g(1) * t3 + g(2) * t5 - g(3) * t16) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t19 * t34 + t33 * t7) * mrSges(7,1) + (-t19 * t33 + t34 * t7) * mrSges(7,2)) - g(3) * ((-t17 * t33 + t34 * t85) * mrSges(7,1) + (-t17 * t34 - t33 * t85) * mrSges(7,2))];
taug  = t6(:);
