% Calculate Gravitation load on the joints for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:46:47
% EndTime: 2019-03-09 02:46:48
% DurationCPUTime: 0.80s
% Computational Cost: add. (394->99), mult. (532->111), div. (0->0), fcn. (515->10), ass. (0->52)
t91 = mrSges(5,1) + mrSges(6,1);
t81 = mrSges(5,2) - mrSges(6,3);
t90 = -mrSges(6,2) - mrSges(5,3);
t30 = sin(qJ(1));
t32 = cos(qJ(1));
t80 = g(1) * t32 + g(2) * t30;
t89 = m(7) * pkin(8) + mrSges(7,3) + t90;
t24 = sin(pkin(10));
t26 = cos(pkin(10));
t29 = sin(qJ(6));
t31 = cos(qJ(6));
t43 = t24 * t29 + t26 * t31;
t44 = t24 * t31 - t26 * t29;
t88 = t43 * mrSges(7,1) + t44 * mrSges(7,2) - t81 * t24 + t91 * t26;
t27 = cos(pkin(9));
t85 = -mrSges(2,1) - m(3) * pkin(1) - t27 * mrSges(3,1) + sin(pkin(9)) * mrSges(3,2);
t23 = pkin(9) + qJ(3);
t21 = sin(t23);
t16 = t21 * qJ(4);
t22 = cos(t23);
t17 = t22 * pkin(3);
t59 = t17 + t16;
t84 = pkin(4) * t26 + qJ(5) * t24;
t82 = m(6) + m(7);
t56 = m(5) + t82;
t79 = -m(3) * qJ(2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t48 = t22 * mrSges(4,1) - t21 * mrSges(4,2);
t77 = t90 * t21 - t48;
t74 = m(7) * pkin(5) + t91;
t71 = pkin(5) * t26;
t68 = g(3) * t21;
t64 = t26 * t32;
t28 = -pkin(7) - qJ(2);
t63 = t28 * t32;
t62 = t30 * t24;
t61 = t30 * t26;
t60 = t32 * t24;
t20 = pkin(2) * t27 + pkin(1);
t55 = -t20 - t17;
t53 = t32 * t20 - t30 * t28;
t52 = t84 * t22 + t59;
t5 = t22 * t62 + t64;
t6 = t22 * t61 - t60;
t50 = t29 * t6 - t31 * t5;
t49 = -t29 * t5 - t31 * t6;
t42 = t59 * t32 + t53;
t41 = -pkin(3) - t84;
t8 = t22 * t64 + t62;
t7 = t22 * t60 - t61;
t2 = t29 * t7 + t31 * t8;
t1 = -t29 * t8 + t31 * t7;
t3 = [(-m(4) * t53 - m(5) * t42 - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t82 * (t8 * pkin(4) + t7 * qJ(5) + t42) - t74 * t8 + t81 * t7 + t79 * t30 + (t89 * t21 - t48 + t85) * t32) * g(2) + (m(5) * t63 - t49 * mrSges(7,1) - t50 * mrSges(7,2) - t82 * (-t6 * pkin(4) - t5 * qJ(5) - t63) + t74 * t6 - t81 * t5 + (m(4) * t28 + t79) * t32 + (m(4) * t20 - m(7) * t55 - (m(7) * (pkin(8) - qJ(4)) + mrSges(7,3)) * t21 + (-m(5) - m(6)) * (t55 - t16) - t77 - t85) * t30) * g(1) (-g(1) * t30 + g(2) * t32) * (m(3) + m(4) + t56) (-m(5) * t59 - m(6) * t52 - m(7) * (-pkin(8) * t21 + t52) + t21 * mrSges(7,3) + (-m(7) * t71 - t88) * t22 + t77) * g(3) + ((mrSges(4,1) - m(7) * (t41 - t71) - m(6) * t41 + m(5) * pkin(3) + t88) * t21 + (-qJ(4) * t56 + mrSges(4,2) + t89) * t22) * t80 (t22 * g(3) - t80 * t21) * t56, t82 * (-g(1) * t7 - g(2) * t5 - t24 * t68) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (-t50 * mrSges(7,1) + t49 * mrSges(7,2)) - (t44 * mrSges(7,1) - t43 * mrSges(7,2)) * t68];
taug  = t3(:);
