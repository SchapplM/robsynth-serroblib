% Calculate Gravitation load on the joints for
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:30:33
% EndTime: 2019-03-09 02:30:33
% DurationCPUTime: 0.48s
% Computational Cost: add. (209->77), mult. (350->91), div. (0->0), fcn. (303->8), ass. (0->44)
t24 = cos(qJ(5));
t72 = m(7) * (pkin(5) * t24 + pkin(4)) + m(6) * pkin(4);
t20 = qJ(5) + qJ(6);
t14 = sin(t20);
t15 = cos(t20);
t21 = sin(qJ(5));
t69 = mrSges(6,1) * t24 + mrSges(7,1) * t15 - mrSges(6,2) * t21 - mrSges(7,2) * t14 + t72;
t68 = m(7) * (-pkin(9) - pkin(8)) - mrSges(7,3) - m(6) * pkin(8) - mrSges(6,3);
t23 = sin(qJ(1));
t26 = cos(qJ(1));
t67 = g(1) * t26 + g(2) * t23;
t61 = m(7) * pkin(5);
t66 = mrSges(6,1) + t61;
t65 = -m(5) - m(6) - m(7);
t42 = -m(4) + t65;
t63 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2) + mrSges(5,3);
t22 = sin(qJ(4));
t25 = cos(qJ(4));
t35 = t22 * mrSges(5,1) + t25 * mrSges(5,2);
t62 = t22 * t72 + t25 * t68 + mrSges(2,1) - mrSges(3,2) + mrSges(4,3) + t35;
t48 = t26 * t15;
t55 = t23 * t14;
t5 = t22 * t55 - t48;
t49 = t26 * t14;
t54 = t23 * t15;
t6 = -t22 * t54 - t49;
t60 = -mrSges(7,1) * t5 + mrSges(7,2) * t6;
t7 = -t22 * t49 - t54;
t8 = t22 * t48 - t55;
t59 = mrSges(7,1) * t7 - mrSges(7,2) * t8;
t56 = g(3) * t25;
t53 = t23 * t21;
t52 = t23 * t24;
t47 = t26 * t21;
t46 = t26 * t24;
t44 = pkin(1) * t26 + qJ(2) * t23;
t43 = qJ(3) * t26 + t44;
t34 = -mrSges(7,1) * t14 - mrSges(7,2) * t15;
t11 = -t22 * t47 - t52;
t9 = t22 * t53 - t46;
t18 = t26 * qJ(2);
t12 = t22 * t46 - t53;
t10 = -t22 * t52 - t47;
t1 = [(t53 * t61 - m(3) * t44 - m(4) * t43 - t12 * mrSges(6,1) - t8 * mrSges(7,1) - t11 * mrSges(6,2) - t7 * mrSges(7,2) + t65 * (-t23 * pkin(7) + t43) + t63 * t23 - t62 * t26) * g(2) + (t47 * t61 - t10 * mrSges(6,1) - t6 * mrSges(7,1) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + t65 * (-t26 * pkin(7) + t18) + (-m(3) - m(4)) * t18 + t63 * t26 + (m(3) * pkin(1) + t42 * (-pkin(1) - qJ(3)) + t62) * t23) * g(1) (-g(1) * t23 + g(2) * t26) * (m(3) - t42) t67 * t42, g(3) * t35 + (t68 * g(3) + t67 * (-mrSges(5,1) - t69)) * t25 + (t69 * g(3) + t67 * (mrSges(5,2) + t68)) * t22 (mrSges(6,2) * t24 + t21 * t66 - t34) * t56 + (-t10 * mrSges(6,2) + t66 * t9 - t60) * g(2) + (t12 * mrSges(6,2) - t11 * t66 - t59) * g(1), -g(1) * t59 - g(2) * t60 - t34 * t56];
taug  = t1(:);
