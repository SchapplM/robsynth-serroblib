% Calculate Gravitation load on the joints for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:22:55
% EndTime: 2019-03-09 02:22:57
% DurationCPUTime: 0.56s
% Computational Cost: add. (328->79), mult. (333->94), div. (0->0), fcn. (291->10), ass. (0->45)
t72 = mrSges(5,2) - m(6) * pkin(8) + m(7) * (-pkin(9) - pkin(8)) - mrSges(6,3) - mrSges(7,3);
t26 = sin(qJ(4));
t29 = cos(qJ(4));
t71 = t26 * mrSges(5,1) + t72 * t29;
t28 = cos(qJ(5));
t70 = -m(6) * pkin(4) - m(7) * (t28 * pkin(5) + pkin(4));
t67 = m(7) * pkin(5);
t23 = qJ(1) + pkin(10);
t18 = sin(t23);
t19 = cos(t23);
t68 = g(1) * t18 - g(2) * t19;
t25 = sin(qJ(5));
t66 = t25 * t67;
t65 = -m(5) - m(6) - m(7);
t64 = -mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t43 = m(4) - t65;
t63 = t70 * t26 + mrSges(3,2) - mrSges(4,3) - t71;
t61 = -mrSges(6,1) - t67;
t24 = qJ(5) + qJ(6);
t20 = sin(t24);
t21 = cos(t24);
t60 = t28 * mrSges(6,1) + t21 * mrSges(7,1) - t25 * mrSges(6,2) - t20 * mrSges(7,2) - t70;
t47 = t26 * t20;
t5 = -t18 * t47 + t19 * t21;
t46 = t26 * t21;
t6 = t18 * t46 + t19 * t20;
t56 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t7 = t18 * t21 + t19 * t47;
t8 = -t18 * t20 + t19 * t46;
t55 = t7 * mrSges(7,1) + t8 * mrSges(7,2);
t51 = g(3) * t29;
t27 = sin(qJ(1));
t50 = t27 * pkin(1);
t30 = cos(qJ(1));
t22 = t30 * pkin(1);
t49 = t19 * t25;
t48 = t25 * t26;
t45 = t26 * t28;
t44 = t19 * pkin(2) + t18 * qJ(3) + t22;
t36 = -mrSges(7,1) * t20 - mrSges(7,2) * t21;
t11 = t18 * t28 + t19 * t48;
t9 = -t18 * t48 + t19 * t28;
t12 = -t18 * t25 + t19 * t45;
t10 = t18 * t45 + t49;
t1 = [(-t49 * t67 - m(3) * t22 - m(4) * t44 - t30 * mrSges(2,1) - t10 * mrSges(6,1) - t6 * mrSges(7,1) + t27 * mrSges(2,2) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + t65 * (t19 * pkin(7) + t44) + t64 * t19 + t63 * t18) * g(2) + (m(3) * t50 + t27 * mrSges(2,1) - t12 * mrSges(6,1) - t8 * mrSges(7,1) + t30 * mrSges(2,2) + t11 * mrSges(6,2) + t7 * mrSges(7,2) - t43 * (t19 * qJ(3) - t50) + (m(4) * pkin(2) + t66 + t65 * (-pkin(2) - pkin(7)) - t64) * t18 + t63 * t19) * g(1) (-m(3) - t43) * g(3), -t68 * t43 (t60 * t26 + t71) * g(3) + t68 * ((-mrSges(5,1) - t60) * t29 + t72 * t26) (mrSges(6,1) * t25 + mrSges(6,2) * t28 - t36 + t66) * t51 + (-t12 * mrSges(6,2) + t61 * t11 - t55) * g(2) + (t10 * mrSges(6,2) + t61 * t9 - t56) * g(1), -g(1) * t56 - g(2) * t55 - t36 * t51];
taug  = t1(:);
