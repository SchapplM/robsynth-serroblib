% Calculate Gravitation load on the joints for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:25:39
% EndTime: 2019-03-09 02:25:40
% DurationCPUTime: 0.55s
% Computational Cost: add. (315->82), mult. (552->100), div. (0->0), fcn. (613->10), ass. (0->46)
t28 = cos(qJ(5));
t84 = m(7) * (t28 * pkin(5) + pkin(4)) + m(6) * pkin(4);
t64 = m(7) * pkin(5);
t74 = -mrSges(6,1) - t64;
t25 = qJ(5) + qJ(6);
t19 = sin(t25);
t20 = cos(t25);
t26 = sin(qJ(5));
t83 = t28 * mrSges(6,1) + t20 * mrSges(7,1) - t26 * mrSges(6,2) - t19 * mrSges(7,2) + t84;
t82 = m(7) * (-pkin(9) - pkin(8)) - mrSges(7,3) - m(6) * pkin(8) - mrSges(6,3);
t27 = sin(qJ(4));
t81 = t82 * t27;
t80 = t19 * mrSges(7,1) + t20 * mrSges(7,2);
t29 = cos(qJ(4));
t79 = t83 * t29 - t81;
t76 = mrSges(2,1) + mrSges(3,1);
t41 = t29 * mrSges(5,1) - t27 * mrSges(5,2);
t75 = -mrSges(4,1) - t41;
t73 = mrSges(2,2) - mrSges(3,3);
t72 = mrSges(4,2) - mrSges(5,3);
t69 = -m(5) - m(6) - m(7);
t67 = m(4) - t69;
t66 = -t28 * mrSges(6,2) + t74 * t26;
t46 = sin(pkin(10));
t47 = cos(pkin(10));
t57 = sin(qJ(1));
t58 = cos(qJ(1));
t11 = -t57 * t46 - t58 * t47;
t12 = t58 * t46 - t57 * t47;
t52 = t20 * t29;
t54 = t19 * t29;
t63 = (-t11 * t20 + t12 * t54) * mrSges(7,1) + (t11 * t19 + t12 * t52) * mrSges(7,2);
t5 = t11 * t54 + t12 * t20;
t6 = -t11 * t52 + t12 * t19;
t62 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t56 = t12 * t26;
t51 = t26 * t29;
t50 = t28 * t29;
t49 = t80 * t27;
t48 = t58 * pkin(1) + t57 * qJ(2);
t45 = t58 * pkin(2) + t48;
t42 = -t57 * pkin(1) + t58 * qJ(2);
t7 = t11 * t51 + t12 * t28;
t34 = -t57 * pkin(2) + t42;
t8 = -t11 * t50 + t56;
t1 = [(-t56 * t64 - m(3) * t48 - m(4) * t45 - t8 * mrSges(6,1) - t6 * mrSges(7,1) - t7 * mrSges(6,2) - t5 * mrSges(7,2) - t76 * t58 + t73 * t57 + t69 * (-t11 * pkin(3) + t12 * pkin(7) + t45) + t72 * t12 + (t84 * t29 - t75 - t81) * t11) * g(2) + (-m(3) * t42 - m(4) * t34 + t73 * t58 + t76 * t57 + t69 * (t12 * pkin(3) + t34) + (t75 - t79) * t12 + (t69 * pkin(7) + t66 + t72 - t80) * t11) * g(1) (-t57 * g(1) + t58 * g(2)) * (m(3) + t67) t67 * g(3) (t41 + t79) * g(3) + (g(1) * t11 + g(2) * t12) * ((-mrSges(5,2) - t82) * t29 + (-mrSges(5,1) - t83) * t27) (t66 * t27 - t49) * g(3) + (-(t11 * t26 + t12 * t50) * mrSges(6,2) - t63 + t74 * (-t11 * t28 + t12 * t51)) * g(2) + (t8 * mrSges(6,2) + t74 * t7 - t62) * g(1), -g(1) * t62 - g(2) * t63 - g(3) * t49];
taug  = t1(:);
