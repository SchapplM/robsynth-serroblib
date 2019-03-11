% Calculate Gravitation load on the joints for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:18:02
% EndTime: 2019-03-09 02:18:03
% DurationCPUTime: 0.44s
% Computational Cost: add. (405->87), mult. (294->96), div. (0->0), fcn. (240->12), ass. (0->50)
t30 = pkin(11) + qJ(4);
t26 = qJ(5) + t30;
t19 = sin(t26);
t20 = cos(t26);
t35 = sin(qJ(6));
t58 = t35 * mrSges(7,2);
t83 = -t19 * t58 + t20 * (-m(7) * pkin(9) - mrSges(7,3));
t76 = t20 * pkin(5) + t19 * pkin(9);
t82 = m(7) * t76;
t81 = -t20 * mrSges(6,1) + (mrSges(6,2) - mrSges(7,3)) * t19;
t78 = m(4) + m(5);
t77 = -m(6) - m(7);
t37 = cos(qJ(6));
t57 = t37 * mrSges(7,1);
t75 = -(t57 - t58) * t20 + t81;
t73 = m(3) + t78;
t33 = cos(pkin(11));
t21 = t33 * pkin(3) + pkin(2);
t22 = sin(t30);
t24 = cos(t30);
t47 = t24 * mrSges(5,1) - t22 * mrSges(5,2);
t72 = m(5) * t21 + mrSges(3,1) + m(4) * pkin(2) + t33 * mrSges(4,1) - sin(pkin(11)) * mrSges(4,2) + t47 - t81;
t34 = -pkin(7) - qJ(3);
t71 = -m(4) * qJ(3) + m(5) * t34 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t70 = pkin(4) * t22;
t18 = pkin(4) * t24;
t36 = sin(qJ(1));
t67 = t36 * pkin(1);
t38 = cos(qJ(1));
t28 = t38 * pkin(1);
t66 = mrSges(6,2) * t20;
t31 = qJ(1) + pkin(10);
t23 = sin(t31);
t62 = t23 * t35;
t61 = t23 * t37;
t25 = cos(t31);
t60 = t25 * t35;
t59 = t25 * t37;
t54 = -t77 + t78;
t52 = t83 * t23;
t51 = t83 * t25;
t40 = m(7) * (-pkin(5) * t19 - t70) - t19 * t57;
t39 = t66 + (m(7) * pkin(5) + mrSges(6,1) + t57) * t19;
t29 = -pkin(8) + t34;
t10 = t18 + t21;
t4 = t20 * t59 + t62;
t3 = -t20 * t60 + t61;
t2 = -t20 * t61 + t60;
t1 = t20 * t62 + t59;
t5 = [(-t38 * mrSges(2,1) - t4 * mrSges(7,1) + t36 * mrSges(2,2) - t3 * mrSges(7,2) + t77 * (t25 * t10 - t23 * t29 + t28) - t73 * t28 + t71 * t23 + (-t72 - t82) * t25) * g(2) + (t36 * mrSges(2,1) - t2 * mrSges(7,1) + t38 * mrSges(2,2) - t1 * mrSges(7,2) + t77 * (-t25 * t29 - t67) + t73 * t67 + t71 * t25 + (m(6) * t10 - m(7) * (-t10 - t76) + t72) * t23) * g(1) (-m(3) - t54) * g(3) (-g(1) * t23 + g(2) * t25) * t54, -g(1) * (t25 * t40 - t51) - g(2) * (t23 * t40 - t52) + (-t47 - m(6) * t18 - m(7) * (t18 + t76) + t75) * g(3) + (m(6) * t70 + mrSges(5,1) * t22 + mrSges(6,1) * t19 + mrSges(5,2) * t24 + t66) * (g(1) * t25 + g(2) * t23) (t75 - t82) * g(3) + (t23 * t39 + t52) * g(2) + (t25 * t39 + t51) * g(1), -g(1) * (t3 * mrSges(7,1) - t4 * mrSges(7,2)) - g(2) * (-t1 * mrSges(7,1) + t2 * mrSges(7,2)) - g(3) * (-mrSges(7,1) * t35 - mrSges(7,2) * t37) * t19];
taug  = t5(:);
