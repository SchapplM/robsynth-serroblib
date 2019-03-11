% Calculate Gravitation load on the joints for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:36:02
% EndTime: 2019-03-08 20:36:04
% DurationCPUTime: 0.83s
% Computational Cost: add. (584->86), mult. (942->122), div. (0->0), fcn. (1079->14), ass. (0->41)
t32 = qJ(5) + qJ(6);
t29 = sin(t32);
t30 = cos(t32);
t38 = sin(qJ(5));
t40 = cos(qJ(5));
t91 = mrSges(5,1) + m(7) * (pkin(5) * t40 + pkin(4)) + t30 * mrSges(7,1) - t29 * mrSges(7,2) + m(6) * pkin(4) + t40 * mrSges(6,1) - t38 * mrSges(6,2);
t77 = mrSges(5,2) + m(7) * (-pkin(10) - pkin(9)) - mrSges(7,3) - m(6) * pkin(9) - mrSges(6,3);
t31 = pkin(12) + qJ(4);
t27 = sin(t31);
t28 = cos(t31);
t36 = cos(pkin(12));
t76 = m(4) * pkin(2) + t36 * mrSges(4,1) - sin(pkin(12)) * mrSges(4,2) - t77 * t27 + mrSges(3,1) + t91 * t28;
t84 = -m(7) * pkin(5) - mrSges(6,1);
t75 = -m(4) * qJ(3) - t29 * mrSges(7,1) - t40 * mrSges(6,2) - t30 * mrSges(7,2) + t84 * t38 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t79 = m(5) + m(6) + m(7);
t34 = sin(pkin(11));
t39 = sin(qJ(2));
t41 = cos(qJ(2));
t62 = cos(pkin(11));
t63 = cos(pkin(6));
t50 = t63 * t62;
t17 = t34 * t39 - t41 * t50;
t18 = t34 * t41 + t39 * t50;
t35 = sin(pkin(6));
t56 = t35 * t62;
t8 = t18 * t28 - t27 * t56;
t73 = (t17 * t30 - t29 * t8) * mrSges(7,1) + (-t17 * t29 - t30 * t8) * mrSges(7,2);
t57 = t34 * t63;
t20 = -t39 * t57 + t41 * t62;
t68 = t34 * t35;
t10 = t20 * t28 + t27 * t68;
t19 = t39 * t62 + t41 * t57;
t72 = (-t10 * t29 + t19 * t30) * mrSges(7,1) + (-t10 * t30 - t19 * t29) * mrSges(7,2);
t67 = t35 * t39;
t14 = t27 * t63 + t28 * t67;
t66 = t35 * t41;
t71 = (-t14 * t29 - t30 * t66) * mrSges(7,1) + (-t14 * t30 + t29 * t66) * mrSges(7,2);
t61 = m(4) + t79;
t37 = -pkin(8) - qJ(3);
t25 = pkin(3) * t36 + pkin(2);
t1 = [(-m(2) - m(3) - t61) * g(3) (-t79 * (-t17 * t25 - t18 * t37) + t75 * t18 + t76 * t17) * g(2) + (-t79 * (-t19 * t25 - t20 * t37) + t75 * t20 + t76 * t19) * g(1) + (-t79 * t25 * t66 + (-t76 * t41 + (t79 * t37 + t75) * t39) * t35) * g(3) (-g(1) * t19 - g(2) * t17 + g(3) * t66) * t61 (t77 * t14 - t91 * (-t27 * t67 + t28 * t63)) * g(3) + (t77 * t8 - t91 * (-t18 * t27 - t28 * t56)) * g(2) + (-t91 * (-t20 * t27 + t28 * t68) + t77 * t10) * g(1) (-(-t14 * t40 + t38 * t66) * mrSges(6,2) - t71 + t84 * (-t14 * t38 - t40 * t66)) * g(3) + (-(-t17 * t38 - t40 * t8) * mrSges(6,2) - t73 + t84 * (t17 * t40 - t38 * t8)) * g(2) + (-(-t10 * t40 - t19 * t38) * mrSges(6,2) - t72 + t84 * (-t10 * t38 + t19 * t40)) * g(1), -g(1) * t72 - g(2) * t73 - g(3) * t71];
taug  = t1(:);
