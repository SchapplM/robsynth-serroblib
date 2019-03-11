% Calculate Gravitation load on the joints for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPPRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:17:43
% EndTime: 2019-03-08 19:17:45
% DurationCPUTime: 0.80s
% Computational Cost: add. (429->77), mult. (1093->122), div. (0->0), fcn. (1322->12), ass. (0->44)
t84 = m(6) + m(7);
t76 = m(5) + t84;
t32 = sin(qJ(6));
t35 = cos(qJ(6));
t78 = -m(7) * pkin(5) - t35 * mrSges(7,1) + t32 * mrSges(7,2) - mrSges(6,1);
t83 = m(7) * pkin(9) - mrSges(6,2) + mrSges(7,3);
t29 = sin(pkin(10));
t31 = cos(pkin(10));
t34 = sin(qJ(2));
t37 = cos(qJ(2));
t62 = cos(pkin(6));
t54 = t37 * t62;
t42 = -t29 * t54 - t31 * t34;
t41 = t42 * pkin(2);
t33 = sin(qJ(5));
t36 = cos(qJ(5));
t81 = -t76 * qJ(4) + t78 * t33 + t83 * t36 + mrSges(4,2) - mrSges(5,3);
t80 = -t29 * t34 + t31 * t54;
t60 = sin(pkin(11));
t61 = cos(pkin(11));
t21 = t34 * t60 - t37 * t61;
t59 = m(4) + t76;
t47 = t62 * t60;
t48 = t62 * t61;
t63 = -t34 * t48 - t37 * t47;
t13 = -t31 * t21 + t29 * t63;
t75 = t29 * t21 + t31 * t63;
t73 = t32 * mrSges(7,1) + t35 * mrSges(7,2) + t84 * pkin(8) + mrSges(4,1) - mrSges(5,2) + mrSges(6,3);
t30 = sin(pkin(6));
t69 = t30 * t33;
t68 = t30 * t36;
t27 = t30 * t37 * pkin(2);
t56 = t34 * t62;
t52 = t80 * pkin(2);
t22 = -t34 * t61 - t37 * t60;
t39 = -t34 * t47 + t37 * t48;
t20 = t22 * t30;
t19 = t21 * t30;
t15 = t19 * t33 + t62 * t36;
t12 = t22 * t31 - t29 * t39;
t9 = t29 * t22 + t31 * t39;
t4 = t31 * t68 + t9 * t33;
t2 = -t12 * t33 + t29 * t68;
t1 = [(-m(2) - m(3) - t59) * g(3) (-(mrSges(3,1) * t37 - mrSges(3,2) * t34) * t30 - m(4) * t27 - t76 * (-t19 * pkin(3) + t27) - t81 * t20 + t73 * t19) * g(3) + (-t80 * mrSges(3,1) - (-t29 * t37 - t31 * t56) * mrSges(3,2) - m(4) * t52 - t76 * (t9 * pkin(3) + t52) - t73 * t9 - t81 * t75) * g(2) + (-t42 * mrSges(3,1) - (t29 * t56 - t31 * t37) * mrSges(3,2) - m(4) * t41 - t73 * t12 + t81 * t13 - t76 * (t12 * pkin(3) + t41)) * g(1) ((-g(1) * t29 + g(2) * t31) * t30 - g(3) * t62) * t59, t76 * (g(1) * t12 + g(2) * t9 - g(3) * t19) (-t83 * t15 + t78 * (t19 * t36 - t62 * t33)) * g(3) + (t83 * t4 + t78 * (t31 * t69 - t36 * t9)) * g(2) + (-t83 * t2 + t78 * (-t12 * t36 - t29 * t69)) * g(1), -g(1) * ((t13 * t35 - t2 * t32) * mrSges(7,1) + (-t13 * t32 - t2 * t35) * mrSges(7,2)) - g(2) * ((t32 * t4 - t35 * t75) * mrSges(7,1) + (t32 * t75 + t35 * t4) * mrSges(7,2)) - g(3) * ((-t15 * t32 - t20 * t35) * mrSges(7,1) + (-t15 * t35 + t20 * t32) * mrSges(7,2))];
taug  = t1(:);
