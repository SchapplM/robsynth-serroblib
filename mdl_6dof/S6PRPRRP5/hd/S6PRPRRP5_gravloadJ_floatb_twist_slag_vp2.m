% Calculate Gravitation load on the joints for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:14:09
% EndTime: 2019-03-08 20:14:11
% DurationCPUTime: 0.78s
% Computational Cost: add. (366->88), mult. (903->132), div. (0->0), fcn. (1029->10), ass. (0->52)
t89 = -mrSges(5,2) + mrSges(7,3);
t81 = -m(5) - m(6) - m(7);
t58 = m(4) - t81;
t88 = mrSges(6,1) + mrSges(7,1);
t82 = -mrSges(6,2) - mrSges(7,2);
t38 = sin(qJ(4));
t41 = cos(qJ(4));
t87 = m(6) * (pkin(4) * t38 - pkin(9) * t41) - t41 * mrSges(6,3);
t40 = cos(qJ(5));
t33 = pkin(5) * t40 + pkin(4);
t36 = -qJ(6) - pkin(9);
t86 = -m(7) * (t33 * t38 + t36 * t41) - t38 * mrSges(5,1) + mrSges(3,2) - mrSges(4,3) + t89 * t41;
t75 = m(7) * pkin(5);
t80 = mrSges(3,1) - mrSges(4,2) + mrSges(5,3);
t37 = sin(qJ(5));
t79 = m(6) * pkin(4) + m(7) * t33 + t82 * t37 + t88 * t40 + mrSges(5,1);
t78 = -t75 - t88;
t77 = -m(6) * pkin(9) + m(7) * t36 - mrSges(6,3) - t89;
t76 = -t58 * qJ(3) + t86 - t87;
t34 = sin(pkin(10));
t39 = sin(qJ(2));
t42 = cos(qJ(2));
t60 = cos(pkin(10));
t61 = cos(pkin(6));
t48 = t61 * t60;
t21 = t34 * t39 - t42 * t48;
t74 = t21 * t37;
t53 = t34 * t61;
t23 = t39 * t60 + t42 * t53;
t73 = t23 * t37;
t35 = sin(pkin(6));
t72 = t34 * t35;
t71 = t35 * t39;
t70 = t35 * t42;
t69 = t37 * t38;
t68 = t37 * t39;
t67 = t37 * t42;
t66 = t38 * t40;
t65 = t39 * t40;
t62 = pkin(2) * t70 + qJ(3) * t71;
t52 = t35 * t60;
t26 = -t38 * t70 + t41 * t61;
t25 = t38 * t61 + t41 * t70;
t24 = -t39 * t53 + t42 * t60;
t22 = t34 * t42 + t39 * t48;
t20 = t23 * pkin(2);
t19 = t21 * pkin(2);
t12 = -t21 * t38 + t41 * t52;
t11 = t21 * t41 + t38 * t52;
t10 = t23 * t38 + t41 * t72;
t9 = -t23 * t41 + t38 * t72;
t1 = [(-m(2) - m(3) - t58) * g(3) (t74 * t75 + m(4) * t19 - t88 * (t22 * t66 - t74) + t81 * (-pkin(8) * t21 - t19) + t82 * (-t21 * t40 - t22 * t69) + t80 * t21 + t76 * t22) * g(2) + (t73 * t75 + m(4) * t20 - t88 * (t24 * t66 - t73) + t82 * (-t23 * t40 - t24 * t69) + t81 * (-pkin(8) * t23 - t20) + t80 * t23 + t76 * t24) * g(1) + (-m(4) * t62 - t87 * t71 + t81 * (pkin(8) * t70 + t62) + (-t88 * (t38 * t65 + t67) + t82 * (-t38 * t68 + t40 * t42) - t67 * t75 - t80 * t42 + t86 * t39) * t35) * g(3) (-g(1) * t23 - g(2) * t21 + g(3) * t70) * t58 (t79 * t25 + t77 * t26) * g(3) + (-t79 * t11 - t77 * t12) * g(2) + (t77 * t10 + t79 * t9) * g(1) (t82 * (-t26 * t40 - t35 * t68) + t78 * (-t26 * t37 + t35 * t65)) * g(3) + (t82 * (t12 * t40 - t22 * t37) + t78 * (t12 * t37 + t22 * t40)) * g(2) + (t82 * (-t10 * t40 - t24 * t37) + t78 * (-t10 * t37 + t24 * t40)) * g(1) (-g(1) * t9 + g(2) * t11 - g(3) * t25) * m(7)];
taug  = t1(:);
