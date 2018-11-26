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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 14:54
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:54:02
% EndTime: 2018-11-23 14:54:03
% DurationCPUTime: 0.66s
% Computational Cost: add. (915->91), mult. (763->125), div. (0->0), fcn. (674->20), ass. (0->60)
t91 = m(6) + m(7);
t46 = sin(qJ(6));
t49 = cos(qJ(6));
t88 = -m(7) * pkin(5) - mrSges(7,1) * t49 + mrSges(7,2) * t46 - mrSges(6,1);
t86 = m(5) + t91;
t90 = m(7) * pkin(9) - mrSges(6,2) + mrSges(7,3);
t41 = pkin(6) - qJ(2);
t31 = sin(t41) / 0.2e1;
t40 = pkin(6) + qJ(2);
t34 = sin(t40);
t89 = t34 / 0.2e1 + t31;
t32 = cos(t40) / 0.2e1;
t38 = cos(t41);
t24 = t38 / 0.2e1 + t32;
t85 = t46 * mrSges(7,1) + t49 * mrSges(7,2) + t91 * pkin(8) + mrSges(4,1) - mrSges(5,2) + mrSges(6,3);
t47 = sin(qJ(5));
t50 = cos(qJ(5));
t84 = -t86 * qJ(4) + t88 * t47 + t90 * t50 + mrSges(4,2) - mrSges(5,3);
t39 = qJ(2) + pkin(11);
t36 = cos(t39);
t42 = sin(pkin(10));
t81 = t42 * t36;
t48 = sin(qJ(2));
t80 = t42 * t48;
t43 = sin(pkin(6));
t79 = t43 * t47;
t78 = t43 * t50;
t44 = cos(pkin(10));
t77 = t44 * t36;
t76 = t44 * t48;
t75 = t89 * pkin(2);
t73 = m(4) + t86;
t71 = pkin(6) - t39;
t70 = pkin(6) + t39;
t22 = t24 * pkin(2);
t69 = -pkin(2) * t80 + t44 * t22;
t67 = cos(t70);
t66 = sin(t71);
t65 = sin(t70);
t64 = -pkin(2) * t76 - t42 * t22;
t63 = cos(t71) / 0.2e1;
t62 = t66 / 0.2e1;
t61 = t65 / 0.2e1;
t56 = -t65 / 0.2e1 + t62;
t55 = t67 / 0.2e1 + t63;
t54 = t61 - t66 / 0.2e1;
t51 = cos(qJ(2));
t45 = cos(pkin(6));
t33 = sin(t39);
t23 = t31 - t34 / 0.2e1;
t21 = t63 - t67 / 0.2e1;
t20 = t62 + t61;
t15 = -t20 * t47 + t45 * t50;
t12 = -t42 * t54 + t77;
t11 = t44 * t33 + t42 * t55;
t9 = t44 * t54 + t81;
t8 = t33 * t42 - t44 * t55;
t4 = t44 * t78 - t8 * t47;
t2 = t11 * t47 + t42 * t78;
t1 = [(-m(2) - m(3) - t73) * g(3) (-t89 * mrSges(3,1) - (t32 - t38 / 0.2e1) * mrSges(3,2) - m(4) * t75 - t86 * (t20 * pkin(3) + t75) + t84 * t21 - t85 * t20) * g(3) + (-(t24 * t44 - t80) * mrSges(3,1) - (t44 * t23 - t42 * t51) * mrSges(3,2) - m(4) * t69 - t86 * (-t8 * pkin(3) + t69) + t85 * t8 + t84 * (-t44 * t56 + t81)) * g(2) + (-(-t24 * t42 - t76) * mrSges(3,1) - (-t42 * t23 - t44 * t51) * mrSges(3,2) - m(4) * t64 - t86 * (-t11 * pkin(3) + t64) + t84 * (t42 * t56 + t77) + t85 * t11) * g(1) (-g(3) * t45 + (-g(1) * t42 + g(2) * t44) * t43) * t73, t86 * (-g(1) * t11 - g(2) * t8 + g(3) * t20) (-t90 * t15 + t88 * (-t20 * t50 - t45 * t47)) * g(3) + (t90 * t4 + t88 * (t44 * t79 + t50 * t8)) * g(2) + (-t90 * t2 + t88 * (t11 * t50 - t42 * t79)) * g(1), -g(1) * ((t12 * t49 - t2 * t46) * mrSges(7,1) + (-t12 * t46 - t2 * t49) * mrSges(7,2)) - g(2) * ((t4 * t46 + t49 * t9) * mrSges(7,1) + (t4 * t49 - t46 * t9) * mrSges(7,2)) - g(3) * ((-t15 * t46 + t21 * t49) * mrSges(7,1) + (-t15 * t49 - t21 * t46) * mrSges(7,2))];
taug  = t1(:);
