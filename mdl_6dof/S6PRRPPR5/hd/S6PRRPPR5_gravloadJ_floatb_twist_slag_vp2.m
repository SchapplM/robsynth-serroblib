% Calculate Gravitation load on the joints for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
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
% Datum: 2018-11-23 15:10
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRPPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:10:39
% EndTime: 2018-11-23 15:10:39
% DurationCPUTime: 0.83s
% Computational Cost: add. (1027->95), mult. (1221->117), div. (0->0), fcn. (1180->16), ass. (0->62)
t102 = m(6) + m(7);
t98 = m(5) + t102;
t114 = -m(6) * qJ(5) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3) + m(7) * (-pkin(9) - qJ(5)) - mrSges(7,3);
t37 = pkin(11) + qJ(6);
t35 = sin(t37);
t36 = cos(t37);
t38 = sin(pkin(11));
t40 = cos(pkin(11));
t113 = -t35 * mrSges(7,1) - t40 * mrSges(6,2) - t36 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t38;
t42 = sin(qJ(3));
t44 = cos(qJ(3));
t112 = pkin(3) * t44 + qJ(4) * t42;
t105 = t113 * t42 + t114 * t44 - mrSges(3,1);
t90 = pkin(4) + pkin(8);
t104 = m(7) * (pkin(5) * t40 + t90) + t36 * mrSges(7,1) - t35 * mrSges(7,2) + m(6) * t90 + t40 * mrSges(6,1) - t38 * mrSges(6,2) + mrSges(5,1) - mrSges(3,2) + mrSges(4,3);
t103 = pkin(3) * t98 - t114;
t39 = sin(pkin(10));
t43 = sin(qJ(2));
t76 = pkin(6) + qJ(2);
t63 = cos(t76) / 0.2e1;
t77 = pkin(6) - qJ(2);
t66 = cos(t77);
t49 = t66 / 0.2e1 + t63;
t79 = cos(pkin(10));
t14 = t39 * t43 - t49 * t79;
t101 = t112 * t14;
t17 = t39 * t49 + t43 * t79;
t100 = t112 * t17;
t64 = sin(t76);
t61 = t64 / 0.2e1;
t65 = sin(t77);
t62 = t65 / 0.2e1;
t27 = t61 + t62;
t99 = t112 * t27;
t96 = -t98 * qJ(4) + t113;
t45 = cos(qJ(2));
t85 = t39 * t45;
t80 = cos(pkin(6));
t78 = sin(pkin(6));
t11 = t14 * pkin(2);
t48 = t62 - t64 / 0.2e1;
t16 = -t48 * t79 + t85;
t71 = pkin(8) * t16 - t11;
t12 = t17 * pkin(2);
t67 = t79 * t45;
t19 = t39 * t48 + t67;
t70 = pkin(8) * t19 - t12;
t26 = t27 * pkin(2);
t28 = t63 - t66 / 0.2e1;
t69 = -pkin(8) * t28 + t26;
t68 = t39 * t78;
t56 = t79 * t78;
t55 = t61 - t65 / 0.2e1;
t21 = -t28 * t44 + t42 * t80;
t20 = -t28 * t42 - t44 * t80;
t18 = -t39 * t55 + t67;
t15 = t55 * t79 + t85;
t6 = t18 * t44 + t42 * t68;
t5 = t18 * t42 - t44 * t68;
t4 = t15 * t44 - t42 * t56;
t3 = t15 * t42 + t44 * t56;
t1 = [(-m(2) - m(3) - m(4) - t98) * g(3) (-m(4) * t69 - m(5) * (t69 + t99) - t102 * (t26 + t99) + t104 * t28 + t105 * t27) * g(3) + (-m(4) * t71 - m(5) * (t71 - t101) - t102 * (-t11 - t101) - t104 * t16 - t105 * t14) * g(2) + (-m(4) * t70 - m(5) * (t70 - t100) - t102 * (-t12 - t100) - t104 * t19 - t105 * t17) * g(1) (t103 * t20 + t96 * t21) * g(3) + (t103 * t3 + t96 * t4) * g(2) + (t103 * t5 + t96 * t6) * g(1), t98 * (-g(1) * t5 - g(2) * t3 - g(3) * t20) t102 * (-g(1) * t6 - g(2) * t4 - g(3) * t21) -g(1) * ((-t17 * t35 + t36 * t5) * mrSges(7,1) + (-t17 * t36 - t35 * t5) * mrSges(7,2)) - g(2) * ((-t14 * t35 + t3 * t36) * mrSges(7,1) + (-t14 * t36 - t3 * t35) * mrSges(7,2)) - g(3) * ((t20 * t36 + t27 * t35) * mrSges(7,1) + (-t20 * t35 + t27 * t36) * mrSges(7,2))];
taug  = t1(:);
