% Calculate Gravitation load on the joints for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
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
% Datum: 2018-11-23 15:52
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:52:11
% EndTime: 2018-11-23 15:52:12
% DurationCPUTime: 0.65s
% Computational Cost: add. (380->87), mult. (321->92), div. (0->0), fcn. (270->12), ass. (0->48)
t19 = sin(pkin(11));
t20 = cos(pkin(11));
t69 = -m(6) * pkin(4) - t20 * mrSges(6,1) + t19 * mrSges(6,2) - mrSges(5,1);
t56 = m(6) + m(7);
t42 = m(5) + t56;
t68 = mrSges(5,2) - mrSges(7,3);
t18 = qJ(1) + pkin(9);
t10 = sin(t18);
t13 = cos(t18);
t64 = g(1) * t13 + g(2) * t10;
t17 = qJ(3) + pkin(10);
t9 = sin(t17);
t67 = t9 * t64;
t12 = cos(t17);
t23 = sin(qJ(3));
t25 = cos(qJ(3));
t66 = -t25 * mrSges(4,1) + t23 * mrSges(4,2) + t69 * t12 + t68 * t9;
t65 = m(6) * qJ(5) + mrSges(6,3);
t63 = -m(4) * pkin(2) - mrSges(3,1) + t66;
t61 = -m(3) - m(4);
t60 = -m(5) - m(7);
t21 = -qJ(4) - pkin(7);
t57 = -m(4) * pkin(7) + m(6) * t21 - t20 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t19;
t16 = pkin(11) + qJ(6);
t8 = sin(t16);
t51 = t13 * t8;
t24 = sin(qJ(1));
t50 = t24 * pkin(1);
t14 = t25 * pkin(3);
t26 = cos(qJ(1));
t15 = t26 * pkin(1);
t49 = t9 * mrSges(6,3);
t22 = -pkin(8) - qJ(5);
t47 = t9 * t22;
t7 = t14 + pkin(2);
t46 = t13 * t7 + t15;
t45 = t10 * t12;
t11 = cos(t16);
t44 = t13 * t11;
t43 = t9 * qJ(5);
t6 = t20 * pkin(5) + pkin(4);
t39 = t12 * t6 - t47;
t34 = m(7) * t6 + t11 * mrSges(7,1) - t8 * mrSges(7,2);
t4 = t10 * t8 + t12 * t44;
t3 = t10 * t11 - t12 * t51;
t2 = -t11 * t45 + t51;
t1 = t8 * t45 + t44;
t5 = [(-m(6) * t46 - t26 * mrSges(2,1) - t4 * mrSges(7,1) + t24 * mrSges(2,2) - t3 * mrSges(7,2) + t60 * (-t10 * t21 + t46) + t61 * t15 + t57 * t10 + (-m(7) * t39 - t65 * t9 + t63) * t13) * g(2) + (t24 * mrSges(2,1) - t2 * mrSges(7,1) + t26 * mrSges(2,2) - t1 * mrSges(7,2) + t60 * (-t13 * t21 - t50) + (m(6) - t61) * t50 + t57 * t13 + (m(5) * t7 - m(6) * (-t7 - t43) + t49 - m(7) * (-t39 - t7) - t63) * t10) * g(1) (-t42 + t61) * g(3) (t34 - t69) * t67 + (-m(5) * t14 - m(6) * (t14 + t43) - t49 - m(7) * (t14 - t47) + t66 - t34 * t12) * g(3) + (mrSges(4,2) * t25 + (m(7) * t22 - t65 + t68) * t12 + (t42 * pkin(3) + mrSges(4,1)) * t23) * t64 (-g(1) * t10 + g(2) * t13) * t42 (t12 * g(3) - t67) * t56, -g(1) * (t3 * mrSges(7,1) - t4 * mrSges(7,2)) - g(2) * (-t1 * mrSges(7,1) + t2 * mrSges(7,2)) - g(3) * (-mrSges(7,1) * t8 - mrSges(7,2) * t11) * t9];
taug  = t5(:);
