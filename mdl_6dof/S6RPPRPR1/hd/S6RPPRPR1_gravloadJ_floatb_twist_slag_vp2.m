% Calculate Gravitation load on the joints for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Datum: 2018-11-23 15:39
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:39:47
% EndTime: 2018-11-23 15:39:48
% DurationCPUTime: 0.39s
% Computational Cost: add. (359->83), mult. (294->90), div. (0->0), fcn. (247->12), ass. (0->44)
t17 = qJ(1) + pkin(9);
t10 = sin(t17);
t13 = cos(t17);
t58 = g(1) * t13 + g(2) * t10;
t57 = -m(3) - m(4);
t56 = -m(5) - m(7);
t51 = m(6) + m(7);
t16 = pkin(10) + qJ(4);
t12 = cos(t16);
t18 = sin(pkin(11));
t20 = cos(pkin(11));
t30 = m(6) * pkin(4) + t20 * mrSges(6,1) - t18 * mrSges(6,2);
t55 = t30 * t12;
t21 = cos(pkin(10));
t9 = sin(t16);
t34 = t12 * mrSges(5,1) - t9 * mrSges(5,2);
t53 = mrSges(3,1) + m(4) * pkin(2) + t21 * mrSges(4,1) - sin(pkin(10)) * mrSges(4,2) + t34 + t9 * mrSges(7,3);
t23 = -pkin(7) - qJ(3);
t52 = -m(4) * qJ(3) + m(6) * t23 - t20 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t18;
t15 = pkin(11) + qJ(6);
t8 = sin(t15);
t47 = t13 * t8;
t24 = sin(qJ(1));
t46 = t24 * pkin(1);
t25 = cos(qJ(1));
t14 = t25 * pkin(1);
t7 = t21 * pkin(3) + pkin(2);
t44 = t13 * t7 + t14;
t43 = t10 * t12;
t11 = cos(t15);
t42 = t13 * t11;
t41 = m(4) + m(5) + t51;
t22 = -pkin(8) - qJ(5);
t40 = -m(7) * t22 + mrSges(7,3);
t39 = m(6) * qJ(5) + mrSges(6,3);
t6 = t20 * pkin(5) + pkin(4);
t36 = t12 * t6 - t9 * t22;
t32 = m(7) * t6 + t11 * mrSges(7,1) - t8 * mrSges(7,2);
t26 = t39 * t9 + t55;
t4 = t10 * t8 + t12 * t42;
t3 = t10 * t11 - t12 * t47;
t2 = -t11 * t43 + t47;
t1 = t43 * t8 + t42;
t5 = [(-m(6) * t44 - t25 * mrSges(2,1) - t4 * mrSges(7,1) + t24 * mrSges(2,2) - t3 * mrSges(7,2) + t56 * (-t10 * t23 + t44) + t57 * t14 + t52 * t10 + (-m(7) * t36 - t26 - t53) * t13) * g(2) + (t24 * mrSges(2,1) - t2 * mrSges(7,1) + t25 * mrSges(2,2) - t1 * mrSges(7,2) + t56 * (-t13 * t23 - t46) + (m(6) - t57) * t46 + t52 * t13 + (m(5) * t7 - m(6) * (-t9 * qJ(5) - t7) + t9 * mrSges(6,3) + t55 - m(7) * (-t36 - t7) + t53) * t10) * g(1) (-m(3) - t41) * g(3) (-g(1) * t10 + g(2) * t13) * t41 (-t26 - t34) * g(3) + (-t40 * g(3) + t58 * (mrSges(5,1) + t30 + t32)) * t9 + (-t32 * g(3) + t58 * (mrSges(5,2) - t39 - t40)) * t12 (t12 * g(3) - t9 * t58) * t51, -g(1) * (t3 * mrSges(7,1) - t4 * mrSges(7,2)) - g(2) * (-t1 * mrSges(7,1) + t2 * mrSges(7,2)) - g(3) * (-mrSges(7,1) * t8 - mrSges(7,2) * t11) * t9];
taug  = t5(:);
