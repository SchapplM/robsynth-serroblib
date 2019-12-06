% Calculate Gravitation load on the joints for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:04
% EndTime: 2019-12-05 17:44:07
% DurationCPUTime: 0.41s
% Computational Cost: add. (216->66), mult. (272->79), div. (0->0), fcn. (249->10), ass. (0->42)
t39 = m(4) + m(5) + m(6);
t58 = -m(3) - t39;
t25 = sin(pkin(9));
t27 = cos(pkin(9));
t48 = t25 * pkin(3);
t24 = pkin(9) + qJ(4);
t18 = sin(t24);
t50 = pkin(4) * t18;
t59 = t58 * qJ(2) - m(5) * t48 - m(6) * (t48 + t50) + mrSges(2,2) - mrSges(3,3) - t25 * mrSges(4,1) - t27 * mrSges(4,2);
t57 = m(6) * pkin(4) + mrSges(5,1);
t17 = t27 * pkin(3) + pkin(2);
t19 = cos(t24);
t26 = sin(pkin(8));
t28 = cos(pkin(8));
t29 = -pkin(6) - qJ(3);
t55 = mrSges(2,1) + (mrSges(3,1) + m(4) * pkin(2) + t27 * mrSges(4,1) - t25 * mrSges(4,2) + m(5) * t17 + m(6) * (pkin(4) * t19 + t17)) * t28 - t58 * pkin(1) + (-mrSges(3,2) + m(4) * qJ(3) + mrSges(4,3) - m(5) * t29 + mrSges(5,3) - m(6) * (-pkin(7) + t29) + mrSges(6,3)) * t26;
t20 = qJ(5) + t24;
t16 = cos(t20);
t31 = cos(qJ(1));
t42 = t31 * t16;
t15 = sin(t20);
t30 = sin(qJ(1));
t47 = t30 * t15;
t5 = t28 * t47 + t42;
t43 = t31 * t15;
t46 = t30 * t16;
t6 = t28 * t46 - t43;
t52 = t5 * mrSges(6,1) + t6 * mrSges(6,2);
t7 = t28 * t43 - t46;
t8 = -t28 * t42 - t47;
t51 = -t7 * mrSges(6,1) + t8 * mrSges(6,2);
t49 = g(1) * t26;
t45 = t30 * t18;
t44 = t30 * t19;
t41 = t31 * t18;
t40 = t31 * t19;
t36 = -mrSges(6,1) * t15 - mrSges(6,2) * t16;
t11 = t28 * t41 - t44;
t9 = t28 * t45 + t40;
t12 = -t28 * t40 - t45;
t10 = t28 * t44 - t41;
t1 = [(t10 * mrSges(5,1) + t6 * mrSges(6,1) - t9 * mrSges(5,2) - t5 * mrSges(6,2) + t55 * t30 + t31 * t59) * g(3) + (-t12 * mrSges(5,1) - t8 * mrSges(6,1) - t11 * mrSges(5,2) - t7 * mrSges(6,2) - t30 * t59 + t55 * t31) * g(2), (g(2) * t31 + g(3) * t30) * t58, (g(1) * t28 + (g(2) * t30 - g(3) * t31) * t26) * t39, (m(6) * t50 + mrSges(5,1) * t18 + mrSges(5,2) * t19 - t36) * t49 + (-t12 * mrSges(5,2) + t57 * t11 - t51) * g(3) + (-t10 * mrSges(5,2) - t57 * t9 - t52) * g(2), -g(2) * t52 - g(3) * t51 - t36 * t49];
taug = t1(:);
