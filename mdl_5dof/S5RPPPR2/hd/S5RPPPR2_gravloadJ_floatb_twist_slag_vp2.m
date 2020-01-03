% Calculate Gravitation load on the joints for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:20
% EndTime: 2020-01-03 11:22:21
% DurationCPUTime: 0.39s
% Computational Cost: add. (166->71), mult. (381->103), div. (0->0), fcn. (411->10), ass. (0->42)
t26 = sin(pkin(8));
t31 = sin(qJ(1));
t45 = cos(pkin(8));
t42 = t31 * t45;
t29 = cos(pkin(7));
t33 = cos(qJ(1));
t48 = t33 * t29;
t13 = t26 * t48 - t42;
t30 = sin(qJ(5));
t32 = cos(qJ(5));
t41 = t33 * t45;
t14 = t31 * t26 + t29 * t41;
t25 = sin(pkin(9));
t28 = cos(pkin(9));
t27 = sin(pkin(7));
t49 = t33 * t27;
t6 = t14 * t28 + t25 * t49;
t61 = -t13 * t32 + t6 * t30;
t60 = t13 * t30 + t6 * t32;
t59 = m(5) + m(6);
t58 = -t29 * mrSges(3,1) + t27 * mrSges(3,2) - mrSges(2,1);
t57 = mrSges(2,2) - mrSges(3,3);
t56 = mrSges(4,2) - mrSges(5,3);
t55 = -m(6) * pkin(6) + mrSges(5,2) - mrSges(6,3);
t52 = t26 * t27;
t51 = t31 * t27;
t50 = t31 * t29;
t47 = t33 * pkin(1) + t31 * qJ(2);
t46 = qJ(3) * t27;
t44 = m(4) + t59;
t40 = t31 * pkin(1) - t33 * qJ(2);
t39 = pkin(2) * t48 + t33 * t46 + t47;
t37 = pkin(2) * t50 + t31 * t46 + t40;
t36 = t14 * pkin(3) + t13 * qJ(4) + t39;
t11 = t26 * t50 + t41;
t12 = -t33 * t26 + t29 * t42;
t34 = t12 * pkin(3) + t11 * qJ(4) + t37;
t10 = t27 * t45 * t28 - t29 * t25;
t4 = t12 * t28 + t25 * t51;
t2 = t11 * t30 + t4 * t32;
t1 = t11 * t32 - t4 * t30;
t3 = [(-m(3) * t40 - m(4) * t37 - t12 * mrSges(4,1) - mrSges(4,3) * t51 - m(5) * t34 - t4 * mrSges(5,1) - m(6) * (t4 * pkin(4) + t34) - t2 * mrSges(6,1) - t1 * mrSges(6,2) - t57 * t33 + t58 * t31 + t55 * (t12 * t25 - t28 * t51) + t56 * t11) * g(3) + (-m(3) * t47 - m(4) * t39 - t14 * mrSges(4,1) - mrSges(4,3) * t49 - m(5) * t36 - t6 * mrSges(5,1) - m(6) * (t6 * pkin(4) + t36) - t60 * mrSges(6,1) + t61 * mrSges(6,2) + t55 * (t14 * t25 - t28 * t49) + t58 * t33 + t57 * t31 + t56 * t13) * g(2), (g(2) * t33 + g(3) * t31) * (m(3) + t44), (g(1) * t29 + (-g(2) * t31 + g(3) * t33) * t27) * t44, t59 * (-g(1) * t52 - g(2) * t11 + g(3) * t13), -g(1) * ((-t10 * t30 + t32 * t52) * mrSges(6,1) + (-t10 * t32 - t30 * t52) * mrSges(6,2)) - g(2) * (t1 * mrSges(6,1) - t2 * mrSges(6,2)) - g(3) * (t61 * mrSges(6,1) + t60 * mrSges(6,2))];
taug = t3(:);
