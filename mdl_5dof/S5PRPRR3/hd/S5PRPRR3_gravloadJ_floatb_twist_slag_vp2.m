% Calculate Gravitation load on the joints for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:46:37
% EndTime: 2019-12-05 15:46:39
% DurationCPUTime: 0.33s
% Computational Cost: add. (181->49), mult. (222->62), div. (0->0), fcn. (190->10), ass. (0->30)
t29 = m(4) + m(5) + m(6);
t50 = t29 * pkin(2) + mrSges(3,1);
t15 = sin(qJ(4));
t17 = cos(qJ(4));
t12 = qJ(4) + qJ(5);
t8 = sin(t12);
t9 = cos(t12);
t48 = m(5) * pkin(3) + t17 * mrSges(5,1) - t15 * mrSges(5,2) + mrSges(4,1) + m(6) * (pkin(4) * t17 + pkin(3)) + t9 * mrSges(6,1) - t8 * mrSges(6,2);
t47 = -m(5) * pkin(6) + m(6) * (-pkin(7) - pkin(6)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t46 = m(6) * pkin(4) + mrSges(5,1);
t11 = qJ(2) + pkin(9);
t6 = sin(t11);
t43 = g(3) * t6;
t14 = cos(pkin(8));
t34 = t14 * t9;
t35 = t14 * t8;
t13 = sin(pkin(8));
t36 = t13 * t9;
t37 = t13 * t8;
t7 = cos(t11);
t42 = (-t7 * t37 - t34) * mrSges(6,1) + (-t7 * t36 + t35) * mrSges(6,2);
t41 = (-t7 * t35 + t36) * mrSges(6,1) + (-t7 * t34 - t37) * mrSges(6,2);
t33 = t13 * t15;
t32 = t13 * t17;
t31 = t14 * t15;
t30 = t14 * t17;
t28 = -mrSges(6,1) * t8 - mrSges(6,2) * t9;
t18 = cos(qJ(2));
t16 = sin(qJ(2));
t1 = [(-m(2) - m(3) - t29) * g(3), (mrSges(3,2) * t16 - t50 * t18 + t47 * t6 - t48 * t7) * g(3) + (g(1) * t14 + g(2) * t13) * (mrSges(3,2) * t18 + t50 * t16 + t47 * t7 + t48 * t6), (-g(1) * t13 + g(2) * t14) * t29, (mrSges(5,2) * t17 + t15 * t46 - t28) * t43 + (-(-t7 * t32 + t31) * mrSges(5,2) - t42 - t46 * (-t7 * t33 - t30)) * g(2) + (-(-t7 * t30 - t33) * mrSges(5,2) - t41 - t46 * (-t7 * t31 + t32)) * g(1), -g(1) * t41 - g(2) * t42 - t28 * t43];
taug = t1(:);
