% Calculate Gravitation load on the joints for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:28
% EndTime: 2019-12-05 15:56:30
% DurationCPUTime: 0.51s
% Computational Cost: add. (289->61), mult. (520->93), div. (0->0), fcn. (582->12), ass. (0->34)
t28 = sin(qJ(5));
t30 = cos(qJ(5));
t57 = -m(6) * pkin(4) - t30 * mrSges(6,1) + t28 * mrSges(6,2) - mrSges(5,1);
t55 = -m(6) * pkin(8) + mrSges(5,2) - mrSges(6,3);
t51 = -m(4) * qJ(3) - t28 * mrSges(6,1) - t30 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t22 = pkin(10) + qJ(4);
t20 = sin(t22);
t21 = cos(t22);
t26 = cos(pkin(10));
t60 = -m(4) * pkin(2) - t26 * mrSges(4,1) + sin(pkin(10)) * mrSges(4,2) - mrSges(3,1) + t57 * t21 + t55 * t20;
t59 = m(5) + m(6);
t53 = m(4) + t59;
t24 = sin(pkin(9));
t25 = sin(pkin(5));
t47 = t24 * t25;
t29 = sin(qJ(2));
t46 = t25 * t29;
t31 = cos(qJ(2));
t45 = t25 * t31;
t44 = cos(pkin(5));
t43 = cos(pkin(9));
t40 = t24 * t44;
t39 = t25 * t43;
t36 = t44 * t43;
t27 = -pkin(7) - qJ(3);
t19 = t26 * pkin(3) + pkin(2);
t14 = -t29 * t40 + t31 * t43;
t13 = t29 * t43 + t31 * t40;
t12 = t24 * t31 + t29 * t36;
t11 = t24 * t29 - t31 * t36;
t8 = t20 * t44 + t21 * t46;
t4 = t14 * t21 + t20 * t47;
t2 = t12 * t21 - t20 * t39;
t1 = [(-m(2) - m(3) - t53) * g(3), (-t59 * (-t11 * t19 - t12 * t27) + t51 * t12 - t60 * t11) * g(2) + (-t59 * (-t13 * t19 - t14 * t27) + t51 * t14 - t60 * t13) * g(1) + (-t59 * t19 * t45 + (t60 * t31 + (t59 * t27 + t51) * t29) * t25) * g(3), t53 * (-g(1) * t13 - g(2) * t11 + g(3) * t45), (t55 * t8 + t57 * (-t20 * t46 + t21 * t44)) * g(3) + (t55 * t2 + t57 * (-t12 * t20 - t21 * t39)) * g(2) + (t55 * t4 + t57 * (-t14 * t20 + t21 * t47)) * g(1), -g(1) * ((t13 * t30 - t4 * t28) * mrSges(6,1) + (-t13 * t28 - t4 * t30) * mrSges(6,2)) - g(2) * ((t11 * t30 - t2 * t28) * mrSges(6,1) + (-t11 * t28 - t2 * t30) * mrSges(6,2)) - g(3) * ((-t8 * t28 - t30 * t45) * mrSges(6,1) + (t28 * t45 - t8 * t30) * mrSges(6,2))];
taug = t1(:);
