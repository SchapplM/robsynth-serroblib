% Calculate Gravitation load on the joints for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPPRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:49
% EndTime: 2019-12-05 14:57:51
% DurationCPUTime: 0.20s
% Computational Cost: add. (113->32), mult. (169->56), div. (0->0), fcn. (159->8), ass. (0->20)
t13 = sin(qJ(5));
t14 = cos(qJ(5));
t25 = m(6) * pkin(4) + t14 * mrSges(6,1) - t13 * mrSges(6,2) + mrSges(5,1);
t24 = -m(6) * pkin(6) + mrSges(5,2) - mrSges(6,3);
t9 = sin(pkin(8));
t22 = t13 * t9;
t21 = t14 * t9;
t10 = sin(pkin(7));
t11 = cos(pkin(8));
t20 = t10 * t11;
t12 = cos(pkin(7));
t19 = t11 * t12;
t18 = m(4) + m(5) + m(6);
t17 = m(3) + t18;
t8 = pkin(9) + qJ(4);
t7 = cos(t8);
t6 = sin(t8);
t4 = t10 * t6 + t7 * t19;
t2 = -t12 * t6 + t7 * t20;
t1 = [(-m(2) - t17) * g(3), (-g(1) * t10 + g(2) * t12) * t17, (g(3) * t11 + (-g(1) * t12 - g(2) * t10) * t9) * t18, (t24 * t7 + t25 * t6) * g(3) * t9 + (t24 * t2 - t25 * (-t12 * t7 - t6 * t20)) * g(2) + (t24 * t4 - t25 * (t10 * t7 - t6 * t19)) * g(1), -g(1) * ((t12 * t21 - t4 * t13) * mrSges(6,1) + (-t12 * t22 - t14 * t4) * mrSges(6,2)) - g(2) * ((t10 * t21 - t2 * t13) * mrSges(6,1) + (-t10 * t22 - t14 * t2) * mrSges(6,2)) - g(3) * ((-t11 * t14 - t7 * t22) * mrSges(6,1) + (t11 * t13 - t7 * t21) * mrSges(6,2))];
taug = t1(:);
