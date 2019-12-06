% Calculate Gravitation load on the joints for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:53:09
% EndTime: 2019-12-05 17:53:10
% DurationCPUTime: 0.27s
% Computational Cost: add. (199->50), mult. (165->49), div. (0->0), fcn. (121->10), ass. (0->30)
t18 = qJ(1) + pkin(8);
t11 = sin(t18);
t51 = t11 * g(2);
t47 = m(5) + m(6);
t46 = m(3) + m(4) + t47;
t50 = t46 * pkin(1) + mrSges(2,1);
t17 = qJ(3) + pkin(9);
t14 = qJ(5) + t17;
t7 = sin(t14);
t8 = cos(t14);
t49 = mrSges(6,1) * t7 + mrSges(6,2) * t8;
t10 = sin(t17);
t12 = cos(t17);
t20 = sin(qJ(3));
t31 = t8 * mrSges(6,1) - t7 * mrSges(6,2);
t48 = t12 * mrSges(5,1) - t20 * mrSges(4,2) - t10 * mrSges(5,2) + t31;
t22 = cos(qJ(3));
t32 = m(5) * pkin(3) + mrSges(4,1);
t45 = -t32 * t20 + m(6) * (-pkin(3) * t20 - pkin(4) * t10) - mrSges(5,1) * t10 - mrSges(4,2) * t22 - mrSges(5,2) * t12;
t15 = t22 * pkin(3);
t35 = pkin(4) * t12 + t15;
t44 = mrSges(3,1) + m(4) * pkin(2) + t22 * mrSges(4,1) + m(5) * (t15 + pkin(2)) + m(6) * (pkin(2) + t35) + t48;
t19 = -qJ(4) - pkin(6);
t42 = m(4) * pkin(6) - m(5) * t19 - m(6) * (-pkin(7) + t19) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t13 = cos(t18);
t37 = g(3) * t13;
t33 = t49 * t51;
t23 = cos(qJ(1));
t21 = sin(qJ(1));
t1 = [(mrSges(2,2) * t23 + t11 * t44 - t13 * t42 + t50 * t21) * g(3) + (-t21 * mrSges(2,2) + t11 * t42 + t13 * t44 + t50 * t23) * g(2), -t46 * g(1), -t33 + (-m(6) * t35 - t32 * t22 - t48) * g(1) + t45 * t51 + (t49 - t45) * t37, t47 * (-g(2) * t13 - g(3) * t11), -g(1) * t31 + t37 * t49 - t33];
taug = t1(:);
