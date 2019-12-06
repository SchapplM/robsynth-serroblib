% Calculate Gravitation load on the joints for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:58:09
% EndTime: 2019-12-05 18:58:10
% DurationCPUTime: 0.26s
% Computational Cost: add. (348->42), mult. (223->44), div. (0->0), fcn. (169->10), ass. (0->30)
t24 = sin(qJ(4));
t22 = qJ(4) + qJ(5);
t17 = sin(t22);
t19 = cos(t22);
t37 = t19 * mrSges(6,1) - t17 * mrSges(6,2);
t53 = -t24 * mrSges(5,2) + t37;
t52 = mrSges(6,1) * t17 + mrSges(6,2) * t19;
t26 = cos(qJ(4));
t51 = m(5) * pkin(3) + m(6) * (t26 * pkin(4) + pkin(3)) + t26 * mrSges(5,1) + mrSges(4,1) + t53;
t50 = -m(5) * pkin(8) - mrSges(5,3) - mrSges(6,3) + mrSges(4,2) + m(6) * (-pkin(9) - pkin(8));
t38 = m(6) * pkin(4) + mrSges(5,1);
t47 = -mrSges(5,2) * t26 - t24 * t38;
t23 = qJ(1) + qJ(2);
t21 = qJ(3) + t23;
t14 = sin(t21);
t46 = t52 * t14;
t15 = cos(t21);
t45 = g(3) * t15;
t39 = m(4) + m(5) + m(6);
t35 = pkin(2) * t39 + mrSges(3,1);
t34 = mrSges(2,1) + (m(3) + t39) * pkin(1);
t32 = -t50 * t14 + t51 * t15;
t31 = t51 * t14 + t50 * t15;
t18 = sin(t23);
t20 = cos(t23);
t30 = -t18 * mrSges(3,2) + t20 * t35 + t32;
t29 = t20 * mrSges(3,2) + t18 * t35 + t31;
t27 = cos(qJ(1));
t25 = sin(qJ(1));
t1 = [(t27 * mrSges(2,2) + t25 * t34 + t29) * g(3) + (-t25 * mrSges(2,2) + t27 * t34 + t30) * g(2), g(2) * t30 + g(3) * t29, g(2) * t32 + g(3) * t31, (t47 * t14 - t46) * g(2) + (-t26 * t38 - t53) * g(1) + (t52 - t47) * t45, -g(1) * t37 - g(2) * t46 + t45 * t52];
taug = t1(:);
