% Calculate Gravitation load on the joints for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:21:59
% EndTime: 2019-12-05 18:22:00
% DurationCPUTime: 0.24s
% Computational Cost: add. (243->37), mult. (173->37), div. (0->0), fcn. (127->8), ass. (0->22)
t21 = cos(qJ(4));
t31 = mrSges(5,1) + mrSges(6,1);
t19 = sin(qJ(4));
t30 = mrSges(5,2) + mrSges(6,2);
t39 = t30 * t19;
t40 = m(5) * pkin(3) + m(6) * (pkin(4) * t21 + pkin(3)) + t31 * t21 + mrSges(4,1) - t39;
t37 = -m(5) * pkin(7) - mrSges(6,3) - mrSges(5,3) + mrSges(4,2) + m(6) * (-qJ(5) - pkin(7));
t17 = qJ(1) + qJ(2);
t29 = m(4) + m(5) + m(6);
t28 = m(6) * pkin(4) + t31;
t27 = t29 * pkin(2) + mrSges(3,1);
t26 = mrSges(2,1) + (m(3) + t29) * pkin(1);
t14 = pkin(8) + t17;
t11 = sin(t14);
t12 = cos(t14);
t15 = sin(t17);
t16 = cos(t17);
t24 = -t15 * mrSges(3,2) - t37 * t11 + t40 * t12 + t27 * t16;
t23 = mrSges(3,2) * t16 + t40 * t11 + t37 * t12 + t27 * t15;
t22 = cos(qJ(1));
t20 = sin(qJ(1));
t1 = [(mrSges(2,2) * t22 + t26 * t20 + t23) * g(3) + (-t20 * mrSges(2,2) + t26 * t22 + t24) * g(2), t24 * g(2) + t23 * g(3), -t29 * g(1), (-t28 * t21 + t39) * g(1) + (-g(2) * t11 + g(3) * t12) * (t28 * t19 + t30 * t21), (-g(2) * t12 - g(3) * t11) * m(6)];
taug = t1(:);
