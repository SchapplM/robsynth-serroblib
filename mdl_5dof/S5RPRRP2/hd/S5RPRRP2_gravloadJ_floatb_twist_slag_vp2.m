% Calculate Gravitation load on the joints for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:01:23
% EndTime: 2019-12-05 18:01:24
% DurationCPUTime: 0.20s
% Computational Cost: add. (224->37), mult. (160->37), div. (0->0), fcn. (117->8), ass. (0->23)
t20 = cos(qJ(4));
t31 = mrSges(5,1) + mrSges(6,1);
t18 = sin(qJ(4));
t30 = mrSges(5,2) + mrSges(6,2);
t39 = t30 * t18;
t40 = m(5) * pkin(3) + m(6) * (pkin(4) * t20 + pkin(3)) + t31 * t20 + mrSges(4,1) - t39;
t37 = -m(5) * pkin(7) - mrSges(6,3) - mrSges(5,3) + mrSges(4,2) + m(6) * (-qJ(5) - pkin(7));
t29 = m(4) + m(5) + m(6);
t16 = qJ(1) + pkin(8);
t28 = m(3) + t29;
t27 = m(6) * pkin(4) + t31;
t26 = t29 * pkin(2) + mrSges(3,1);
t25 = t28 * pkin(1) + mrSges(2,1);
t15 = qJ(3) + t16;
t10 = sin(t15);
t11 = cos(t15);
t23 = -t37 * t10 + t40 * t11;
t22 = t40 * t10 + t37 * t11;
t21 = cos(qJ(1));
t19 = sin(qJ(1));
t14 = cos(t16);
t13 = sin(t16);
t1 = [(mrSges(2,2) * t21 + mrSges(3,2) * t14 + t26 * t13 + t25 * t19 + t22) * g(3) + (-t19 * mrSges(2,2) - t13 * mrSges(3,2) + t26 * t14 + t25 * t21 + t23) * g(2), -t28 * g(1), t23 * g(2) + t22 * g(3), (-t27 * t20 + t39) * g(1) + (-g(2) * t10 + g(3) * t11) * (t27 * t18 + t30 * t20), (-g(2) * t11 - g(3) * t10) * m(6)];
taug = t1(:);
