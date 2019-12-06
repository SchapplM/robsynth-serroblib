% Calculate Gravitation load on the joints for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:10
% EndTime: 2019-12-05 18:18:11
% DurationCPUTime: 0.20s
% Computational Cost: add. (242->37), mult. (156->37), div. (0->0), fcn. (114->10), ass. (0->22)
t19 = pkin(9) + qJ(5);
t14 = sin(t19);
t15 = cos(t19);
t42 = mrSges(6,1) * t15 - mrSges(6,2) * t14;
t22 = cos(pkin(9));
t41 = -sin(pkin(9)) * mrSges(5,2) + m(5) * pkin(3) + m(6) * (pkin(4) * t22 + pkin(3)) + t22 * mrSges(5,1) + mrSges(4,1) + t42;
t40 = m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3) - m(6) * (-pkin(7) - qJ(4));
t39 = m(5) + m(6);
t20 = qJ(1) + qJ(2);
t33 = m(4) + t39;
t30 = pkin(2) * t33 + mrSges(3,1);
t29 = mrSges(2,1) + (m(3) + t33) * pkin(1);
t16 = pkin(8) + t20;
t11 = sin(t16);
t12 = cos(t16);
t17 = sin(t20);
t18 = cos(t20);
t27 = -t17 * mrSges(3,2) + t40 * t11 + t41 * t12 + t30 * t18;
t26 = mrSges(3,2) * t18 + t41 * t11 - t40 * t12 + t30 * t17;
t25 = cos(qJ(1));
t24 = sin(qJ(1));
t1 = [(mrSges(2,2) * t25 + t24 * t29 + t26) * g(3) + (-t24 * mrSges(2,2) + t25 * t29 + t27) * g(2), g(2) * t27 + g(3) * t26, -t33 * g(1), t39 * (-g(2) * t12 - g(3) * t11), -g(1) * t42 + (-g(2) * t11 + g(3) * t12) * (mrSges(6,1) * t14 + mrSges(6,2) * t15)];
taug = t1(:);
