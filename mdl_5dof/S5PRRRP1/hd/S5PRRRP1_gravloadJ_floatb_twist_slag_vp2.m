% Calculate Gravitation load on the joints for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:39:47
% EndTime: 2019-12-05 16:39:49
% DurationCPUTime: 0.20s
% Computational Cost: add. (213->40), mult. (145->42), div. (0->0), fcn. (105->6), ass. (0->23)
t24 = cos(qJ(4));
t23 = sin(qJ(4));
t33 = mrSges(5,2) + mrSges(6,2);
t29 = t33 * t23;
t34 = mrSges(5,1) + mrSges(6,1);
t39 = -t34 * t24 - mrSges(4,1) + t29;
t37 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t21 = pkin(8) + qJ(2);
t20 = qJ(3) + t21;
t15 = sin(t20);
t16 = cos(t20);
t32 = t16 * pkin(3) + t15 * pkin(7);
t31 = m(4) + m(5) + m(6);
t17 = pkin(4) * t24 + pkin(3);
t22 = -qJ(5) - pkin(7);
t30 = -t15 * t22 + t16 * t17;
t28 = m(6) * pkin(4) + t34;
t26 = t37 * t15 + t39 * t16;
t25 = (-m(5) * pkin(7) + m(6) * t22 + t37) * t16 + (m(5) * pkin(3) + m(6) * t17 - t39) * t15;
t19 = cos(t21);
t18 = sin(t21);
t14 = pkin(2) * t19;
t1 = [(-m(2) - m(3) - t31) * g(3), (mrSges(3,2) * t18 - m(5) * (t14 + t32) - m(6) * (t14 + t30) + (-m(4) * pkin(2) - mrSges(3,1)) * t19 + t26) * g(2) + (mrSges(3,2) * t19 + (t31 * pkin(2) + mrSges(3,1)) * t18 + t25) * g(1), (-m(5) * t32 - m(6) * t30 + t26) * g(2) + t25 * g(1), (-t28 * t24 + t29) * g(3) + (g(1) * t16 + g(2) * t15) * (t28 * t23 + t33 * t24), (-g(1) * t15 + g(2) * t16) * m(6)];
taug = t1(:);
