% Calculate Gravitation load on the joints for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:04
% EndTime: 2019-12-05 17:06:04
% DurationCPUTime: 0.15s
% Computational Cost: add. (254->40), mult. (128->41), div. (0->0), fcn. (90->8), ass. (0->25)
t36 = mrSges(5,2) - mrSges(6,3);
t21 = sin(qJ(5));
t22 = cos(qJ(5));
t35 = t22 * mrSges(6,1) - t21 * mrSges(6,2);
t34 = -mrSges(5,1) - t35;
t33 = m(5) + m(6);
t20 = pkin(9) + qJ(2);
t19 = qJ(3) + t20;
t16 = qJ(4) + t19;
t11 = sin(t16);
t12 = cos(t16);
t32 = t12 * pkin(4) + t11 * pkin(8);
t15 = cos(t19);
t10 = pkin(3) * t15;
t29 = m(4) + t33;
t28 = t10 + t32;
t26 = t36 * t11 + t34 * t12;
t14 = sin(t19);
t25 = -t15 * mrSges(4,1) + mrSges(4,2) * t14 + t26;
t24 = (m(6) * pkin(4) - t34) * t11 + (-m(6) * pkin(8) + t36) * t12;
t23 = mrSges(4,2) * t15 + (t33 * pkin(3) + mrSges(4,1)) * t14 + t24;
t18 = cos(t20);
t17 = sin(t20);
t13 = pkin(2) * t18;
t1 = [(-m(2) - m(3) - t29) * g(3), (mrSges(3,2) * t17 - m(5) * (t10 + t13) - m(6) * (t13 + t28) + (-m(4) * pkin(2) - mrSges(3,1)) * t18 + t25) * g(2) + (mrSges(3,2) * t18 + (t29 * pkin(2) + mrSges(3,1)) * t17 + t23) * g(1), (-m(5) * t10 - m(6) * t28 + t25) * g(2) + t23 * g(1), (-m(6) * t32 + t26) * g(2) + t24 * g(1), -g(3) * t35 + (g(1) * t12 + g(2) * t11) * (mrSges(6,1) * t21 + mrSges(6,2) * t22)];
taug = t1(:);
