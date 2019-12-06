% Calculate Gravitation load on the joints for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:21:29
% EndTime: 2019-12-05 16:21:31
% DurationCPUTime: 0.43s
% Computational Cost: add. (206->70), mult. (281->91), div. (0->0), fcn. (247->10), ass. (0->32)
t17 = qJ(3) + pkin(9);
t14 = qJ(5) + t17;
t10 = cos(t14);
t12 = sin(t17);
t13 = cos(t17);
t23 = cos(qJ(3));
t15 = t23 * pkin(3);
t21 = sin(qJ(3));
t7 = pkin(4) * t13 + t15;
t9 = sin(t14);
t56 = mrSges(3,1) + m(5) * (t15 + pkin(2)) + t13 * mrSges(5,1) - t12 * mrSges(5,2) + m(4) * pkin(2) + t23 * mrSges(4,1) - t21 * mrSges(4,2) + m(6) * (pkin(2) + t7) + t10 * mrSges(6,1) - t9 * mrSges(6,2);
t20 = -qJ(4) - pkin(6);
t55 = mrSges(3,2) + m(6) * (-pkin(7) + t20) - mrSges(6,3) + m(5) * t20 - mrSges(5,3) - m(4) * pkin(6) - mrSges(4,3);
t18 = sin(pkin(8));
t19 = cos(pkin(8));
t54 = g(1) * t19 + g(2) * t18;
t50 = m(5) + m(6);
t53 = -m(5) * pkin(3) - mrSges(4,1);
t24 = cos(qJ(2));
t39 = t24 * t10;
t43 = t18 * t24;
t49 = (-t19 * t10 - t9 * t43) * mrSges(6,1) + (-t18 * t39 + t19 * t9) * mrSges(6,2);
t42 = t19 * t24;
t48 = (t18 * t10 - t9 * t42) * mrSges(6,1) + (-t18 * t9 - t19 * t39) * mrSges(6,2);
t22 = sin(qJ(2));
t45 = g(3) * t22;
t44 = t21 * pkin(3);
t41 = t21 * t24;
t40 = t23 * t24;
t34 = -mrSges(6,1) * t9 - mrSges(6,2) * t10;
t6 = -pkin(4) * t12 - t44;
t1 = [(-m(2) - m(3) - m(4) - t50) * g(3), (-t56 * g(3) + t54 * t55) * t24 + (t55 * g(3) + t54 * t56) * t22, (m(5) * t44 - m(6) * t6 + mrSges(4,1) * t21 + mrSges(5,1) * t12 + mrSges(4,2) * t23 + mrSges(5,2) * t13 - t34) * t45 + (-(-t18 * t40 + t19 * t21) * mrSges(4,2) - (-t12 * t43 - t19 * t13) * mrSges(5,1) - (t19 * t12 - t13 * t43) * mrSges(5,2) - m(6) * (-t19 * t7 + t6 * t43) - t49 + t53 * (-t18 * t41 - t19 * t23)) * g(2) + (-(-t18 * t21 - t19 * t40) * mrSges(4,2) - (-t12 * t42 + t18 * t13) * mrSges(5,1) - (-t18 * t12 - t13 * t42) * mrSges(5,2) - m(6) * (t18 * t7 + t6 * t42) - t48 + t53 * (t18 * t23 - t19 * t41)) * g(1), (t24 * g(3) - t22 * t54) * t50, -g(1) * t48 - g(2) * t49 - t34 * t45];
taug = t1(:);
