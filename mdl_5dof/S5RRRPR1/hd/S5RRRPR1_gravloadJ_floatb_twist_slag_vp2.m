% Calculate Gravitation load on the joints for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:37:47
% EndTime: 2019-12-05 18:37:48
% DurationCPUTime: 0.31s
% Computational Cost: add. (288->62), mult. (254->58), div. (0->0), fcn. (193->10), ass. (0->33)
t24 = sin(qJ(2));
t23 = qJ(2) + qJ(3);
t17 = pkin(9) + t23;
t12 = sin(t17);
t13 = cos(t17);
t18 = sin(t23);
t19 = cos(t23);
t15 = qJ(5) + t17;
t10 = cos(t15);
t9 = sin(t15);
t42 = t10 * mrSges(6,1) - t9 * mrSges(6,2);
t32 = -t19 * mrSges(4,1) - t13 * mrSges(5,1) + t18 * mrSges(4,2) + t12 * mrSges(5,2) - t42;
t60 = -t24 * mrSges(3,2) - t32;
t25 = sin(qJ(1));
t27 = cos(qJ(1));
t57 = g(1) * t27 + g(2) * t25;
t26 = cos(qJ(2));
t21 = t26 * pkin(2);
t14 = pkin(3) * t19;
t46 = t14 + t21;
t8 = pkin(4) * t13;
t45 = t8 + t46;
t56 = mrSges(2,1) + m(4) * (t21 + pkin(1)) + m(5) * (pkin(1) + t46) + m(3) * pkin(1) + t26 * mrSges(3,1) + m(6) * (pkin(1) + t45) + t60;
t28 = -pkin(7) - pkin(6);
t22 = -qJ(4) + t28;
t55 = mrSges(2,2) + m(6) * (-pkin(8) + t22) - mrSges(6,3) + m(5) * t22 - mrSges(5,3) + m(4) * t28 - mrSges(4,3) - m(3) * pkin(6) - mrSges(3,3);
t54 = pkin(3) * t18;
t51 = t24 * pkin(2);
t43 = m(4) * pkin(2) + mrSges(3,1);
t3 = -pkin(4) * t12 - t54;
t37 = mrSges(6,1) * t9 + mrSges(6,2) * t10;
t31 = mrSges(5,1) * t12 + mrSges(4,2) * t19 + mrSges(5,2) * t13 + t37;
t1 = [(t25 * t55 - t27 * t56) * g(2) + (t25 * t56 + t27 * t55) * g(1), (-m(5) * t46 - m(6) * t45 - t43 * t26 - t60) * g(3) + t57 * (-m(5) * (-t51 - t54) - m(6) * (t3 - t51) + mrSges(4,1) * t18 + mrSges(3,2) * t26 + t43 * t24 + t31), (-m(5) * t14 - m(6) * (t8 + t14) + t32) * g(3) + t57 * (-m(6) * t3 + (m(5) * pkin(3) + mrSges(4,1)) * t18 + t31), (m(5) + m(6)) * (-g(1) * t25 + g(2) * t27), -g(3) * t42 + t37 * t57];
taug = t1(:);
