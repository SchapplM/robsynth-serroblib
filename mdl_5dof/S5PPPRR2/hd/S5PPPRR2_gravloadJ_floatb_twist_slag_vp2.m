% Calculate Gravitation load on the joints for
% S5PPPRR2
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
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPPRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:25
% EndTime: 2019-12-05 14:59:27
% DurationCPUTime: 0.22s
% Computational Cost: add. (99->37), mult. (251->68), div. (0->0), fcn. (272->10), ass. (0->27)
t18 = sin(qJ(5));
t20 = cos(qJ(5));
t32 = -m(6) * pkin(4) - t20 * mrSges(6,1) + t18 * mrSges(6,2) - mrSges(5,1);
t31 = -m(6) * pkin(6) + mrSges(5,2) - mrSges(6,3);
t12 = sin(pkin(9));
t13 = sin(pkin(8));
t30 = t12 * t13;
t19 = sin(qJ(4));
t29 = t13 * t19;
t21 = cos(qJ(4));
t28 = t13 * t21;
t14 = sin(pkin(7));
t16 = cos(pkin(8));
t27 = t14 * t16;
t17 = cos(pkin(7));
t26 = t16 * t17;
t25 = m(4) + m(5) + m(6);
t24 = m(3) + t25;
t15 = cos(pkin(9));
t10 = t15 * t28 - t16 * t19;
t8 = t12 * t14 + t15 * t26;
t7 = t12 * t26 - t14 * t15;
t6 = -t12 * t17 + t15 * t27;
t5 = t12 * t27 + t15 * t17;
t4 = t17 * t29 + t21 * t8;
t2 = t14 * t29 + t21 * t6;
t1 = [(-m(2) - t24) * g(3), (-g(1) * t14 + g(2) * t17) * t24, (t16 * g(3) + (-g(1) * t17 - g(2) * t14) * t13) * t25, (t32 * (-t15 * t29 - t16 * t21) + t31 * t10) * g(3) + (t31 * t2 + t32 * (t14 * t28 - t6 * t19)) * g(2) + (t31 * t4 + t32 * (t17 * t28 - t8 * t19)) * g(1), -g(1) * ((-t18 * t4 + t20 * t7) * mrSges(6,1) + (-t18 * t7 - t20 * t4) * mrSges(6,2)) - g(2) * ((-t18 * t2 + t20 * t5) * mrSges(6,1) + (-t18 * t5 - t2 * t20) * mrSges(6,2)) - g(3) * ((-t10 * t18 + t20 * t30) * mrSges(6,1) + (-t10 * t20 - t18 * t30) * mrSges(6,2))];
taug = t1(:);
