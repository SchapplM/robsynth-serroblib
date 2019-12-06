% Calculate Gravitation load on the joints for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(1,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:08:38
% EndTime: 2019-12-05 18:08:39
% DurationCPUTime: 0.33s
% Computational Cost: add. (122->56), mult. (297->82), div. (0->0), fcn. (293->8), ass. (0->35)
t27 = mrSges(5,2) - mrSges(6,3);
t28 = mrSges(4,2) - mrSges(5,3);
t10 = sin(qJ(4));
t14 = cos(qJ(4));
t44 = -mrSges(5,1) * t14 + t27 * t10;
t11 = sin(qJ(3));
t43 = t28 * t11;
t9 = sin(qJ(5));
t39 = t11 * t9;
t13 = cos(qJ(5));
t37 = t11 * t13;
t36 = t11 * t14;
t16 = cos(qJ(1));
t35 = t11 * t16;
t12 = sin(qJ(1));
t34 = t12 * t10;
t15 = cos(qJ(3));
t33 = t14 * t15;
t32 = t15 * mrSges(4,1);
t31 = t15 * t13;
t30 = t16 * t10;
t29 = t16 * t14;
t26 = m(3) + m(4) + m(5) + m(6);
t25 = -mrSges(2,1) - mrSges(3,1) - t32;
t22 = t9 * t36 + t31;
t21 = -t13 * t36 + t15 * t9;
t20 = t13 * mrSges(6,1) - t9 * mrSges(6,2) + mrSges(5,1);
t19 = -t26 * qJ(2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t6 = t15 * t29 + t34;
t5 = -t12 * t14 + t15 * t30;
t4 = t12 * t33 - t30;
t3 = t15 * t34 + t29;
t2 = t6 * t13 + t9 * t35;
t1 = t13 * t35 - t6 * t9;
t7 = [(-t6 * mrSges(5,1) - t2 * mrSges(6,1) - t1 * mrSges(6,2) + t27 * t5 + (t25 + t43) * t16 + t19 * t12) * g(2) + (-t27 * t3 + t20 * t4 + ((t9 * mrSges(6,1) + t13 * mrSges(6,2) - t28) * t11 - t25) * t12 + t19 * t16) * g(1), (-t12 * g(1) + t16 * g(2)) * t26, (-t32 - (t14 * t31 + t39) * mrSges(6,1) - (-t9 * t33 + t37) * mrSges(6,2) + t44 * t15 + t43) * g(3) + (g(1) * t16 + g(2) * t12) * (-t21 * mrSges(6,1) - t22 * mrSges(6,2) + t28 * t15 + (mrSges(4,1) - t44) * t11), (t20 * t3 + t27 * t4) * g(2) + (t20 * t5 + t27 * t6) * g(1) + (t20 * t10 + t27 * t14) * g(3) * t11, -g(1) * (t1 * mrSges(6,1) - mrSges(6,2) * t2) - g(2) * ((t12 * t37 - t4 * t9) * mrSges(6,1) + (-t12 * t39 - t4 * t13) * mrSges(6,2)) - g(3) * (-t22 * mrSges(6,1) + t21 * mrSges(6,2))];
taug = t7(:);
