% Calculate Gravitation load on the joints for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:14
% EndTime: 2019-12-05 17:51:15
% DurationCPUTime: 0.32s
% Computational Cost: add. (261->58), mult. (185->69), div. (0->0), fcn. (158->10), ass. (0->34)
t52 = mrSges(4,2) - mrSges(5,3);
t20 = sin(pkin(9));
t21 = cos(pkin(9));
t51 = mrSges(4,1) - m(6) * (-pkin(4) * t21 - pkin(3)) + t21 * mrSges(5,1) + (m(6) * pkin(7) - mrSges(5,2) + mrSges(6,3)) * t20;
t49 = m(5) + m(6);
t19 = qJ(1) + pkin(8);
t18 = qJ(3) + t19;
t14 = sin(t18);
t15 = cos(t18);
t24 = cos(qJ(5));
t22 = sin(qJ(5));
t40 = t21 * t22;
t7 = -t14 * t24 + t15 * t40;
t39 = t21 * t24;
t8 = -t14 * t22 - t15 * t39;
t48 = -t8 * mrSges(6,1) - t7 * mrSges(6,2) - t52 * t14 + t51 * t15;
t5 = t14 * t40 + t15 * t24;
t6 = t14 * t39 - t15 * t22;
t47 = t6 * mrSges(6,1) - t5 * mrSges(6,2) + t51 * t14 + t52 * t15;
t23 = sin(qJ(1));
t46 = pkin(1) * t23;
t25 = cos(qJ(1));
t45 = pkin(1) * t25;
t44 = pkin(3) * t14;
t43 = t15 * pkin(3);
t38 = t14 * qJ(4);
t16 = sin(t19);
t36 = -pkin(2) * t16 - t46;
t17 = cos(t19);
t35 = -pkin(2) * t17 - t45;
t11 = t15 * qJ(4);
t31 = t11 + t36;
t29 = t35 - t38;
t1 = [(t23 * mrSges(2,1) + mrSges(2,2) * t25 + m(3) * t46 + mrSges(3,1) * t16 + mrSges(3,2) * t17 - m(4) * t36 - m(5) * (t31 - t44) - m(6) * t31 + t47) * g(3) + (mrSges(2,1) * t25 - t23 * mrSges(2,2) + m(3) * t45 + t17 * mrSges(3,1) - t16 * mrSges(3,2) - m(4) * t35 - m(5) * (t29 - t43) - m(6) * t29 + t48) * g(2), (-m(3) - m(4) - t49) * g(1), (-m(5) * (t11 - t44) - m(6) * t11 + t47) * g(3) + (-m(5) * (-t38 - t43) + m(6) * t38 + t48) * g(2), t49 * (-g(2) * t15 - g(3) * t14), -g(2) * (mrSges(6,1) * t5 + mrSges(6,2) * t6) - g(3) * (-mrSges(6,1) * t7 + mrSges(6,2) * t8) - g(1) * (-mrSges(6,1) * t22 - mrSges(6,2) * t24) * t20];
taug = t1(:);
