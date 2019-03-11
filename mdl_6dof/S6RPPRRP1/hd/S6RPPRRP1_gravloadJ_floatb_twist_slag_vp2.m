% Calculate Gravitation load on the joints for
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:57:37
% EndTime: 2019-03-09 01:57:38
% DurationCPUTime: 0.57s
% Computational Cost: add. (363->75), mult. (324->84), div. (0->0), fcn. (278->10), ass. (0->42)
t54 = -mrSges(6,1) - mrSges(7,1);
t53 = mrSges(6,2) + mrSges(7,2);
t59 = mrSges(6,3) + mrSges(7,3);
t19 = sin(qJ(5));
t21 = cos(qJ(5));
t7 = t21 * pkin(5) + pkin(4);
t58 = -m(6) * pkin(4) - m(7) * t7 + t53 * t19 + t54 * t21;
t17 = -qJ(6) - pkin(8);
t57 = -m(6) * pkin(8) + m(7) * t17 - t59;
t14 = qJ(1) + pkin(9);
t11 = cos(t14);
t9 = sin(t14);
t56 = g(1) * t11 + g(2) * t9;
t47 = m(7) * pkin(5);
t55 = m(3) + m(4);
t52 = -m(5) - m(6) - m(7);
t51 = t47 - t54;
t16 = cos(pkin(10));
t13 = pkin(10) + qJ(4);
t10 = cos(t13);
t8 = sin(t13);
t29 = t10 * mrSges(5,1) - t8 * mrSges(5,2);
t49 = m(4) * pkin(2) + t16 * mrSges(4,1) - sin(pkin(10)) * mrSges(4,2) + mrSges(3,1) + t29 + t59 * t8;
t48 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t20 = sin(qJ(1));
t43 = t20 * pkin(1);
t22 = cos(qJ(1));
t12 = t22 * pkin(1);
t40 = t9 * t19;
t39 = t9 * t21;
t38 = t11 * t19;
t37 = t11 * t21;
t36 = m(4) - t52;
t32 = t10 * pkin(4) + t8 * pkin(8);
t31 = t10 * t7 - t8 * t17;
t1 = t10 * t40 + t37;
t3 = -t10 * t38 + t39;
t18 = -pkin(7) - qJ(3);
t6 = t16 * pkin(3) + pkin(2);
t4 = t10 * t37 + t40;
t2 = -t10 * t39 + t38;
t5 = [(-t40 * t47 - t22 * mrSges(2,1) + t20 * mrSges(2,2) + t54 * t4 + t52 * (t11 * t6 - t9 * t18 + t12) - t53 * t3 - t55 * t12 + t48 * t9 + (-m(6) * t32 - m(7) * t31 - t49) * t11) * g(2) + (-t38 * t47 + t20 * mrSges(2,1) + t22 * mrSges(2,2) + t55 * t43 + t52 * (-t11 * t18 - t43) + t54 * t2 - t53 * t1 + t48 * t11 + (m(5) * t6 - m(6) * (-t32 - t6) - m(7) * (-t31 - t6) + t49) * t9) * g(1) (-m(3) - t36) * g(3) (-g(1) * t9 + g(2) * t11) * t36, -g(3) * t29 + (t57 * g(3) + t56 * (mrSges(5,1) - t58)) * t8 + (t58 * g(3) + t56 * (mrSges(5,2) + t57)) * t10 (t51 * t19 + t53 * t21) * g(3) * t8 + (t51 * t1 - t53 * t2) * g(2) + (-t51 * t3 + t53 * t4) * g(1) (g(3) * t10 - t56 * t8) * m(7)];
taug  = t5(:);
