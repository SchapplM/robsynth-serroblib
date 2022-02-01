% Calculate Gravitation load on the joints for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m [6x1]
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
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:08
% EndTime: 2022-01-23 09:14:09
% DurationCPUTime: 0.22s
% Computational Cost: add. (185->41), mult. (145->41), div. (0->0), fcn. (105->10), ass. (0->25)
t14 = pkin(9) + qJ(4);
t10 = qJ(5) + t14;
t3 = sin(t10);
t4 = cos(t10);
t26 = t4 * mrSges(6,1) - mrSges(6,2) * t3;
t6 = sin(t14);
t40 = -t6 * mrSges(5,2) + t26;
t28 = m(4) + m(5) + m(6);
t38 = -m(3) - t28;
t39 = pkin(1) * t38 - mrSges(2,1);
t15 = qJ(1) + pkin(8);
t7 = sin(t15);
t9 = cos(t15);
t37 = g(1) * t9 + g(2) * t7;
t17 = cos(pkin(9));
t5 = t17 * pkin(3) + pkin(2);
t8 = cos(t14);
t35 = mrSges(3,1) + m(4) * pkin(2) + t17 * mrSges(4,1) - sin(pkin(9)) * mrSges(4,2) + m(5) * t5 + t8 * mrSges(5,1) + m(6) * (pkin(4) * t8 + t5) + t40;
t18 = -pkin(6) - qJ(3);
t34 = -m(4) * qJ(3) + m(5) * t18 - m(6) * (pkin(7) - t18) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t27 = m(6) * pkin(4) + mrSges(5,1);
t25 = mrSges(6,1) * t3 + mrSges(6,2) * t4;
t20 = cos(qJ(1));
t19 = sin(qJ(1));
t1 = [(t19 * mrSges(2,2) + t39 * t20 + t34 * t7 - t35 * t9) * g(2) + (mrSges(2,2) * t20 - t39 * t19 + t34 * t9 + t35 * t7) * g(1), t38 * g(3), (-g(1) * t7 + g(2) * t9) * t28, (-t27 * t8 - t40) * g(3) + t37 * (mrSges(5,2) * t8 + t27 * t6 + t25), -g(3) * t26 + t37 * t25];
taug = t1(:);
