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
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:55:43
% EndTime: 2020-01-03 11:55:44
% DurationCPUTime: 0.20s
% Computational Cost: add. (242->55), mult. (156->56), div. (0->0), fcn. (114->10), ass. (0->31)
t31 = pkin(9) + qJ(5);
t24 = sin(t31);
t25 = cos(t31);
t54 = mrSges(6,1) * t25 - mrSges(6,2) * t24;
t53 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t34 = cos(pkin(9));
t52 = sin(pkin(9)) * mrSges(5,2) - mrSges(4,1) - mrSges(5,1) * t34 - t54;
t51 = m(5) + m(6);
t32 = qJ(1) + qJ(2);
t27 = sin(t32);
t22 = pkin(2) * t27;
t28 = cos(t32);
t23 = pkin(2) * t28;
t36 = sin(qJ(1));
t29 = t36 * pkin(1);
t47 = t22 + t29;
t26 = pkin(8) + t32;
t19 = sin(t26);
t20 = cos(t26);
t21 = pkin(4) * t34 + pkin(3);
t35 = -pkin(7) - qJ(4);
t46 = t19 * t21 + t20 * t35 + t22;
t45 = t20 * pkin(3) + t19 * qJ(4) + t23;
t44 = -m(3) * pkin(1) - mrSges(2,1);
t43 = -t19 * t35 + t20 * t21 + t23;
t39 = -t28 * mrSges(3,1) + t27 * mrSges(3,2) + t53 * t19 + t52 * t20;
t38 = -t27 * mrSges(3,1) - t28 * mrSges(3,2) + (m(5) * qJ(4) - t53) * t20 + t52 * t19;
t37 = cos(qJ(1));
t30 = t37 * pkin(1);
t14 = t19 * pkin(3);
t1 = [(-mrSges(2,2) * t37 - m(4) * t47 - m(5) * (t14 + t47) - m(6) * (t29 + t46) + t44 * t36 + t38) * g(3) + (t36 * mrSges(2,2) - m(4) * (t23 + t30) - m(5) * (t30 + t45) - m(6) * (t30 + t43) + t44 * t37 + t39) * g(2), (-m(4) * t22 - m(5) * (t14 + t22) - m(6) * t46 + t38) * g(3) + (-m(4) * t23 - m(5) * t45 - m(6) * t43 + t39) * g(2), (-m(4) - t51) * g(1), t51 * (g(2) * t20 + g(3) * t19), -g(1) * t54 + (g(2) * t19 - g(3) * t20) * (mrSges(6,1) * t24 + mrSges(6,2) * t25)];
taug = t1(:);
