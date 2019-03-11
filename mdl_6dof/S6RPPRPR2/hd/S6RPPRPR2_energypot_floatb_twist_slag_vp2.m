% Calculate potential energy for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRPR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRPR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:41:32
% EndTime: 2019-03-09 01:41:33
% DurationCPUTime: 0.47s
% Computational Cost: add. (235->84), mult. (173->74), div. (0->0), fcn. (135->10), ass. (0->39)
t53 = -mrSges(5,1) + mrSges(6,2);
t52 = mrSges(5,2) - mrSges(6,3);
t51 = -m(1) - m(2);
t50 = -m(3) - m(4);
t49 = m(6) + m(7);
t19 = qJ(1) + pkin(9);
t13 = cos(t19);
t18 = pkin(10) + qJ(4);
t10 = sin(t18);
t39 = qJ(5) * t10;
t12 = cos(t18);
t42 = t13 * t12;
t48 = pkin(4) * t42 + t13 * t39;
t20 = sin(pkin(10));
t21 = cos(pkin(10));
t47 = -m(4) * pkin(2) - t21 * mrSges(4,1) + t20 * mrSges(4,2) + t52 * t10 + t53 * t12 - mrSges(3,1);
t46 = -m(4) * qJ(3) - mrSges(6,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t11 = sin(t19);
t45 = t11 * t12;
t23 = sin(qJ(6));
t44 = t11 * t23;
t25 = cos(qJ(6));
t43 = t11 * t25;
t41 = t13 * t23;
t40 = t13 * t25;
t38 = pkin(6) + r_base(3);
t24 = sin(qJ(1));
t37 = t24 * pkin(1) + r_base(2);
t26 = cos(qJ(1));
t36 = t26 * pkin(1) + r_base(1);
t9 = t21 * pkin(3) + pkin(2);
t35 = t13 * t9 + t36;
t14 = qJ(2) + t38;
t22 = -pkin(7) - qJ(3);
t34 = t11 * t9 + t13 * t22 + t37;
t33 = t20 * pkin(3) + t14;
t29 = pkin(4) * t45 + t11 * t39 + t34;
t28 = -t11 * t22 + t35;
t1 = (-m(1) * r_base(3) - m(2) * t38 - m(5) * t33 - t20 * mrSges(4,1) - t21 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - t49 * (t10 * pkin(4) + t33) + t50 * t14 + (t23 * mrSges(7,1) + t25 * mrSges(7,2) + t49 * qJ(5) - t52) * t12 + (-m(7) * pkin(8) - mrSges(7,3) + t53) * t10) * g(3) + (-mrSges(1,2) - t24 * mrSges(2,1) - t26 * mrSges(2,2) - m(5) * t34 - m(6) * t29 - m(7) * (pkin(8) * t45 + t29) - (t10 * t44 - t40) * mrSges(7,1) - (t10 * t43 + t41) * mrSges(7,2) - mrSges(7,3) * t45 + t51 * r_base(2) + t50 * t37 + (m(7) * pkin(5) - t46) * t13 + t47 * t11) * g(2) + (-mrSges(1,1) - t26 * mrSges(2,1) + t24 * mrSges(2,2) - m(5) * t28 - m(6) * (t28 + t48) - m(7) * (pkin(8) * t42 + t35 + t48) - (t10 * t41 + t43) * mrSges(7,1) - (t10 * t40 - t44) * mrSges(7,2) - mrSges(7,3) * t42 + t51 * r_base(1) + t50 * t36 + t47 * t13 + (-m(7) * (pkin(5) - t22) + t46) * t11) * g(1);
U  = t1;
