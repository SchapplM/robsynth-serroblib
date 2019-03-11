% Calculate potential energy for
% S6RRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:01:48
% EndTime: 2019-03-09 22:01:49
% DurationCPUTime: 0.49s
% Computational Cost: add. (238->84), mult. (187->74), div. (0->0), fcn. (149->10), ass. (0->39)
t54 = -mrSges(5,1) + mrSges(6,2);
t53 = mrSges(5,2) - mrSges(6,3);
t52 = -m(2) - m(3);
t51 = m(6) + m(7);
t27 = cos(qJ(1));
t21 = qJ(2) + qJ(3);
t16 = qJ(4) + t21;
t11 = sin(t16);
t40 = qJ(5) * t11;
t12 = cos(t16);
t43 = t27 * t12;
t50 = pkin(4) * t43 + t27 * t40;
t49 = -m(1) - m(4) + t52;
t26 = cos(qJ(2));
t13 = t26 * pkin(2) + pkin(1);
t14 = sin(t21);
t15 = cos(t21);
t23 = sin(qJ(2));
t48 = -m(3) * pkin(1) - m(4) * t13 - t26 * mrSges(3,1) - t15 * mrSges(4,1) + t23 * mrSges(3,2) + t14 * mrSges(4,2) + t53 * t11 + t54 * t12 - mrSges(2,1);
t28 = -pkin(8) - pkin(7);
t47 = -m(3) * pkin(7) + m(4) * t28 - mrSges(6,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t24 = sin(qJ(1));
t46 = t24 * t12;
t22 = sin(qJ(6));
t45 = t24 * t22;
t25 = cos(qJ(6));
t44 = t24 * t25;
t42 = t27 * t22;
t41 = t27 * t25;
t19 = pkin(6) + r_base(3);
t3 = pkin(3) * t15 + t13;
t39 = t27 * t3 + r_base(1);
t20 = -pkin(9) + t28;
t38 = t27 * t20 + t24 * t3 + r_base(2);
t37 = t23 * pkin(2) + t19;
t36 = pkin(3) * t14 + t37;
t35 = pkin(4) * t46 + t24 * t40 + t38;
t32 = -t24 * t20 + t39;
t1 = (-m(1) * r_base(3) - m(4) * t37 - m(5) * t36 - t23 * mrSges(3,1) - t14 * mrSges(4,1) - t26 * mrSges(3,2) - t15 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - t51 * (t11 * pkin(4) + t36) + t52 * t19 + (t22 * mrSges(7,1) + t25 * mrSges(7,2) + t51 * qJ(5) - t53) * t12 + (-m(7) * pkin(10) - mrSges(7,3) + t54) * t11) * g(3) + (-mrSges(1,2) - m(5) * t38 - m(6) * t35 - m(7) * (pkin(10) * t46 + t35) - (t11 * t45 - t41) * mrSges(7,1) - (t11 * t44 + t42) * mrSges(7,2) - mrSges(7,3) * t46 + t49 * r_base(2) + (m(7) * pkin(5) - t47) * t27 + t48 * t24) * g(2) + (-mrSges(1,1) - m(5) * t32 - m(6) * (t32 + t50) - m(7) * (pkin(10) * t43 + t39 + t50) - (t11 * t42 + t44) * mrSges(7,1) - (t11 * t41 - t45) * mrSges(7,2) - mrSges(7,3) * t43 + t49 * r_base(1) + t48 * t27 + (-m(7) * (pkin(5) - t20) + t47) * t24) * g(1);
U  = t1;
