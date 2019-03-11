% Calculate potential energy for
% S6RPPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRP2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:00:11
% EndTime: 2019-03-09 02:00:11
% DurationCPUTime: 0.41s
% Computational Cost: add. (259->77), mult. (199->69), div. (0->0), fcn. (169->10), ass. (0->37)
t57 = -m(1) - m(2);
t56 = -m(3) - m(4);
t55 = -m(6) - m(7);
t54 = -mrSges(6,3) - mrSges(7,2);
t23 = pkin(10) + qJ(4);
t15 = sin(t23);
t17 = cos(t23);
t25 = sin(pkin(10));
t26 = cos(pkin(10));
t53 = -m(4) * pkin(2) - t26 * mrSges(4,1) - t17 * mrSges(5,1) + t25 * mrSges(4,2) + t15 * mrSges(5,2) - mrSges(3,1);
t52 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t51 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t50 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t49 = pkin(4) * t17;
t24 = qJ(1) + pkin(9);
t16 = sin(t24);
t48 = t16 * t15;
t28 = sin(qJ(5));
t47 = t16 * t28;
t30 = cos(qJ(5));
t46 = t16 * t30;
t18 = cos(t24);
t45 = t18 * t15;
t44 = t18 * t28;
t43 = t18 * t30;
t42 = pkin(6) + r_base(3);
t29 = sin(qJ(1));
t41 = t29 * pkin(1) + r_base(2);
t31 = cos(qJ(1));
t40 = t31 * pkin(1) + r_base(1);
t19 = qJ(2) + t42;
t14 = t26 * pkin(3) + pkin(2);
t27 = -pkin(7) - qJ(3);
t39 = t16 * t14 + t18 * t27 + t41;
t38 = t25 * pkin(3) + t19;
t35 = t18 * t14 - t16 * t27 + t40;
t1 = (-m(1) * r_base(3) - m(2) * t42 - m(5) * t38 - t25 * mrSges(4,1) - t26 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t55 * (t15 * pkin(4) - t17 * pkin(8) + t38) + t56 * t19 + (-mrSges(5,2) - t54) * t17 + (t28 * t51 + t30 * t52 - mrSges(5,1)) * t15) * g(3) + (-m(5) * t39 - t29 * mrSges(2,1) - t31 * mrSges(2,2) - mrSges(1,2) + t57 * r_base(2) + t54 * t48 + t56 * t41 + t55 * (pkin(8) * t48 + t16 * t49 + t39) + t52 * (t17 * t46 - t44) + t51 * (t17 * t47 + t43) + t50 * t18 + t53 * t16) * g(2) + (-m(5) * t35 - t31 * mrSges(2,1) + t29 * mrSges(2,2) - mrSges(1,1) + t57 * r_base(1) + t54 * t45 + t56 * t40 + t55 * (pkin(8) * t45 + t18 * t49 + t35) + t52 * (t17 * t43 + t47) + t51 * (t17 * t44 - t46) + t53 * t18 - t50 * t16) * g(1);
U  = t1;
