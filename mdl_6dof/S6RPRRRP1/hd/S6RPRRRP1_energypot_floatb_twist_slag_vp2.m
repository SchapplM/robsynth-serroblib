% Calculate potential energy for
% S6RPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP1_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:55:30
% EndTime: 2019-03-09 05:55:31
% DurationCPUTime: 0.46s
% Computational Cost: add. (259->65), mult. (199->50), div. (0->0), fcn. (169->10), ass. (0->30)
t25 = sin(qJ(5));
t28 = cos(qJ(5));
t51 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t52 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t62 = t25 * t51 + t28 * t52 - mrSges(5,1);
t59 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t53 = -m(6) - m(7);
t58 = -m(5) + t53;
t57 = -m(4) * pkin(7) + t52 * t25 - t51 * t28 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t24 = qJ(3) + qJ(4);
t18 = sin(t24);
t19 = cos(t24);
t26 = sin(qJ(3));
t29 = cos(qJ(3));
t56 = -m(4) * pkin(2) - t29 * mrSges(4,1) + t26 * mrSges(4,2) - mrSges(3,1) + (t53 * pkin(9) + t59) * t18 + (t53 * pkin(4) + t62) * t19;
t55 = -m(1) - m(2);
t54 = -m(3) - m(4);
t42 = pkin(6) + r_base(3);
t27 = sin(qJ(1));
t41 = t27 * pkin(1) + r_base(2);
t30 = cos(qJ(1));
t40 = t30 * pkin(1) + r_base(1);
t17 = qJ(2) + t42;
t38 = t26 * pkin(3) + t17;
t31 = -pkin(8) - pkin(7);
t23 = qJ(1) + pkin(10);
t16 = cos(t23);
t15 = sin(t23);
t14 = t29 * pkin(3) + pkin(2);
t1 = (-m(1) * r_base(3) - m(2) * t42 - m(5) * t38 - t26 * mrSges(4,1) - t29 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t53 * (t18 * pkin(4) - t19 * pkin(9) + t38) + t54 * t17 - t59 * t19 + t62 * t18) * g(3) + (-t27 * mrSges(2,1) - t30 * mrSges(2,2) + t54 * t41 + t55 * r_base(2) - mrSges(1,2) + t58 * (t15 * t14 + t16 * t31 + t41) - t57 * t16 + t56 * t15) * g(2) + (-t30 * mrSges(2,1) + t27 * mrSges(2,2) + t54 * t40 + t55 * r_base(1) - mrSges(1,1) + t58 * (t16 * t14 - t15 * t31 + t40) + t57 * t15 + t56 * t16) * g(1);
U  = t1;
