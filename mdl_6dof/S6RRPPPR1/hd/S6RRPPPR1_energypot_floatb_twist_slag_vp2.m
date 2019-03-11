% Calculate potential energy for
% S6RRPPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPPR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPPR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR1_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR1_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:05:58
% EndTime: 2019-03-09 08:05:59
% DurationCPUTime: 0.62s
% Computational Cost: add. (260->78), mult. (288->71), div. (0->0), fcn. (278->10), ass. (0->40)
t64 = -mrSges(4,2) + mrSges(6,2) + mrSges(5,3) - mrSges(7,3);
t22 = qJ(2) + pkin(9);
t19 = sin(t22);
t20 = cos(t22);
t63 = pkin(3) * t20 + qJ(4) * t19;
t28 = sin(qJ(2));
t31 = cos(qJ(2));
t62 = -m(3) * pkin(1) - t31 * mrSges(3,1) - t20 * mrSges(4,1) + t28 * mrSges(3,2) - mrSges(2,1) + (m(7) * pkin(8) - t64) * t19;
t61 = -m(2) - m(3);
t60 = -m(6) - m(7);
t24 = sin(pkin(10));
t25 = cos(pkin(10));
t59 = (pkin(4) * t25 + qJ(5) * t24) * t19;
t58 = -m(1) + t61;
t55 = m(3) * pkin(7) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3);
t27 = sin(qJ(6));
t30 = cos(qJ(6));
t53 = -t27 * mrSges(7,1) - t30 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t52 = -m(7) * pkin(5) - mrSges(7,1) * t30 + mrSges(7,2) * t27 - mrSges(5,1) - mrSges(6,1);
t32 = cos(qJ(1));
t48 = t24 * t32;
t47 = t25 * t32;
t29 = sin(qJ(1));
t46 = t29 * t24;
t45 = t29 * t25;
t23 = pkin(6) + r_base(3);
t43 = pkin(2) * t28 + t23;
t18 = pkin(2) * t31 + pkin(1);
t26 = -qJ(3) - pkin(7);
t42 = t18 * t29 + t26 * t32 + r_base(2);
t41 = pkin(3) * t19 + t43;
t39 = t18 * t32 - t26 * t29 + r_base(1);
t38 = t29 * t63 + t42;
t37 = -qJ(4) * t20 + t41;
t36 = t32 * t63 + t39;
t6 = t20 * t47 + t46;
t5 = t20 * t48 - t45;
t4 = t20 * t45 - t48;
t3 = t20 * t46 + t47;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - t28 * mrSges(3,1) - t31 * mrSges(3,2) - m(4) * t43 - m(5) * t37 - m(6) * (t37 + t59) - m(7) * (t41 + t59) + t61 * t23 + (-m(7) * (pkin(8) - qJ(4)) + t64) * t20 + (t24 * t53 + t25 * t52 - mrSges(4,1)) * t19) * g(3) + (-m(4) * t42 - m(5) * t38 - mrSges(1,2) + t60 * (pkin(4) * t4 + qJ(5) * t3 + t38) + t52 * t4 + t53 * t3 + t58 * r_base(2) + t55 * t32 + t62 * t29) * g(2) + (-m(4) * t39 - m(5) * t36 - mrSges(1,1) + t60 * (pkin(4) * t6 + qJ(5) * t5 + t36) + t52 * t6 + t53 * t5 + t58 * r_base(1) - t55 * t29 + t62 * t32) * g(1);
U  = t1;
