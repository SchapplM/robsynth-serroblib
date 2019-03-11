% Calculate potential energy for
% S6RRPPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR8_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRR8_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR8_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR8_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR8_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:22:45
% EndTime: 2019-03-09 09:22:46
% DurationCPUTime: 0.75s
% Computational Cost: add. (227->85), mult. (376->79), div. (0->0), fcn. (388->10), ass. (0->38)
t63 = mrSges(5,2) + mrSges(4,3) - mrSges(6,3) - mrSges(7,3) - mrSges(3,2);
t26 = sin(qJ(2));
t29 = cos(qJ(2));
t31 = -pkin(9) - pkin(8);
t62 = -t29 * mrSges(3,1) - mrSges(2,1) + (m(6) * pkin(8) - m(7) * t31 - t63) * t26;
t61 = -m(1) - m(2);
t60 = -m(5) - m(6);
t23 = sin(pkin(10));
t24 = cos(pkin(10));
t59 = (pkin(3) * t24 + qJ(4) * t23) * t26;
t58 = mrSges(2,2) - mrSges(3,3);
t22 = qJ(5) + qJ(6);
t15 = sin(t22);
t16 = cos(t22);
t28 = cos(qJ(5));
t54 = -t15 * mrSges(7,1) - mrSges(6,2) * t28 - t16 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t25 = sin(qJ(5));
t53 = -m(7) * (pkin(5) * t25 + qJ(4)) - mrSges(6,1) * t25 + t54;
t52 = -m(6) * pkin(4) - m(7) * (pkin(5) * t28 + pkin(4)) - t28 * mrSges(6,1) - t16 * mrSges(7,1) + t25 * mrSges(6,2) + t15 * mrSges(7,2) - mrSges(4,1) - mrSges(5,1);
t27 = sin(qJ(1));
t47 = t27 * t29;
t30 = cos(qJ(1));
t46 = t29 * t30;
t45 = qJ(3) * t26;
t21 = pkin(6) + r_base(3);
t44 = t26 * pkin(2) + t21;
t42 = t30 * pkin(1) + t27 * pkin(7) + r_base(1);
t39 = t27 * pkin(1) - t30 * pkin(7) + r_base(2);
t38 = pkin(2) * t46 + t30 * t45 + t42;
t37 = -t29 * qJ(3) + t44;
t6 = t23 * t27 + t24 * t46;
t36 = t6 * pkin(3) + t38;
t35 = pkin(2) * t47 + t27 * t45 + t39;
t4 = -t23 * t30 + t24 * t47;
t34 = t4 * pkin(3) + t35;
t5 = t23 * t46 - t27 * t24;
t3 = t23 * t47 + t24 * t30;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - m(4) * t37 - m(5) * (t37 + t59) + (-m(6) - m(7)) * (t44 + t59) + (-m(2) - m(3)) * t21 + (-m(6) * (pkin(8) - qJ(3)) - m(7) * (-qJ(3) - t31) + t63) * t29 + (t52 * t24 - mrSges(3,1) + ((-m(7) * pkin(5) - mrSges(6,1)) * t25 + t54) * t23) * t26) * g(3) + (-m(3) * t39 - m(4) * t35 - m(7) * t34 - mrSges(1,2) + t61 * r_base(2) + t60 * (t3 * qJ(4) + t34) - t58 * t30 + t52 * t4 + t53 * t3 + t62 * t27) * g(2) + (-m(3) * t42 - m(4) * t38 - m(7) * t36 - mrSges(1,1) + t61 * r_base(1) + t60 * (t5 * qJ(4) + t36) + t52 * t6 + t53 * t5 + t58 * t27 + t62 * t30) * g(1);
U  = t1;
