% Calculate potential energy for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
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
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPP4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRPP4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:58:44
% EndTime: 2019-03-09 09:58:45
% DurationCPUTime: 0.65s
% Computational Cost: add. (204->89), mult. (266->83), div. (0->0), fcn. (240->8), ass. (0->41)
t62 = -mrSges(3,1) + mrSges(4,2);
t61 = mrSges(3,2) - mrSges(4,3);
t60 = -m(1) - m(2);
t59 = m(4) + m(5);
t58 = -m(6) - m(7);
t25 = sin(qJ(1));
t24 = sin(qJ(2));
t43 = qJ(3) * t24;
t27 = cos(qJ(2));
t46 = t25 * t27;
t57 = pkin(2) * t46 + t25 * t43;
t56 = t61 * t24 + t62 * t27 - mrSges(2,1);
t55 = mrSges(4,1) + mrSges(3,3) - mrSges(2,2);
t54 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t53 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t52 = -m(5) * pkin(8) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t20 = qJ(4) + pkin(9);
t14 = sin(t20);
t51 = t14 * t25;
t28 = cos(qJ(1));
t50 = t24 * t28;
t15 = cos(t20);
t49 = t25 * t15;
t23 = sin(qJ(4));
t48 = t25 * t23;
t26 = cos(qJ(4));
t47 = t25 * t26;
t45 = t26 * t28;
t44 = t27 * t28;
t21 = pkin(6) + r_base(3);
t42 = t25 * pkin(1) + r_base(2);
t41 = t23 * t50;
t40 = t24 * t48;
t39 = t24 * pkin(2) + t21;
t37 = t28 * pkin(1) + t25 * pkin(7) + r_base(1);
t36 = t42 + t57;
t33 = -pkin(7) * t28 + t42;
t32 = pkin(2) * t44 + t28 * t43 + t37;
t22 = -qJ(5) - pkin(8);
t13 = pkin(4) * t26 + pkin(3);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - t59 * t39 + t58 * (-t22 * t24 + t39) + (-m(2) - m(3)) * t21 + (t23 * mrSges(5,1) + t26 * mrSges(5,2) + t58 * (-pkin(4) * t23 - qJ(3)) - t53 * t15 + t54 * t14 + t59 * qJ(3) - t61) * t27 + (t52 + t62) * t24) * g(3) + (-mrSges(1,2) - m(3) * t33 - m(4) * (t33 + t57) - m(5) * t36 - (t40 - t45) * mrSges(5,1) - t24 * t47 * mrSges(5,2) + t60 * r_base(2) + t58 * (-t22 * t46 + pkin(4) * t40 + (-pkin(7) - t13) * t28 + t36) - t54 * (-t15 * t28 + t24 * t51) + t53 * (t14 * t28 + t24 * t49) + t52 * t46 + (-m(5) * (-pkin(3) - pkin(7)) - t23 * mrSges(5,2) + t55) * t28 + t56 * t25) * g(2) + (-mrSges(1,1) - m(3) * t37 - (t41 + t47) * mrSges(5,1) - (t24 * t45 - t48) * mrSges(5,2) + t60 * r_base(1) - t59 * t32 + t58 * (pkin(4) * t41 + t25 * t13 - t22 * t44 + t32) - t54 * (t14 * t50 + t49) - t53 * (-t15 * t50 + t51) + t52 * t44 + t56 * t28 + (-m(5) * pkin(3) - t55) * t25) * g(1);
U  = t1;
