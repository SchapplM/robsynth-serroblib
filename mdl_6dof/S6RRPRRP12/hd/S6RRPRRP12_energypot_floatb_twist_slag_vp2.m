% Calculate potential energy for
% S6RRPRRP12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP12_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRP12_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP12_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP12_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP12_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:50:15
% EndTime: 2019-03-09 12:50:16
% DurationCPUTime: 0.65s
% Computational Cost: add. (204->89), mult. (266->83), div. (0->0), fcn. (240->8), ass. (0->41)
t62 = -mrSges(3,1) + mrSges(4,2);
t61 = mrSges(3,2) - mrSges(4,3);
t60 = -m(1) - m(2);
t59 = m(4) + m(5);
t58 = -m(6) - m(7);
t24 = sin(qJ(1));
t23 = sin(qJ(2));
t43 = qJ(3) * t23;
t26 = cos(qJ(2));
t47 = t24 * t26;
t57 = pkin(2) * t47 + t24 * t43;
t56 = mrSges(5,2) * t23;
t55 = t61 * t23 + t62 * t26 - mrSges(2,1);
t54 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t53 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t52 = -m(5) * pkin(8) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t22 = sin(qJ(4));
t51 = -t22 * mrSges(5,2) + mrSges(4,1) - mrSges(2,2) + mrSges(3,3);
t50 = t23 * t24;
t27 = cos(qJ(1));
t49 = t23 * t27;
t25 = cos(qJ(4));
t48 = t24 * t25;
t46 = t25 * t27;
t45 = t26 * t27;
t28 = -pkin(9) - pkin(8);
t44 = t26 * t28;
t20 = pkin(6) + r_base(3);
t42 = t24 * pkin(1) + r_base(2);
t41 = t22 * t50;
t40 = t22 * t49;
t39 = t23 * pkin(2) + t20;
t37 = t27 * pkin(1) + t24 * pkin(7) + r_base(1);
t36 = t42 + t57;
t33 = -t27 * pkin(7) + t42;
t32 = pkin(2) * t45 + t27 * t43 + t37;
t21 = qJ(4) + qJ(5);
t15 = cos(t21);
t14 = sin(t21);
t13 = pkin(4) * t25 + pkin(3);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - t59 * t39 + t58 * (-t23 * t28 + t39) + (-m(2) - m(3)) * t20 + (t22 * mrSges(5,1) + t25 * mrSges(5,2) + t58 * (-pkin(4) * t22 - qJ(3)) - t53 * t15 + t54 * t14 + t59 * qJ(3) - t61) * t26 + (t52 + t62) * t23) * g(3) + (-mrSges(1,2) - m(3) * t33 - m(4) * (t33 + t57) - m(5) * t36 - (t41 - t46) * mrSges(5,1) - t48 * t56 + t60 * r_base(2) + t58 * (-t24 * t44 + pkin(4) * t41 + (-pkin(7) - t13) * t27 + t36) - t54 * (t14 * t50 - t15 * t27) + t53 * (t14 * t27 + t15 * t50) + t52 * t47 + (-m(5) * (-pkin(3) - pkin(7)) + t51) * t27 + t55 * t24) * g(2) + (-mrSges(1,1) - m(3) * t37 - (t40 + t48) * mrSges(5,1) - t46 * t56 + t60 * r_base(1) - t59 * t32 + t58 * (pkin(4) * t40 + t24 * t13 - t27 * t44 + t32) - t54 * (t14 * t49 + t15 * t24) - t53 * (t14 * t24 - t15 * t49) + t52 * t45 + t55 * t27 + (-m(5) * pkin(3) - t51) * t24) * g(1);
U  = t1;
