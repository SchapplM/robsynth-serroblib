% Calculate potential energy for
% S6RRPRRP11
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
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP11_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRP11_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP11_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP11_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:45:18
% EndTime: 2019-03-09 12:45:18
% DurationCPUTime: 0.61s
% Computational Cost: add. (197->84), mult. (251->78), div. (0->0), fcn. (221->8), ass. (0->36)
t21 = qJ(4) + qJ(5);
t12 = sin(t21);
t22 = sin(qJ(4));
t46 = pkin(4) * t22;
t58 = -m(6) * t46 - m(7) * (pkin(5) * t12 + t46) + mrSges(3,2) - mrSges(4,3);
t28 = -pkin(9) - pkin(8);
t57 = m(6) * t28 + m(7) * (-qJ(6) + t28) - mrSges(3,1) + mrSges(4,2) - mrSges(6,3) - mrSges(7,3);
t56 = -m(1) - m(2);
t24 = sin(qJ(1));
t23 = sin(qJ(2));
t39 = qJ(3) * t23;
t26 = cos(qJ(2));
t42 = t24 * t26;
t55 = pkin(2) * t42 + t24 * t39;
t54 = mrSges(5,2) * t23;
t53 = mrSges(6,1) + mrSges(7,1);
t52 = mrSges(6,2) + mrSges(7,2);
t51 = -m(4) - m(6) - m(7);
t50 = m(5) - t51;
t49 = -m(5) * pkin(8) - mrSges(5,3);
t48 = t58 * t23 + t57 * t26 - mrSges(2,1);
t25 = cos(qJ(4));
t11 = t25 * pkin(4) + pkin(3);
t13 = cos(t21);
t47 = m(6) * t11 + m(7) * (pkin(5) * t13 + t11) - t22 * mrSges(5,2) + mrSges(4,1) - mrSges(2,2) + mrSges(3,3);
t45 = t23 * t24;
t27 = cos(qJ(1));
t44 = t23 * t27;
t43 = t24 * t25;
t41 = t25 * t27;
t40 = t26 * t27;
t20 = pkin(6) + r_base(3);
t38 = t24 * pkin(1) + r_base(2);
t36 = t27 * pkin(1) + t24 * pkin(7) + r_base(1);
t33 = -t27 * pkin(7) + t38;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t20 - t50 * (t23 * pkin(2) + t20) + (t22 * mrSges(5,1) + t25 * mrSges(5,2) + t50 * qJ(3) + t53 * t12 + t52 * t13 - t58) * t26 + (t49 + t57) * t23) * g(3) + (-mrSges(1,2) - m(3) * t33 - m(5) * (pkin(8) * t42 + t38 + t55) - (t22 * t45 - t41) * mrSges(5,1) - t43 * t54 - mrSges(5,3) * t42 + t56 * r_base(2) - t53 * (t12 * t45 - t13 * t27) - t52 * (t12 * t27 + t13 * t45) + t51 * (t33 + t55) + (-m(5) * (-pkin(3) - pkin(7)) + t47) * t27 + t48 * t24) * g(2) + (-mrSges(1,1) - m(3) * t36 - (t22 * t44 + t43) * mrSges(5,1) - t41 * t54 + t56 * r_base(1) + t49 * t40 - t53 * (t12 * t44 + t13 * t24) - t52 * (-t12 * t24 + t13 * t44) - t50 * (pkin(2) * t40 + t27 * t39 + t36) + (-m(5) * pkin(3) - t47) * t24 + t48 * t27) * g(1);
U  = t1;
