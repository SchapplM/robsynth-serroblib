% Calculate potential energy for
% S6RRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPPR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:28:18
% EndTime: 2019-03-09 15:28:18
% DurationCPUTime: 0.43s
% Computational Cost: add. (210->74), mult. (205->58), div. (0->0), fcn. (167->8), ass. (0->32)
t54 = -m(7) * pkin(9) - mrSges(4,1) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t19 = sin(qJ(6));
t22 = cos(qJ(6));
t53 = t22 * mrSges(7,1) - t19 * mrSges(7,2) + mrSges(6,1) - mrSges(4,2) + mrSges(5,3);
t52 = -m(2) - m(3);
t51 = -m(6) - m(7);
t24 = cos(qJ(1));
t18 = qJ(2) + qJ(3);
t13 = sin(t18);
t43 = qJ(4) * t13;
t14 = cos(t18);
t44 = t14 * t24;
t50 = pkin(3) * t44 + t24 * t43;
t49 = -m(1) + t52;
t20 = sin(qJ(2));
t23 = cos(qJ(2));
t47 = -m(3) * pkin(1) - t23 * mrSges(3,1) + t20 * mrSges(3,2) - mrSges(2,1) + t54 * t14 + (-m(7) * pkin(5) - t53) * t13;
t46 = m(3) * pkin(7) - t19 * mrSges(7,1) - t22 * mrSges(7,2) - mrSges(2,2) + mrSges(5,2) + mrSges(3,3) + mrSges(4,3) - mrSges(6,3);
t21 = sin(qJ(1));
t45 = t14 * t21;
t17 = pkin(6) + r_base(3);
t11 = t23 * pkin(2) + pkin(1);
t41 = t24 * t11 + r_base(1);
t40 = t20 * pkin(2) + t17;
t25 = -pkin(8) - pkin(7);
t39 = t21 * t11 + t24 * t25 + r_base(2);
t38 = t13 * pkin(3) + t40;
t36 = pkin(3) * t45 + t21 * t43 + t39;
t30 = -t21 * t25 + t41;
t28 = -t14 * qJ(4) + t38;
t9 = t13 * pkin(4);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - t20 * mrSges(3,1) - t23 * mrSges(3,2) - m(4) * t40 - m(5) * t28 - m(6) * (t28 + t9) - m(7) * (t38 + t9) + t52 * t17 + (-m(7) * (-pkin(5) - qJ(4)) + t53) * t14 + t54 * t13) * g(3) + (-m(4) * t39 - m(5) * t36 - mrSges(1,2) + t51 * (pkin(4) * t45 + t24 * qJ(5) + t36) + t49 * r_base(2) + t46 * t24 + t47 * t21) * g(2) + (-mrSges(1,1) - m(4) * t30 - m(5) * (t30 + t50) + t51 * (pkin(4) * t44 + t41 + t50) + t49 * r_base(1) + t47 * t24 + (t51 * (-qJ(5) - t25) - t46) * t21) * g(1);
U  = t1;
