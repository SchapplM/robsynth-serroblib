% Calculate potential energy for
% S6RRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRRP4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:11:54
% EndTime: 2019-03-10 01:11:55
% DurationCPUTime: 0.54s
% Computational Cost: add. (260->77), mult. (241->68), div. (0->0), fcn. (215->10), ass. (0->36)
t65 = m(5) * pkin(9) - mrSges(4,2) + mrSges(7,2) + mrSges(5,3) + mrSges(6,3);
t24 = sin(qJ(4));
t27 = cos(qJ(4));
t64 = -m(5) * pkin(3) - t27 * mrSges(5,1) + t24 * mrSges(5,2) - mrSges(4,1);
t23 = qJ(2) + qJ(3);
t17 = sin(t23);
t19 = cos(t23);
t25 = sin(qJ(2));
t28 = cos(qJ(2));
t30 = -pkin(10) - pkin(9);
t59 = -m(6) - m(7);
t62 = -m(3) * pkin(1) - t28 * mrSges(3,1) + t25 * mrSges(3,2) + t64 * t19 - mrSges(2,1) + (-t59 * t30 - t65) * t17;
t61 = -m(2) - m(3);
t60 = -m(4) - m(5);
t57 = -m(1) + t61;
t53 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t52 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t51 = m(3) * pkin(7) + mrSges(5,1) * t24 + mrSges(5,2) * t27 - mrSges(2,2) + mrSges(3,3) + mrSges(4,3);
t50 = pkin(4) * t24;
t26 = sin(qJ(1));
t46 = t19 * t26;
t29 = cos(qJ(1));
t45 = t19 * t29;
t22 = qJ(4) + qJ(5);
t18 = cos(t22);
t44 = t26 * t18;
t21 = pkin(6) + r_base(3);
t14 = pkin(2) * t28 + pkin(1);
t43 = t29 * t14 + r_base(1);
t42 = t25 * pkin(2) + t21;
t31 = -pkin(8) - pkin(7);
t41 = t26 * t14 + t29 * t31 + r_base(2);
t38 = -t26 * t31 + t43;
t16 = sin(t22);
t13 = pkin(4) * t27 + pkin(3);
t1 = (-m(1) * r_base(3) - t25 * mrSges(3,1) - t28 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t60 * t42 + t59 * (t17 * t13 + t19 * t30 + t42) + t61 * t21 + t65 * t19 + (t52 * t16 + t53 * t18 + t64) * t17) * g(3) + (-mrSges(1,2) + t60 * t41 + t59 * (t13 * t46 - t29 * t50 + t41) + t53 * (-t16 * t29 + t19 * t44) + t52 * (t16 * t46 + t18 * t29) + t57 * r_base(2) + t51 * t29 + t62 * t26) * g(2) + (-m(4) * t38 - m(5) * t43 - mrSges(1,1) + t59 * (t13 * t45 + t26 * t50 + t38) + t53 * (t16 * t26 + t18 * t45) + t52 * (t16 * t45 - t44) + t57 * r_base(1) + (m(5) * t31 - t51) * t26 + t62 * t29) * g(1);
U  = t1;
