% Calculate potential energy for
% S6RPRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPPR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPPR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:44:06
% EndTime: 2019-03-09 02:44:06
% DurationCPUTime: 0.42s
% Computational Cost: add. (213->76), mult. (191->60), div. (0->0), fcn. (153->8), ass. (0->33)
t53 = -m(7) * pkin(8) - mrSges(4,1) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t19 = sin(qJ(6));
t22 = cos(qJ(6));
t52 = t22 * mrSges(7,1) - t19 * mrSges(7,2) + mrSges(6,1) - mrSges(4,2) + mrSges(5,3);
t51 = -m(1) - m(2);
t50 = m(6) + m(7);
t18 = qJ(1) + pkin(9);
t11 = sin(t18);
t20 = sin(qJ(3));
t42 = qJ(4) * t20;
t23 = cos(qJ(3));
t44 = t11 * t23;
t49 = pkin(3) * t44 + t11 * t42;
t12 = cos(t18);
t48 = pkin(4) * t44 + t12 * qJ(5);
t46 = -mrSges(3,1) + t53 * t23 + (-m(7) * pkin(5) - t52) * t20;
t45 = t19 * mrSges(7,1) + t22 * mrSges(7,2) + mrSges(3,2) - mrSges(5,2) - mrSges(4,3) + mrSges(6,3);
t43 = t12 * t23;
t41 = pkin(6) + r_base(3);
t21 = sin(qJ(1));
t40 = t21 * pkin(1) + r_base(2);
t24 = cos(qJ(1));
t39 = t24 * pkin(1) + r_base(1);
t38 = t11 * pkin(2) + t40;
t13 = qJ(2) + t41;
t37 = t12 * pkin(2) + t11 * pkin(7) + t39;
t36 = t20 * pkin(3) + t13;
t30 = -t12 * pkin(7) + t38;
t29 = pkin(3) * t43 + t12 * t42 + t37;
t27 = -t23 * qJ(4) + t36;
t26 = t30 + t49;
t14 = t20 * pkin(4);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t41 - mrSges(2,3) - mrSges(3,3) - m(5) * t27 - m(6) * (t14 + t27) - m(7) * (t14 + t36) + (-m(3) - m(4)) * t13 + (-m(7) * (-pkin(5) - qJ(4)) + t52) * t23 + t53 * t20) * g(3) + (-mrSges(1,2) - t21 * mrSges(2,1) - t24 * mrSges(2,2) - m(3) * t40 - m(4) * t30 - m(5) * t26 - m(6) * (t26 + t48) - m(7) * (t38 + t48 + t49) + t51 * r_base(2) + (m(7) * pkin(7) - t45) * t12 + t46 * t11) * g(2) + (-m(3) * t39 - m(4) * t37 - m(5) * t29 - t24 * mrSges(2,1) + t21 * mrSges(2,2) - mrSges(1,1) + t51 * r_base(1) - t50 * (pkin(4) * t43 + t29) + t46 * t12 + (t50 * qJ(5) + t45) * t11) * g(1);
U  = t1;
