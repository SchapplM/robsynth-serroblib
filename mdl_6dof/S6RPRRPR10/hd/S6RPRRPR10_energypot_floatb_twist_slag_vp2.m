% Calculate potential energy for
% S6RPRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR10_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPR10_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR10_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR10_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR10_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:34:17
% EndTime: 2019-03-09 05:34:17
% DurationCPUTime: 0.60s
% Computational Cost: add. (174->80), mult. (276->71), div. (0->0), fcn. (266->8), ass. (0->38)
t58 = -mrSges(4,2) + mrSges(6,2) + mrSges(5,3);
t55 = -m(6) - m(7);
t57 = -m(5) + t55;
t56 = -m(1) - m(2);
t54 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t53 = m(7) * pkin(9) + mrSges(7,3);
t22 = sin(qJ(3));
t26 = cos(qJ(3));
t52 = -t22 * mrSges(4,1) + t58 * t26 + mrSges(2,2) - mrSges(3,3);
t20 = sin(qJ(6));
t24 = cos(qJ(6));
t51 = t20 * mrSges(7,1) + t24 * mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t50 = -m(7) * pkin(5) - t24 * mrSges(7,1) + t20 * mrSges(7,2) - mrSges(5,1) - mrSges(6,1);
t49 = pkin(3) * t22;
t21 = sin(qJ(4));
t27 = cos(qJ(1));
t48 = t21 * t27;
t23 = sin(qJ(1));
t47 = t23 * t21;
t25 = cos(qJ(4));
t46 = t23 * t25;
t45 = t23 * t26;
t42 = t27 * t25;
t19 = pkin(6) + r_base(3);
t41 = pkin(8) * t45;
t40 = t23 * pkin(1) + r_base(2);
t39 = pkin(2) + t19;
t37 = t23 * pkin(7) + t40;
t36 = t27 * pkin(1) + t23 * qJ(2) + r_base(1);
t35 = t27 * t26 * pkin(8) + t37;
t34 = t27 * pkin(7) + t36;
t31 = t23 * t49 + t34;
t3 = t22 * t47 - t42;
t4 = t22 * t46 + t48;
t28 = t4 * pkin(4) + qJ(5) * t3 + t31;
t6 = -t22 * t42 + t47;
t5 = t22 * t48 + t46;
t1 = (-m(1) * r_base(3) - m(4) * t39 - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t19 + (t53 - t58) * t22 + t57 * (t26 * pkin(3) + t22 * pkin(8) + t39) + (t55 * (pkin(4) * t25 + qJ(5) * t21) - t51 * t21 + t50 * t25 - mrSges(4,1)) * t26) * g(3) + (-m(3) * t40 - m(4) * t37 - m(5) * t35 - mrSges(1,2) + t56 * r_base(2) + t55 * (t6 * pkin(4) - t5 * qJ(5) + t35) + t50 * t6 + t51 * t5 + t54 * t23 + (t53 * t26 + t57 * (-qJ(2) - t49) + (m(3) + m(4)) * qJ(2) - t52) * t27) * g(2) + (-mrSges(1,1) - m(3) * t36 - m(4) * t34 - m(5) * (t31 - t41) - m(6) * (t28 - t41) - m(7) * t28 - (m(7) * (-pkin(8) + pkin(9)) + mrSges(7,3)) * t45 + t56 * r_base(1) + t50 * t4 - t51 * t3 + t54 * t27 + t52 * t23) * g(1);
U  = t1;
