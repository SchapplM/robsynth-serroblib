% Calculate potential energy for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
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
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPPR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPPR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR1_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR1_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:37:44
% EndTime: 2019-03-09 02:37:44
% DurationCPUTime: 0.42s
% Computational Cost: add. (255->68), mult. (186->53), div. (0->0), fcn. (152->12), ass. (0->29)
t19 = sin(pkin(11));
t20 = cos(pkin(11));
t16 = pkin(11) + qJ(6);
t6 = sin(t16);
t9 = cos(t16);
t50 = -m(6) * pkin(4) - t20 * mrSges(6,1) + t19 * mrSges(6,2) - mrSges(5,1) - m(7) * (t20 * pkin(5) + pkin(4)) - t9 * mrSges(7,1) + t6 * mrSges(7,2);
t49 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3) + m(7) * (-pkin(8) - qJ(5)) - mrSges(7,3);
t48 = -m(1) - m(2);
t47 = -m(3) - m(4);
t46 = m(5) + m(6) + m(7);
t17 = qJ(3) + pkin(10);
t10 = cos(t17);
t23 = sin(qJ(3));
t25 = cos(qJ(3));
t7 = sin(t17);
t43 = -m(4) * pkin(2) - t25 * mrSges(4,1) + t23 * mrSges(4,2) + t50 * t10 + t49 * t7 - mrSges(3,1);
t42 = m(4) * pkin(7) + t6 * mrSges(7,1) + t20 * mrSges(6,2) + t9 * mrSges(7,2) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + (m(7) * pkin(5) + mrSges(6,1)) * t19;
t40 = pkin(6) + r_base(3);
t24 = sin(qJ(1));
t39 = t24 * pkin(1) + r_base(2);
t26 = cos(qJ(1));
t38 = t26 * pkin(1) + r_base(1);
t12 = qJ(2) + t40;
t21 = -qJ(4) - pkin(7);
t18 = qJ(1) + pkin(9);
t11 = cos(t18);
t8 = sin(t18);
t5 = t25 * pkin(3) + pkin(2);
t1 = (-m(1) * r_base(3) - m(2) * t40 - t23 * mrSges(4,1) - t25 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t47 * t12 + t50 * t7 - t46 * (t23 * pkin(3) + t12) - t49 * t10) * g(3) + (-t24 * mrSges(2,1) - t26 * mrSges(2,2) - mrSges(1,2) + t48 * r_base(2) + t47 * t39 - t46 * (t11 * t21 + t8 * t5 + t39) + t43 * t8 + t42 * t11) * g(2) + (-t26 * mrSges(2,1) + t24 * mrSges(2,2) - mrSges(1,1) + t48 * r_base(1) + t47 * t38 - t46 * (t11 * t5 + t38) + (t46 * t21 - t42) * t8 + t43 * t11) * g(1);
U  = t1;
