% Calculate potential energy for
% S6RPRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:36:42
% EndTime: 2019-03-09 03:36:43
% DurationCPUTime: 0.42s
% Computational Cost: add. (255->68), mult. (186->53), div. (0->0), fcn. (152->12), ass. (0->29)
t18 = qJ(5) + qJ(6);
t11 = sin(t18);
t12 = cos(t18);
t20 = sin(qJ(5));
t23 = cos(qJ(5));
t50 = -mrSges(5,1) - m(6) * pkin(4) - t23 * mrSges(6,1) + t20 * mrSges(6,2) - m(7) * (t23 * pkin(5) + pkin(4)) - t12 * mrSges(7,1) + t11 * mrSges(7,2);
t49 = mrSges(5,2) + m(7) * (-pkin(9) - pkin(8)) - mrSges(7,3) - m(6) * pkin(8) - mrSges(6,3);
t48 = -m(1) - m(2);
t47 = -m(3) - m(4);
t46 = m(5) + m(6) + m(7);
t21 = sin(qJ(3));
t24 = cos(qJ(3));
t16 = qJ(3) + pkin(11);
t6 = sin(t16);
t8 = cos(t16);
t43 = -m(4) * pkin(2) - t24 * mrSges(4,1) + t21 * mrSges(4,2) + t49 * t6 + t50 * t8 - mrSges(3,1);
t42 = m(4) * pkin(7) + t11 * mrSges(7,1) + t23 * mrSges(6,2) + t12 * mrSges(7,2) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + (m(7) * pkin(5) + mrSges(6,1)) * t20;
t40 = pkin(6) + r_base(3);
t22 = sin(qJ(1));
t39 = t22 * pkin(1) + r_base(2);
t25 = cos(qJ(1));
t38 = t25 * pkin(1) + r_base(1);
t10 = qJ(2) + t40;
t19 = -qJ(4) - pkin(7);
t17 = qJ(1) + pkin(10);
t9 = cos(t17);
t7 = sin(t17);
t5 = t24 * pkin(3) + pkin(2);
t1 = (-m(1) * r_base(3) - m(2) * t40 - t21 * mrSges(4,1) - t24 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t47 * t10 - t46 * (t21 * pkin(3) + t10) - t49 * t8 + t50 * t6) * g(3) + (-t22 * mrSges(2,1) - t25 * mrSges(2,2) - mrSges(1,2) + t48 * r_base(2) + t47 * t39 - t46 * (t9 * t19 + t7 * t5 + t39) + t42 * t9 + t43 * t7) * g(2) + (-t25 * mrSges(2,1) + t22 * mrSges(2,2) - mrSges(1,1) + t48 * r_base(1) + t47 * t38 - t46 * (t9 * t5 + t38) + t43 * t9 + (t46 * t19 - t42) * t7) * g(1);
U  = t1;
