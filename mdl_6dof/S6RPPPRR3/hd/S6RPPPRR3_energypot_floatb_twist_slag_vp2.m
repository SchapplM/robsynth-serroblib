% Calculate potential energy for
% S6RPPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPPRR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPPRR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:13
% EndTime: 2019-03-09 01:33:14
% DurationCPUTime: 0.41s
% Computational Cost: add. (196->68), mult. (243->55), div. (0->0), fcn. (247->10), ass. (0->31)
t24 = sin(qJ(6));
t25 = cos(qJ(6));
t51 = m(7) * pkin(5) + t25 * mrSges(7,1) - t24 * mrSges(7,2) + mrSges(6,1);
t50 = -m(7) * pkin(8) + mrSges(6,2) - mrSges(7,3);
t49 = -m(1) - m(2);
t48 = -m(4) - m(5);
t47 = -m(6) - m(7);
t46 = -mrSges(2,1) - mrSges(3,1);
t45 = mrSges(2,2) - mrSges(3,3);
t19 = pkin(10) + qJ(5);
t11 = sin(t19);
t12 = cos(t19);
t21 = sin(pkin(10));
t22 = cos(pkin(10));
t43 = m(5) * pkin(3) + t22 * mrSges(5,1) - t21 * mrSges(5,2) - t50 * t11 + t51 * t12 + mrSges(4,1);
t42 = m(5) * qJ(4) + t24 * mrSges(7,1) + t25 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t41 = cos(qJ(1));
t40 = sin(qJ(1));
t39 = cos(pkin(9));
t38 = sin(pkin(9));
t20 = pkin(6) + r_base(3);
t13 = -qJ(3) + t20;
t37 = t41 * pkin(1) + t40 * qJ(2) + r_base(1);
t36 = t41 * pkin(2) + t37;
t31 = t40 * pkin(1) - t41 * qJ(2) + r_base(2);
t28 = t40 * pkin(2) + t31;
t23 = -pkin(7) - qJ(4);
t10 = t22 * pkin(4) + pkin(3);
t6 = t41 * t38 - t40 * t39;
t5 = -t40 * t38 - t41 * t39;
t1 = (-m(1) * r_base(3) + t21 * mrSges(5,1) + t22 * mrSges(5,2) - mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + t47 * (-t21 * pkin(4) + t13) + (-m(2) - m(3)) * t20 + t48 * t13 + t50 * t12 + t51 * t11) * g(3) + (-m(3) * t31 - mrSges(1,2) + t49 * r_base(2) - t45 * t41 + t46 * t40 + t48 * t28 + t47 * (-t6 * t10 + t5 * t23 + t28) + t43 * t6 + t42 * t5) * g(2) + (-m(3) * t37 - mrSges(1,1) + t49 * r_base(1) + t46 * t41 + t45 * t40 + t48 * t36 + t47 * (-t5 * t10 - t6 * t23 + t36) - t42 * t6 + t43 * t5) * g(1);
U  = t1;
