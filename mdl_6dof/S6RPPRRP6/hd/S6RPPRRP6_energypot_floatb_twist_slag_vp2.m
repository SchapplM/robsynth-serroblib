% Calculate potential energy for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRP6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP6_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:10:07
% EndTime: 2019-03-09 02:10:07
% DurationCPUTime: 0.39s
% Computational Cost: add. (141->69), mult. (190->62), div. (0->0), fcn. (160->6), ass. (0->29)
t47 = -m(1) - m(2);
t46 = -m(6) - m(7);
t45 = mrSges(6,3) + mrSges(7,2);
t18 = sin(qJ(4));
t21 = cos(qJ(4));
t44 = -t18 * mrSges(5,1) - t21 * mrSges(5,2) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t43 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2) + mrSges(5,3);
t42 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t41 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t19 = sin(qJ(1));
t40 = t19 * t18;
t20 = cos(qJ(5));
t39 = t19 * t20;
t38 = t19 * t21;
t22 = cos(qJ(1));
t37 = t22 * t18;
t36 = t22 * t20;
t35 = t22 * t21;
t16 = pkin(6) + r_base(3);
t34 = pkin(2) + t16;
t33 = t22 * pkin(1) + t19 * qJ(2) + r_base(1);
t32 = pkin(3) + t34;
t31 = t22 * qJ(3) + t33;
t29 = t19 * pkin(1) - t22 * qJ(2) + r_base(2);
t27 = t19 * qJ(3) + t29;
t26 = -t19 * pkin(7) + t31;
t25 = t22 * pkin(7) + t27;
t17 = sin(qJ(5));
t1 = (-m(1) * r_base(3) - m(4) * t34 - m(5) * t32 - mrSges(3,1) - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) + t46 * (t21 * pkin(4) + t18 * pkin(8) + t32) + (-m(2) - m(3)) * t16 + (t41 * t17 + t42 * t20 - mrSges(5,1)) * t21 + (mrSges(5,2) - t45) * t18) * g(3) + (-m(3) * t29 - m(4) * t27 - m(5) * t25 - mrSges(1,2) + t47 * r_base(2) + t45 * t38 + t46 * (pkin(4) * t40 - pkin(8) * t38 + t25) + t42 * (t22 * t17 + t18 * t39) + t41 * (t17 * t40 - t36) - t43 * t22 + t44 * t19) * g(2) + (-m(3) * t33 - m(4) * t31 - m(5) * t26 - mrSges(1,1) + t47 * r_base(1) + t46 * (pkin(4) * t37 - pkin(8) * t35 + t26) + t42 * (-t19 * t17 + t18 * t36) + t45 * t35 + t41 * (t17 * t37 + t39) + t44 * t22 + t43 * t19) * g(1);
U  = t1;
