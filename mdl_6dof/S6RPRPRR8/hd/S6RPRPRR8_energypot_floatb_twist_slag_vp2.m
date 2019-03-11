% Calculate potential energy for
% S6RPRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR8_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRR8_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR8_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR8_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR8_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:57:39
% EndTime: 2019-03-09 03:57:40
% DurationCPUTime: 0.51s
% Computational Cost: add. (194->66), mult. (194->48), div. (0->0), fcn. (160->10), ass. (0->24)
t15 = sin(qJ(5));
t18 = cos(qJ(5));
t13 = qJ(5) + qJ(6);
t5 = sin(t13);
t6 = cos(t13);
t45 = mrSges(5,1) + m(6) * pkin(4) + mrSges(6,1) * t18 - mrSges(6,2) * t15 + m(7) * (pkin(5) * t18 + pkin(4)) + mrSges(7,1) * t6 - mrSges(7,2) * t5;
t44 = mrSges(5,2) - m(6) * pkin(8) + m(7) * (-pkin(9) - pkin(8)) - mrSges(6,3) - mrSges(7,3);
t43 = m(5) + m(6) + m(7);
t16 = sin(qJ(3));
t19 = cos(qJ(3));
t11 = qJ(3) + pkin(10);
t3 = sin(t11);
t4 = cos(t11);
t42 = mrSges(4,1) * t16 + mrSges(4,2) * t19 + t45 * t3 + t44 * t4 - mrSges(2,2) + mrSges(3,3);
t41 = -m(1) - m(2);
t40 = m(3) + m(4);
t35 = -m(4) * pkin(7) - t5 * mrSges(7,1) - mrSges(6,2) * t18 - t6 * mrSges(7,2) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t15 + t43 * (-qJ(4) - pkin(7));
t34 = pkin(3) * t16;
t12 = pkin(6) + r_base(3);
t32 = pkin(2) + t12;
t17 = sin(qJ(1));
t20 = cos(qJ(1));
t31 = t20 * pkin(1) + t17 * qJ(2) + r_base(1);
t1 = (-m(1) * r_base(3) - m(4) * t32 - t19 * mrSges(4,1) + t16 * mrSges(4,2) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t12 - t45 * t4 - t43 * (t19 * pkin(3) + t32) + t44 * t3) * g(3) + (-mrSges(1,2) + t41 * r_base(2) + (-t43 - t40) * (t17 * pkin(1) + r_base(2)) + (-t43 * (-qJ(2) - t34) + t40 * qJ(2) + t42) * t20 + t35 * t17) * g(2) + (-mrSges(1,1) + t41 * r_base(1) - t40 * t31 - t43 * (t17 * t34 + t31) + t35 * t20 - t42 * t17) * g(1);
U  = t1;
