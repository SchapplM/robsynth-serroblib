% Calculate potential energy for
% S6RPPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPPRR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPPRR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:31:37
% EndTime: 2019-03-09 01:31:37
% DurationCPUTime: 0.43s
% Computational Cost: add. (210->71), mult. (155->53), div. (0->0), fcn. (117->10), ass. (0->28)
t43 = m(6) + m(7);
t18 = sin(qJ(6));
t20 = cos(qJ(6));
t42 = -m(7) * pkin(5) - t20 * mrSges(7,1) + t18 * mrSges(7,2) - mrSges(6,1);
t41 = -m(1) - m(2);
t40 = m(4) + m(5);
t38 = -m(7) * pkin(8) - mrSges(7,3);
t15 = sin(pkin(10));
t16 = cos(pkin(10));
t13 = pkin(10) + qJ(5);
t5 = sin(t13);
t7 = cos(t13);
t37 = -t15 * mrSges(5,1) - t16 * mrSges(5,2) - t7 * mrSges(6,2) + t42 * t5 + mrSges(3,2) - mrSges(4,3);
t36 = -m(5) * qJ(4) - t18 * mrSges(7,1) - t20 * mrSges(7,2) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3) + t43 * (-pkin(7) - qJ(4));
t35 = pkin(4) * t15;
t34 = pkin(6) + r_base(3);
t19 = sin(qJ(1));
t33 = t19 * pkin(1) + r_base(2);
t21 = cos(qJ(1));
t32 = t21 * pkin(1) + r_base(1);
t30 = -qJ(3) - t35;
t9 = qJ(2) + t34;
t14 = qJ(1) + pkin(9);
t6 = sin(t14);
t8 = cos(t14);
t29 = t8 * pkin(2) + t6 * qJ(3) + t32;
t28 = pkin(3) + t9;
t1 = (-m(1) * r_base(3) - m(2) * t34 - m(5) * t28 - t16 * mrSges(5,1) + t15 * mrSges(5,2) - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + (-m(3) - m(4)) * t9 + t42 * t7 - t43 * (t16 * pkin(4) + t28) + (mrSges(6,2) + t38) * t5) * g(3) + (-m(3) * t33 - t19 * mrSges(2,1) - t21 * mrSges(2,2) - mrSges(1,2) + t41 * r_base(2) + (-t43 - t40) * (t6 * pkin(2) + t33) + (-m(6) * t30 - m(7) * (pkin(8) * t7 + t30) - t7 * mrSges(7,3) + t40 * qJ(3) - t37) * t8 + t36 * t6) * g(2) + (-m(3) * t32 - t21 * mrSges(2,1) + t19 * mrSges(2,2) - mrSges(1,1) + t41 * r_base(1) - t40 * t29 - t43 * (t6 * t35 + t29) + t36 * t8 + (-t38 * t7 + t37) * t6) * g(1);
U  = t1;
