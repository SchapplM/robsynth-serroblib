% Calculate potential energy for
% S6RPPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRR7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRR7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR7_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR7_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:32:53
% EndTime: 2019-03-09 02:32:53
% DurationCPUTime: 0.46s
% Computational Cost: add. (193->71), mult. (170->53), div. (0->0), fcn. (132->10), ass. (0->28)
t43 = m(6) + m(7);
t19 = sin(qJ(6));
t21 = cos(qJ(6));
t42 = -m(7) * pkin(5) - mrSges(7,1) * t21 + mrSges(7,2) * t19 - mrSges(6,1);
t41 = -m(1) - m(2);
t39 = m(3) + m(4) + m(5);
t38 = -m(7) * pkin(9) - mrSges(7,3);
t16 = sin(pkin(10));
t17 = cos(pkin(10));
t35 = t16 * pkin(3);
t14 = pkin(10) + qJ(4);
t8 = qJ(5) + t14;
t4 = sin(t8);
t5 = cos(t8);
t6 = sin(t14);
t7 = cos(t14);
t37 = -m(5) * t35 - t16 * mrSges(4,1) - t6 * mrSges(5,1) - t17 * mrSges(4,2) - t7 * mrSges(5,2) - t5 * mrSges(6,2) + t4 * t42 + mrSges(2,2) - mrSges(3,3);
t18 = -pkin(7) - qJ(3);
t36 = -m(4) * qJ(3) + m(5) * t18 - t19 * mrSges(7,1) - t21 * mrSges(7,2) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3) + t43 * (-pkin(8) + t18);
t2 = pkin(4) * t6 + t35;
t34 = -qJ(2) - t2;
t15 = pkin(6) + r_base(3);
t32 = pkin(2) + t15;
t20 = sin(qJ(1));
t22 = cos(qJ(1));
t31 = pkin(1) * t22 + qJ(2) * t20 + r_base(1);
t29 = pkin(3) * t17 + t32;
t1 = (-m(1) * r_base(3) - m(4) * t32 - m(5) * t29 - t17 * mrSges(4,1) - t7 * mrSges(5,1) + t16 * mrSges(4,2) + t6 * mrSges(5,2) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + t42 * t5 - t43 * (pkin(4) * t7 + t29) + (mrSges(6,2) + t38) * t4 + (-m(2) - m(3)) * t15) * g(3) + (-mrSges(1,2) + t41 * r_base(2) + (-t39 - t43) * (pkin(1) * t20 + r_base(2)) + (-m(6) * t34 - m(7) * (pkin(9) * t5 + t34) - t5 * mrSges(7,3) + t39 * qJ(2) - t37) * t22 + t36 * t20) * g(2) + (-mrSges(1,1) + t41 * r_base(1) - t43 * (t2 * t20 + t31) - t39 * t31 + t36 * t22 + (-t38 * t5 + t37) * t20) * g(1);
U  = t1;
