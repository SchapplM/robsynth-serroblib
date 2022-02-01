% Calculate potential energy for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR3_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:04
% EndTime: 2022-01-23 09:14:05
% DurationCPUTime: 0.27s
% Computational Cost: add. (153->58), mult. (109->44), div. (0->0), fcn. (73->10), ass. (0->24)
t34 = -m(1) - m(2);
t33 = -m(3) - m(4);
t32 = -m(5) - m(6) + t33;
t18 = sin(pkin(9));
t19 = cos(pkin(9));
t16 = pkin(9) + qJ(4);
t9 = qJ(5) + t16;
t2 = sin(t9);
t3 = cos(t9);
t4 = t19 * pkin(3) + pkin(2);
t5 = sin(t16);
t7 = cos(t16);
t31 = -m(4) * pkin(2) - t19 * mrSges(4,1) + t18 * mrSges(4,2) - mrSges(3,1) - m(6) * (pkin(4) * t7 + t4) - t3 * mrSges(6,1) + t2 * mrSges(6,2) - m(5) * t4 - t7 * mrSges(5,1) + t5 * mrSges(5,2);
t20 = -pkin(6) - qJ(3);
t30 = m(4) * qJ(3) - m(5) * t20 + m(6) * (pkin(7) - t20) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t29 = pkin(5) + r_base(3);
t10 = qJ(2) + t29;
t26 = t18 * pkin(3) + t10;
t22 = cos(qJ(1));
t21 = sin(qJ(1));
t17 = qJ(1) + pkin(8);
t8 = cos(t17);
t6 = sin(t17);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t29 - mrSges(2,3) - mrSges(3,3) - t18 * mrSges(4,1) - t19 * mrSges(4,2) - m(5) * t26 - t5 * mrSges(5,1) - t7 * mrSges(5,2) - m(6) * (pkin(4) * t5 + t26) - t2 * mrSges(6,1) - t3 * mrSges(6,2) + t33 * t10) * g(3) + (-t21 * mrSges(2,1) - t22 * mrSges(2,2) - mrSges(1,2) + t34 * r_base(2) + t32 * (t21 * pkin(1) + r_base(2)) + t30 * t8 + t31 * t6) * g(2) + (-t22 * mrSges(2,1) + t21 * mrSges(2,2) - mrSges(1,1) + t34 * r_base(1) + t31 * t8 + t32 * (t22 * pkin(1) + r_base(1)) - t30 * t6) * g(1);
U = t1;
