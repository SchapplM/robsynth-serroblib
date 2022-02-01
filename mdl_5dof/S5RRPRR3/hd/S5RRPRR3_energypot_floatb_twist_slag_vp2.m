% Calculate potential energy for
% S5RRPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR3_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:33:55
% EndTime: 2022-01-20 10:33:55
% DurationCPUTime: 0.26s
% Computational Cost: add. (152->57), mult. (86->44), div. (0->0), fcn. (50->10), ass. (0->25)
t34 = -m(1) - m(2);
t33 = -m(5) - m(6);
t17 = sin(qJ(5));
t19 = cos(qJ(5));
t32 = -m(6) * pkin(4) - t19 * mrSges(6,1) + t17 * mrSges(6,2) - mrSges(5,1);
t31 = m(6) * pkin(8) - mrSges(5,2) + mrSges(6,3);
t16 = qJ(1) + qJ(2);
t30 = pkin(5) + r_base(3);
t18 = sin(qJ(1));
t29 = t18 * pkin(1) + r_base(2);
t20 = cos(qJ(1));
t28 = t20 * pkin(1) + r_base(1);
t27 = pkin(6) + t30;
t12 = sin(t16);
t26 = pkin(2) * t12 + t29;
t13 = cos(t16);
t25 = pkin(2) * t13 + t28;
t11 = pkin(9) + t16;
t22 = qJ(3) + t27;
t10 = qJ(4) + t11;
t7 = cos(t11);
t6 = sin(t11);
t4 = cos(t10);
t3 = sin(t10);
t1 = (-m(1) * r_base(3) - m(2) * t30 - m(3) * t27 - m(4) * t22 - mrSges(6,1) * t17 - t19 * mrSges(6,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) + t33 * (pkin(7) + t22)) * g(3) + (-m(3) * t29 - m(4) * t26 - t18 * mrSges(2,1) - t12 * mrSges(3,1) - t6 * mrSges(4,1) - mrSges(2,2) * t20 - t13 * mrSges(3,2) - t7 * mrSges(4,2) - mrSges(1,2) + t34 * r_base(2) + t33 * (pkin(3) * t6 + t26) + t31 * t4 + t32 * t3) * g(2) + (-m(3) * t28 - m(4) * t25 - mrSges(2,1) * t20 - t13 * mrSges(3,1) - t7 * mrSges(4,1) + t18 * mrSges(2,2) + t12 * mrSges(3,2) + t6 * mrSges(4,2) - mrSges(1,1) + t34 * r_base(1) + t32 * t4 + t33 * (pkin(3) * t7 + t25) - t31 * t3) * g(1);
U = t1;
