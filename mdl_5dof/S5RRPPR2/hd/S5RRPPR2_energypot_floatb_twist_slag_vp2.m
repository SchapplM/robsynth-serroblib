% Calculate potential energy for
% S5RRPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPPR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:22
% EndTime: 2022-01-20 10:05:23
% DurationCPUTime: 0.33s
% Computational Cost: add. (170->57), mult. (117->45), div. (0->0), fcn. (85->10), ass. (0->25)
t17 = sin(qJ(5));
t19 = cos(qJ(5));
t39 = -m(6) * pkin(4) - t19 * mrSges(6,1) + t17 * mrSges(6,2) - mrSges(5,1);
t38 = -m(6) * pkin(7) + mrSges(5,2) - mrSges(6,3);
t37 = -m(1) - m(2);
t36 = m(5) + m(6);
t15 = sin(pkin(9));
t16 = cos(pkin(9));
t35 = t38 * t15 + t39 * t16 - mrSges(4,1);
t34 = -t17 * mrSges(6,1) - t19 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t14 = qJ(1) + qJ(2);
t32 = pkin(5) + r_base(3);
t18 = sin(qJ(1));
t31 = t18 * pkin(1) + r_base(2);
t20 = cos(qJ(1));
t30 = t20 * pkin(1) + r_base(1);
t29 = pkin(6) + t32;
t10 = sin(t14);
t28 = pkin(2) * t10 + t31;
t11 = cos(t14);
t27 = pkin(2) * t11 + t30;
t9 = pkin(8) + t14;
t5 = cos(t9);
t4 = sin(t9);
t1 = (-m(1) * r_base(3) - m(2) * t32 - m(3) * t29 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) + (-m(4) - t36) * (qJ(3) + t29) - t38 * t16 + t39 * t15) * g(3) + (-m(3) * t31 - m(4) * t28 - t18 * mrSges(2,1) - t10 * mrSges(3,1) - t20 * mrSges(2,2) - t11 * mrSges(3,2) - mrSges(1,2) + t37 * r_base(2) - t36 * (t4 * pkin(3) + t28) + (t36 * qJ(4) - t34) * t5 + t35 * t4) * g(2) + (-m(3) * t30 - m(4) * t27 - t20 * mrSges(2,1) - t11 * mrSges(3,1) + t18 * mrSges(2,2) + t10 * mrSges(3,2) - mrSges(1,1) + t37 * r_base(1) - t36 * (t5 * pkin(3) + t4 * qJ(4) + t27) + t35 * t5 + t34 * t4) * g(1);
U = t1;
