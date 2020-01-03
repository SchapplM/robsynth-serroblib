% Calculate potential energy for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% m_mdh [6x1]
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
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP8_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRP8_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP8_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:15
% EndTime: 2019-12-31 20:03:16
% DurationCPUTime: 0.49s
% Computational Cost: add. (126->68), mult. (191->62), div. (0->0), fcn. (167->6), ass. (0->33)
t52 = -mrSges(3,1) - mrSges(4,1);
t51 = mrSges(3,2) - mrSges(4,3);
t46 = -m(4) - m(5);
t50 = -m(6) + t46;
t19 = sin(qJ(2));
t22 = cos(qJ(2));
t21 = cos(qJ(4));
t18 = sin(qJ(4));
t38 = t19 * t18;
t27 = t22 * t21 + t38;
t43 = -mrSges(5,2) - mrSges(6,2);
t44 = -mrSges(5,1) - mrSges(6,1);
t6 = -t22 * t18 + t19 * t21;
t49 = t51 * t19 + t52 * t22 + t27 * t44 + t6 * t43 - mrSges(2,1);
t48 = m(5) * pkin(3);
t47 = -m(1) - m(2);
t20 = sin(qJ(1));
t35 = qJ(3) * t19;
t37 = t20 * t22;
t45 = pkin(2) * t37 + t20 * t35;
t39 = mrSges(4,2) + mrSges(3,3) - mrSges(2,2) - mrSges(5,3) - mrSges(6,3);
t23 = cos(qJ(1));
t36 = t23 * t22;
t16 = pkin(5) + r_base(3);
t34 = pkin(4) * t38;
t33 = t20 * pkin(1) + r_base(2);
t32 = t19 * pkin(2) + t16;
t31 = t23 * pkin(1) + t20 * pkin(6) + r_base(1);
t30 = t33 + t45;
t26 = -t23 * pkin(6) + t33;
t17 = -qJ(5) - pkin(7);
t11 = t21 * pkin(4) + pkin(3);
t1 = (-m(1) * r_base(3) - m(6) * t32 - mrSges(1,3) - mrSges(2,3) + t44 * t6 - t43 * t27 + t46 * (-t22 * qJ(3) + t32) + (-m(6) * (-pkin(4) * t18 - qJ(3)) - t51) * t22 + (-m(6) * t11 - t48 + t52) * t19 + (-m(2) - m(3)) * t16) * g(3) + (-mrSges(1,2) - m(3) * t26 - m(4) * (t26 + t45) - m(5) * (pkin(3) * t37 + t30) - m(6) * (t11 * t37 + t30) + t47 * r_base(2)) * g(2) + (-t36 * t48 - m(3) * t31 - mrSges(1,1) + t47 * r_base(1) + t50 * (pkin(2) * t36 + t31)) * g(1) + ((-m(5) * (-pkin(6) + pkin(7)) - m(6) * (-pkin(6) - t17) + t39) * g(2) + (t50 * t35 - m(6) * (t11 * t22 + t34) + t49) * g(1)) * t23 + ((-m(6) * t34 + t49) * g(2) + (m(5) * pkin(7) - m(6) * t17 - t39) * g(1)) * t20;
U = t1;
