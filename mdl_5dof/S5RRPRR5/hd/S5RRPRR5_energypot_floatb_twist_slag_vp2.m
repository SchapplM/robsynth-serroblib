% Calculate potential energy for
% S5RRPRR5
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
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:26
% EndTime: 2020-01-03 12:03:26
% DurationCPUTime: 0.44s
% Computational Cost: add. (153->58), mult. (109->44), div. (0->0), fcn. (73->10), ass. (0->24)
t33 = -m(1) - m(2);
t32 = -m(3) - m(4);
t31 = -m(5) - m(6) + t32;
t17 = sin(pkin(9));
t18 = cos(pkin(9));
t15 = pkin(9) + qJ(4);
t7 = qJ(5) + t15;
t2 = sin(t7);
t3 = cos(t7);
t4 = t18 * pkin(3) + pkin(2);
t5 = sin(t15);
t6 = cos(t15);
t30 = m(4) * pkin(2) + t18 * mrSges(4,1) - t17 * mrSges(4,2) + mrSges(3,1) + m(6) * (pkin(4) * t6 + t4) + t3 * mrSges(6,1) - t2 * mrSges(6,2) + m(5) * t4 + t6 * mrSges(5,1) - t5 * mrSges(5,2);
t19 = -pkin(7) - qJ(3);
t29 = m(4) * qJ(3) - m(5) * t19 - m(6) * (-pkin(8) + t19) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t28 = pkin(5) + r_base(1);
t12 = pkin(6) + t28;
t25 = t17 * pkin(3) + t12;
t21 = cos(qJ(1));
t20 = sin(qJ(1));
t16 = qJ(1) + qJ(2);
t9 = cos(t16);
t8 = sin(t16);
t1 = (mrSges(2,1) * t21 - t20 * mrSges(2,2) - mrSges(1,3) + t33 * r_base(3) + t30 * t9 + t31 * (-pkin(1) * t21 + r_base(3)) + t29 * t8) * g(3) + (-t20 * mrSges(2,1) - mrSges(2,2) * t21 - mrSges(1,2) + t33 * r_base(2) + t31 * (t20 * pkin(1) + r_base(2)) + t29 * t9 - t30 * t8) * g(2) + (-m(1) * r_base(1) - mrSges(1,1) - m(2) * t28 - mrSges(2,3) - mrSges(3,3) - mrSges(4,1) * t17 - mrSges(4,2) * t18 - m(5) * t25 - t5 * mrSges(5,1) - t6 * mrSges(5,2) - m(6) * (pkin(4) * t5 + t25) - t2 * mrSges(6,1) - t3 * mrSges(6,2) + t32 * t12) * g(1);
U = t1;
