% Calculate potential energy for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR8_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRR8_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR8_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR8_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:39
% EndTime: 2019-12-31 19:05:39
% DurationCPUTime: 0.29s
% Computational Cost: add. (124->55), mult. (153->44), div. (0->0), fcn. (141->8), ass. (0->21)
t35 = -m(1) - m(2);
t34 = -mrSges(2,1) - mrSges(3,1);
t33 = mrSges(2,2) - mrSges(3,3);
t32 = -m(4) - m(5) - m(6);
t17 = sin(qJ(4));
t18 = cos(qJ(4));
t16 = qJ(4) + qJ(5);
t7 = sin(t16);
t8 = cos(t16);
t31 = m(5) * pkin(3) + t18 * mrSges(5,1) - t17 * mrSges(5,2) + mrSges(4,1) + m(6) * (pkin(4) * t18 + pkin(3)) + t8 * mrSges(6,1) - t7 * mrSges(6,2);
t30 = m(5) * pkin(7) - m(6) * (-pkin(8) - pkin(7)) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t29 = cos(qJ(3));
t28 = sin(qJ(1));
t27 = sin(qJ(3));
t15 = pkin(5) + r_base(3);
t19 = cos(qJ(1));
t26 = t19 * pkin(1) + t28 * qJ(2) + r_base(1);
t24 = t28 * pkin(1) - t19 * qJ(2) + r_base(2);
t2 = t19 * t27 - t28 * t29;
t1 = -t19 * t29 - t28 * t27;
t3 = (-m(1) * r_base(3) + t7 * mrSges(6,1) + mrSges(5,2) * t18 + t8 * mrSges(6,2) - mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + (m(6) * pkin(4) + mrSges(5,1)) * t17 + (-m(2) - m(3)) * t15 + t32 * (-pkin(6) + t15)) * g(3) + (-m(3) * t24 - mrSges(1,2) + t35 * r_base(2) + t34 * t28 - t33 * t19 + t32 * (t28 * pkin(2) + t24) + t31 * t2 + t30 * t1) * g(2) + (-m(3) * t26 - mrSges(1,1) + t35 * r_base(1) + t33 * t28 + t34 * t19 + t32 * (t19 * pkin(2) + t26) - t30 * t2 + t31 * t1) * g(1);
U = t3;
