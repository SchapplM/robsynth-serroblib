% Calculate potential energy for
% S5RRPRP5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRP5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRP5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP5_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP5_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:53:56
% EndTime: 2019-12-31 19:53:57
% DurationCPUTime: 0.30s
% Computational Cost: add. (157->55), mult. (129->41), div. (0->0), fcn. (93->8), ass. (0->24)
t36 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t35 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t34 = -m(2) - m(3);
t33 = -m(5) - m(6);
t32 = -m(1) - m(4) + t34;
t15 = qJ(2) + pkin(8);
t10 = cos(t15);
t18 = sin(qJ(2));
t20 = cos(qJ(2));
t11 = qJ(4) + t15;
t6 = sin(t11);
t7 = cos(t11);
t8 = t20 * pkin(2) + pkin(1);
t9 = sin(t15);
t31 = -m(3) * pkin(1) - m(4) * t8 - t20 * mrSges(3,1) - t10 * mrSges(4,1) + t18 * mrSges(3,2) + t9 * mrSges(4,2) - t35 * t6 + t36 * t7 - mrSges(2,1);
t17 = -qJ(3) - pkin(6);
t30 = m(3) * pkin(6) - m(4) * t17 - mrSges(2,2) + mrSges(6,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t16 = pkin(5) + r_base(3);
t28 = t18 * pkin(2) + t16;
t21 = cos(qJ(1));
t19 = sin(qJ(1));
t14 = -pkin(7) + t17;
t3 = pkin(3) * t10 + t8;
t1 = (-m(1) * r_base(3) - m(4) * t28 - mrSges(3,1) * t18 - t9 * mrSges(4,1) - mrSges(3,2) * t20 - t10 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) + t33 * (pkin(3) * t9 + t28) + t35 * t7 + t36 * t6 + t34 * t16) * g(3) + (-mrSges(1,2) + t33 * (t21 * t14 + t19 * t3 + r_base(2)) + t32 * r_base(2) + t30 * t21 + t31 * t19) * g(2) + (-mrSges(1,1) + t33 * (-t19 * t14 + t21 * t3 + r_base(1)) + t32 * r_base(1) + t31 * t21 - t30 * t19) * g(1);
U = t1;
