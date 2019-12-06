% Calculate potential energy for
% S5RRPRR6
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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRR6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR6_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:35:37
% EndTime: 2019-12-05 18:35:37
% DurationCPUTime: 0.42s
% Computational Cost: add. (175->56), mult. (153->44), div. (0->0), fcn. (125->10), ass. (0->22)
t15 = sin(qJ(4));
t17 = cos(qJ(4));
t11 = qJ(4) + qJ(5);
t5 = sin(t11);
t7 = cos(t11);
t38 = -mrSges(4,1) - m(5) * pkin(3) - t17 * mrSges(5,1) + t15 * mrSges(5,2) - m(6) * (t17 * pkin(4) + pkin(3)) - t7 * mrSges(6,1) + t5 * mrSges(6,2);
t37 = mrSges(4,2) - m(5) * pkin(7) + m(6) * (-pkin(8) - pkin(7)) - mrSges(5,3) - mrSges(6,3);
t13 = sin(pkin(9));
t14 = cos(pkin(9));
t36 = t37 * t13 + t38 * t14 - mrSges(3,1);
t35 = -m(1) - m(2);
t32 = m(4) + m(5) + m(6);
t31 = -t5 * mrSges(6,1) - t17 * mrSges(5,2) - t7 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3) + (-m(6) * pkin(4) - mrSges(5,1)) * t15;
t29 = pkin(5) + r_base(1);
t18 = cos(qJ(1));
t28 = t18 * pkin(1) + r_base(3);
t16 = sin(qJ(1));
t27 = -t16 * pkin(1) + r_base(2);
t12 = qJ(1) + qJ(2);
t8 = cos(t12);
t6 = sin(t12);
t1 = (-m(3) * t28 - t18 * mrSges(2,1) + t16 * mrSges(2,2) - mrSges(1,3) + t35 * r_base(3) - t32 * (t8 * pkin(2) + t6 * qJ(3) + t28) + t36 * t8 + t31 * t6) * g(3) + (-m(3) * t27 + t16 * mrSges(2,1) + t18 * mrSges(2,2) - mrSges(1,2) + t35 * r_base(2) - t32 * (t8 * qJ(3) + t27) + t31 * t8 + (t32 * pkin(2) - t36) * t6) * g(2) + (-m(1) * r_base(1) - m(2) * t29 - mrSges(1,1) - mrSges(2,3) - mrSges(3,3) + (-m(3) - t32) * (pkin(6) + t29) - t37 * t14 + t38 * t13) * g(1);
U = t1;
