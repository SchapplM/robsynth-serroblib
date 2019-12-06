% Calculate potential energy for
% S5RPRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRPR4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR4_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:53:05
% EndTime: 2019-12-05 17:53:05
% DurationCPUTime: 0.30s
% Computational Cost: add. (153->58), mult. (109->44), div. (0->0), fcn. (73->10), ass. (0->24)
t33 = -m(1) - m(2);
t32 = -m(3) - m(4);
t31 = -m(5) - m(6) + t32;
t18 = sin(qJ(3));
t15 = qJ(3) + pkin(9);
t9 = qJ(5) + t15;
t2 = sin(t9);
t20 = cos(qJ(3));
t3 = cos(t9);
t4 = t20 * pkin(3) + pkin(2);
t5 = sin(t15);
t7 = cos(t15);
t30 = m(4) * pkin(2) + t20 * mrSges(4,1) - t18 * mrSges(4,2) + mrSges(3,1) + m(6) * (pkin(4) * t7 + t4) + t3 * mrSges(6,1) - t2 * mrSges(6,2) + m(5) * t4 + t7 * mrSges(5,1) - t5 * mrSges(5,2);
t17 = -qJ(4) - pkin(6);
t29 = -m(4) * pkin(6) + m(5) * t17 + m(6) * (-pkin(7) + t17) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t28 = pkin(5) + r_base(1);
t10 = qJ(2) + t28;
t25 = t18 * pkin(3) + t10;
t21 = cos(qJ(1));
t19 = sin(qJ(1));
t16 = qJ(1) + pkin(8);
t8 = cos(t16);
t6 = sin(t16);
t1 = (-mrSges(2,1) * t21 + t19 * mrSges(2,2) - mrSges(1,3) + t33 * r_base(3) - t30 * t8 + t31 * (t21 * pkin(1) + r_base(3)) + t29 * t6) * g(3) + (t19 * mrSges(2,1) + mrSges(2,2) * t21 - mrSges(1,2) + t33 * r_base(2) + t31 * (-pkin(1) * t19 + r_base(2)) + t29 * t8 + t30 * t6) * g(2) + (-m(1) * r_base(1) - mrSges(1,1) - m(2) * t28 - mrSges(2,3) - mrSges(3,3) - mrSges(4,1) * t18 - mrSges(4,2) * t20 - m(5) * t25 - t5 * mrSges(5,1) - t7 * mrSges(5,2) - m(6) * (pkin(4) * t5 + t25) - t2 * mrSges(6,1) - t3 * mrSges(6,2) + t32 * t10) * g(1);
U = t1;
