% Calculate potential energy for
% S5PRRRP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP1_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:39:42
% EndTime: 2019-12-05 16:39:42
% DurationCPUTime: 0.29s
% Computational Cost: add. (148->52), mult. (97->37), div. (0->0), fcn. (61->8), ass. (0->22)
t32 = -m(5) - m(6);
t31 = mrSges(5,2) + mrSges(6,2);
t30 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t29 = -m(1) - m(2);
t28 = -m(4) + t32;
t16 = sin(qJ(4));
t17 = cos(qJ(4));
t27 = t32 * pkin(3) + t31 * t16 + t30 * t17 - mrSges(4,1);
t26 = m(5) * pkin(7) - m(6) * (-qJ(5) - pkin(7)) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t12 = pkin(8) + qJ(2);
t13 = sin(pkin(8));
t25 = t13 * pkin(1) + r_base(2);
t14 = cos(pkin(8));
t24 = t14 * pkin(1) + r_base(1);
t23 = qJ(1) + r_base(3);
t20 = pkin(5) + t23;
t9 = qJ(3) + t12;
t8 = cos(t12);
t7 = sin(t12);
t4 = cos(t9);
t3 = sin(t9);
t1 = (-m(1) * r_base(3) - m(2) * t23 - m(3) * t20 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) - t31 * t17 + t28 * (pkin(6) + t20) + t30 * t16) * g(3) + (-m(3) * t25 - mrSges(2,1) * t13 - t7 * mrSges(3,1) - mrSges(2,2) * t14 - t8 * mrSges(3,2) - mrSges(1,2) + t29 * r_base(2) + t28 * (pkin(2) * t7 + t25) + t26 * t4 + t27 * t3) * g(2) + (-m(3) * t24 - mrSges(2,1) * t14 - t8 * mrSges(3,1) + mrSges(2,2) * t13 + t7 * mrSges(3,2) - mrSges(1,1) + t29 * r_base(1) + t27 * t4 + t28 * (pkin(2) * t8 + t24) - t26 * t3) * g(1);
U = t1;
