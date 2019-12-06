% Calculate potential energy for
% S5RPRPR1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRPR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR1_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR1_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:46:59
% EndTime: 2019-12-05 17:46:59
% DurationCPUTime: 0.29s
% Computational Cost: add. (120->55), mult. (115->41), div. (0->0), fcn. (79->8), ass. (0->20)
t30 = -m(1) - m(2);
t29 = m(3) + m(4) + m(5) + m(6);
t15 = sin(qJ(3));
t17 = cos(qJ(3));
t12 = qJ(3) + pkin(8);
t6 = qJ(5) + t12;
t2 = sin(t6);
t26 = pkin(3) * t15;
t3 = cos(t6);
t4 = sin(t12);
t5 = cos(t12);
t28 = m(5) * t26 + m(6) * (pkin(4) * t4 + t26) - mrSges(2,2) + mrSges(3,3) + mrSges(4,1) * t15 + mrSges(4,2) * t17 + t2 * mrSges(6,1) + t3 * mrSges(6,2) + t4 * mrSges(5,1) + t5 * mrSges(5,2);
t14 = -qJ(4) - pkin(6);
t27 = -m(4) * pkin(6) + m(5) * t14 + m(6) * (-pkin(7) + t14) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t13 = pkin(5) + r_base(3);
t24 = pkin(2) + t13;
t22 = t17 * pkin(3) + t24;
t18 = cos(qJ(1));
t16 = sin(qJ(1));
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - mrSges(3,1) - m(4) * t24 - t17 * mrSges(4,1) + t15 * mrSges(4,2) - m(5) * t22 - t5 * mrSges(5,1) + t4 * mrSges(5,2) - m(6) * (pkin(4) * t5 + t22) - t3 * mrSges(6,1) + t2 * mrSges(6,2) + (-m(2) - m(3)) * t13) * g(3) + (-mrSges(1,2) + t30 * r_base(2) - t29 * (t16 * pkin(1) + r_base(2)) + (qJ(2) * t29 + t28) * t18 + t27 * t16) * g(2) + (-mrSges(1,1) + t30 * r_base(1) - t29 * (t18 * pkin(1) + t16 * qJ(2) + r_base(1)) + t27 * t18 - t28 * t16) * g(1);
U = t1;
