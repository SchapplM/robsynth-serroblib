% Calculate potential energy for
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPPR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR1_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:00
% EndTime: 2020-01-03 11:20:01
% DurationCPUTime: 0.44s
% Computational Cost: add. (175->54), mult. (153->42), div. (0->0), fcn. (125->10), ass. (0->22)
t11 = sin(pkin(9));
t13 = cos(pkin(9));
t9 = pkin(9) + qJ(5);
t3 = sin(t9);
t5 = cos(t9);
t36 = -mrSges(4,1) - m(5) * pkin(3) - t13 * mrSges(5,1) + t11 * mrSges(5,2) - m(6) * (t13 * pkin(4) + pkin(3)) - t5 * mrSges(6,1) + t3 * mrSges(6,2);
t35 = mrSges(4,2) - m(5) * qJ(4) + m(6) * (-pkin(6) - qJ(4)) - mrSges(5,3) - mrSges(6,3);
t30 = m(4) + m(5) + m(6);
t12 = sin(pkin(8));
t14 = cos(pkin(8));
t34 = t35 * t12 + t36 * t14 - mrSges(3,1);
t33 = -m(1) - m(2);
t29 = -m(3) - t30;
t28 = t3 * mrSges(6,1) + t13 * mrSges(5,2) + t5 * mrSges(6,2) - mrSges(3,2) + mrSges(4,3) + (m(6) * pkin(4) + mrSges(5,1)) * t11 + t30 * qJ(3);
t26 = pkin(5) + r_base(1);
t16 = sin(qJ(1));
t25 = t16 * pkin(1) + r_base(2);
t17 = cos(qJ(1));
t10 = qJ(1) + pkin(7);
t6 = cos(t10);
t4 = sin(t10);
t1 = (t17 * mrSges(2,1) - t16 * mrSges(2,2) - mrSges(1,3) + t33 * r_base(3) + (t30 * pkin(2) - t34) * t6 + t29 * (-t17 * pkin(1) + r_base(3)) + t28 * t4) * g(3) + (-m(3) * t25 - t16 * mrSges(2,1) - t17 * mrSges(2,2) - mrSges(1,2) + t33 * r_base(2) - t30 * (t4 * pkin(2) + t25) + t28 * t6 + t34 * t4) * g(2) + (-m(1) * r_base(1) - m(2) * t26 - mrSges(1,1) - mrSges(2,3) - mrSges(3,3) + t29 * (qJ(2) + t26) - t35 * t14 + t36 * t12) * g(1);
U = t1;
