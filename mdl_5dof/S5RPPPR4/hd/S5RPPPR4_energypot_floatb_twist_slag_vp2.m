% Calculate potential energy for
% S5RPPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPPR4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR4_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:04
% EndTime: 2019-12-31 17:45:04
% DurationCPUTime: 0.28s
% Computational Cost: add. (135->55), mult. (101->40), div. (0->0), fcn. (65->8), ass. (0->20)
t30 = -m(1) - m(2);
t29 = m(4) + m(5) + m(6);
t13 = sin(pkin(8));
t14 = cos(pkin(8));
t11 = pkin(8) + qJ(5);
t4 = sin(t11);
t6 = cos(t11);
t28 = t4 * mrSges(6,1) + t14 * mrSges(5,2) + t6 * mrSges(6,2) - mrSges(3,2) + mrSges(4,3) + (m(6) * pkin(4) + mrSges(5,1)) * t13;
t27 = -m(5) * qJ(4) + m(6) * (-pkin(6) - qJ(4)) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t25 = pkin(5) + r_base(3);
t16 = sin(qJ(1));
t24 = t16 * pkin(1) + r_base(2);
t17 = cos(qJ(1));
t23 = t17 * pkin(1) + r_base(1);
t8 = qJ(2) + t25;
t20 = pkin(3) + t8;
t12 = qJ(1) + pkin(7);
t7 = cos(t12);
t5 = sin(t12);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t25 - mrSges(2,3) - mrSges(3,3) - mrSges(4,1) - m(5) * t20 - t14 * mrSges(5,1) + t13 * mrSges(5,2) - m(6) * (t14 * pkin(4) + t20) - t6 * mrSges(6,1) + t4 * mrSges(6,2) + (-m(3) - m(4)) * t8) * g(3) + (-m(3) * t24 - t16 * mrSges(2,1) - mrSges(2,2) * t17 - mrSges(1,2) + t30 * r_base(2) - t29 * (t5 * pkin(2) + t24) + (qJ(3) * t29 + t28) * t7 + t27 * t5) * g(2) + (-m(3) * t23 - mrSges(2,1) * t17 + t16 * mrSges(2,2) - mrSges(1,1) + t30 * r_base(1) - t29 * (t7 * pkin(2) + t5 * qJ(3) + t23) + t27 * t7 - t28 * t5) * g(1);
U = t1;
