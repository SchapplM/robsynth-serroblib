% Calculate potential energy for
% S5RPRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR5_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:03
% EndTime: 2019-12-05 18:16:03
% DurationCPUTime: 0.28s
% Computational Cost: add. (154->56), mult. (97->43), div. (0->0), fcn. (61->10), ass. (0->22)
t30 = -m(1) - m(2);
t29 = -m(4) - m(5) - m(6);
t13 = qJ(4) + qJ(5);
t10 = cos(t13);
t14 = sin(qJ(4));
t16 = cos(qJ(4));
t9 = sin(t13);
t28 = m(5) * pkin(3) + t16 * mrSges(5,1) - t14 * mrSges(5,2) + mrSges(4,1) + m(6) * (pkin(4) * t16 + pkin(3)) + t10 * mrSges(6,1) - t9 * mrSges(6,2);
t27 = -m(5) * pkin(7) + m(6) * (-pkin(8) - pkin(7)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t26 = pkin(5) + r_base(1);
t12 = qJ(1) + pkin(9);
t17 = cos(qJ(1));
t25 = t17 * pkin(1) + r_base(3);
t23 = qJ(2) + t26;
t15 = sin(qJ(1));
t22 = -pkin(1) * t15 + r_base(2);
t8 = qJ(3) + t12;
t7 = cos(t12);
t6 = sin(t12);
t3 = cos(t8);
t2 = sin(t8);
t1 = (-m(3) * t25 - mrSges(2,1) * t17 - t7 * mrSges(3,1) + mrSges(2,2) * t15 + t6 * mrSges(3,2) - mrSges(1,3) + t30 * r_base(3) - t28 * t3 + t29 * (pkin(2) * t7 + t25) + t27 * t2) * g(3) + (-m(3) * t22 + mrSges(2,1) * t15 + t6 * mrSges(3,1) + mrSges(2,2) * t17 + t7 * mrSges(3,2) - mrSges(1,2) + t30 * r_base(2) + t29 * (-pkin(2) * t6 + t22) + t27 * t3 + t28 * t2) * g(2) + (-m(1) * r_base(1) - m(2) * t26 - m(3) * t23 - t9 * mrSges(6,1) - mrSges(5,2) * t16 - t10 * mrSges(6,2) - mrSges(1,1) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) + (-m(6) * pkin(4) - mrSges(5,1)) * t14 + t29 * (pkin(6) + t23)) * g(1);
U = t1;
