% Calculate potential energy for
% S5PRPRR5
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRPRR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR5_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR5_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:53:20
% EndTime: 2019-12-05 15:53:20
% DurationCPUTime: 0.40s
% Computational Cost: add. (170->56), mult. (192->45), div. (0->0), fcn. (168->10), ass. (0->22)
t16 = sin(pkin(9));
t18 = cos(pkin(9));
t15 = pkin(9) + qJ(4);
t8 = qJ(5) + t15;
t3 = sin(t8);
t4 = cos(t8);
t5 = t18 * pkin(3) + pkin(2);
t6 = sin(t15);
t7 = cos(t15);
t44 = -m(4) * pkin(2) - t18 * mrSges(4,1) + t16 * mrSges(4,2) - mrSges(3,1) - m(6) * (pkin(4) * t7 + t5) - t4 * mrSges(6,1) + t3 * mrSges(6,2) - m(5) * t5 - t7 * mrSges(5,1) + t6 * mrSges(5,2);
t20 = -pkin(6) - qJ(3);
t43 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) + m(6) * (-pkin(7) + t20) - mrSges(6,3) + m(5) * t20 - mrSges(5,3);
t42 = -m(1) - m(2);
t41 = m(3) + m(4) + m(5) + m(6);
t21 = sin(qJ(2));
t22 = cos(qJ(2));
t37 = t43 * t21 + t44 * t22 - mrSges(2,1);
t35 = t16 * pkin(3);
t36 = m(5) * t35 + m(6) * (pkin(4) * t6 + t35) - mrSges(2,2) + mrSges(3,3) + t16 * mrSges(4,1) + t18 * mrSges(4,2) + t3 * mrSges(6,1) + t4 * mrSges(6,2) + t6 * mrSges(5,1) + t7 * mrSges(5,2);
t19 = cos(pkin(8));
t17 = sin(pkin(8));
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(2) - t41) * (qJ(1) + r_base(3)) - t43 * t22 + t44 * t21) * g(3) + (-mrSges(1,2) + t42 * r_base(2) - t41 * (t17 * pkin(1) + r_base(2)) + (t41 * pkin(5) + t36) * t19 + t37 * t17) * g(2) + (-mrSges(1,1) + t42 * r_base(1) - t41 * (t19 * pkin(1) + t17 * pkin(5) + r_base(1)) + t37 * t19 - t36 * t17) * g(1);
U = t1;
