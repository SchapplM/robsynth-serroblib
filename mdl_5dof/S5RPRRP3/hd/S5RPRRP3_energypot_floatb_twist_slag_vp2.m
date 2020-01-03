% Calculate potential energy for
% S5RPRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRP3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP3_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:46:49
% EndTime: 2020-01-03 11:46:49
% DurationCPUTime: 0.33s
% Computational Cost: add. (147->53), mult. (109->37), div. (0->0), fcn. (73->8), ass. (0->22)
t34 = m(5) + m(6);
t33 = -mrSges(5,2) - mrSges(6,2);
t32 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t31 = -m(1) - m(2);
t30 = -m(3) - m(4);
t28 = -t34 + t30;
t14 = sin(qJ(3));
t16 = cos(qJ(3));
t13 = qJ(3) + qJ(4);
t6 = sin(t13);
t7 = cos(t13);
t27 = m(4) * pkin(2) + t16 * mrSges(4,1) - t14 * mrSges(4,2) + mrSges(3,1) + t34 * (t16 * pkin(3) + pkin(2)) + t32 * t7 + t33 * t6;
t18 = -pkin(7) - pkin(6);
t26 = m(4) * pkin(6) - m(5) * t18 - m(6) * (-qJ(5) + t18) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t25 = pkin(5) + r_base(1);
t5 = qJ(2) + t25;
t17 = cos(qJ(1));
t15 = sin(qJ(1));
t12 = qJ(1) + pkin(8);
t4 = cos(t12);
t3 = sin(t12);
t1 = (mrSges(2,1) * t17 - mrSges(2,2) * t15 - mrSges(1,3) + t31 * r_base(3) + t27 * t4 + t28 * (-t17 * pkin(1) + r_base(3)) + t26 * t3) * g(3) + (-mrSges(2,1) * t15 - mrSges(2,2) * t17 - mrSges(1,2) + t31 * r_base(2) + t28 * (t15 * pkin(1) + r_base(2)) + t26 * t4 - t27 * t3) * g(2) + (-m(1) * r_base(1) - m(2) * t25 - mrSges(4,1) * t14 - mrSges(4,2) * t16 - mrSges(1,1) - mrSges(2,3) - mrSges(3,3) + t33 * t7 - t34 * (t14 * pkin(3) + t5) - t32 * t6 + t30 * t5) * g(1);
U = t1;
