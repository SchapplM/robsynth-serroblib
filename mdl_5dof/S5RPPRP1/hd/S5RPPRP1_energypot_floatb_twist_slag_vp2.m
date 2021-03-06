% Calculate potential energy for
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:04
% EndTime: 2020-01-03 11:25:04
% DurationCPUTime: 0.47s
% Computational Cost: add. (165->67), mult. (153->62), div. (0->0), fcn. (125->8), ass. (0->29)
t40 = -mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t39 = m(6) * pkin(4);
t38 = -m(1) - m(2);
t37 = -mrSges(5,1) - mrSges(6,1);
t36 = -mrSges(3,2) + mrSges(4,3);
t35 = mrSges(5,2) + mrSges(6,2);
t34 = -m(4) - m(5) - m(6);
t12 = sin(pkin(8));
t13 = cos(pkin(8));
t33 = t13 * mrSges(4,1) + t40 * t12 + mrSges(3,1);
t15 = sin(qJ(4));
t11 = qJ(1) + pkin(7);
t7 = sin(t11);
t32 = t7 * t15;
t8 = cos(t11);
t31 = t8 * t15;
t28 = t13 * t15;
t17 = cos(qJ(4));
t27 = t13 * t17;
t26 = pkin(5) + r_base(1);
t16 = sin(qJ(1));
t25 = t16 * pkin(1) + r_base(2);
t18 = cos(qJ(1));
t24 = -t18 * pkin(1) + r_base(3);
t23 = pkin(3) * t13 + pkin(6) * t12;
t14 = -qJ(5) - pkin(6);
t6 = t17 * pkin(4) + pkin(3);
t22 = -t12 * t14 + t13 * t6;
t1 = (t32 * t39 - m(3) * t24 + t18 * mrSges(2,1) - t16 * mrSges(2,2) - mrSges(1,3) + t38 * r_base(3) + t36 * t7 + t37 * (-t27 * t8 - t32) - t35 * (-t7 * t17 + t28 * t8) + t34 * (-t7 * qJ(3) + t24) + (m(4) * pkin(2) - m(5) * (-pkin(2) - t23) - m(6) * (-pkin(2) - t22) + t33) * t8) * g(3) + (t31 * t39 - m(3) * t25 - t16 * mrSges(2,1) - t18 * mrSges(2,2) - mrSges(1,2) + t38 * r_base(2) + t36 * t8 + t34 * (t7 * pkin(2) - t8 * qJ(3) + t25) + t37 * (t27 * t7 - t31) - t35 * (-t8 * t17 - t28 * t7) + (-m(5) * t23 - m(6) * t22 - t33) * t7) * g(2) + (-m(1) * r_base(1) - m(2) * t26 - mrSges(1,1) - mrSges(2,3) - mrSges(3,3) + (-m(3) + t34) * (qJ(2) + t26) + (m(5) * pkin(6) - m(6) * t14 + t40) * t13 + (-m(5) * pkin(3) - m(6) * t6 + t35 * t15 + t37 * t17 - mrSges(4,1)) * t12) * g(1);
U = t1;
