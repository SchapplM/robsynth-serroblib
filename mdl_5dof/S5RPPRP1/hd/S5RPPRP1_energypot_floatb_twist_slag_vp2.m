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
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:35:31
% EndTime: 2019-12-05 17:35:31
% DurationCPUTime: 0.44s
% Computational Cost: add. (165->67), mult. (153->62), div. (0->0), fcn. (125->8), ass. (0->29)
t42 = -mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t41 = m(6) * pkin(4);
t40 = -m(1) - m(2);
t39 = -mrSges(5,1) - mrSges(6,1);
t38 = mrSges(3,2) - mrSges(4,3);
t37 = mrSges(5,2) + mrSges(6,2);
t36 = -m(4) - m(5) - m(6);
t14 = sin(pkin(8));
t15 = cos(pkin(8));
t35 = t15 * mrSges(4,1) + t42 * t14 + mrSges(3,1);
t17 = sin(qJ(4));
t13 = qJ(1) + pkin(7);
t9 = sin(t13);
t34 = t9 * t17;
t10 = cos(t13);
t33 = t10 * t17;
t30 = t15 * t17;
t19 = cos(qJ(4));
t29 = t15 * t19;
t28 = pkin(5) + r_base(1);
t20 = cos(qJ(1));
t27 = t20 * pkin(1) + r_base(3);
t18 = sin(qJ(1));
t26 = -t18 * pkin(1) + r_base(2);
t24 = pkin(3) * t15 + pkin(6) * t14;
t16 = -qJ(5) - pkin(6);
t8 = t19 * pkin(4) + pkin(3);
t23 = -t14 * t16 + t15 * t8;
t1 = (-t34 * t41 - m(3) * t27 - t20 * mrSges(2,1) + t18 * mrSges(2,2) - mrSges(1,3) + t40 * r_base(3) + t38 * t9 + t39 * (t10 * t29 + t34) - t37 * (-t10 * t30 + t9 * t19) + t36 * (t10 * pkin(2) + t9 * qJ(3) + t27) + (-m(5) * t24 - m(6) * t23 - t35) * t10) * g(3) + (-t33 * t41 - m(3) * t26 + t18 * mrSges(2,1) + t20 * mrSges(2,2) - mrSges(1,2) + t40 * r_base(2) + t36 * (t10 * qJ(3) + t26) + t39 * (-t29 * t9 + t33) + t38 * t10 - t37 * (t10 * t19 + t30 * t9) + (m(4) * pkin(2) - m(5) * (-pkin(2) - t24) - m(6) * (-pkin(2) - t23) + t35) * t9) * g(2) + (-m(1) * r_base(1) - m(2) * t28 - mrSges(1,1) - mrSges(2,3) - mrSges(3,3) + (-m(3) + t36) * (qJ(2) + t28) + (m(5) * pkin(6) - m(6) * t16 + t42) * t15 + (-m(5) * pkin(3) - m(6) * t8 + t37 * t17 + t39 * t19 - mrSges(4,1)) * t14) * g(1);
U = t1;
