% Calculate potential energy for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRP6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPRP6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP6_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:54:54
% EndTime: 2019-12-31 17:54:54
% DurationCPUTime: 0.35s
% Computational Cost: add. (118->58), mult. (122->42), div. (0->0), fcn. (86->6), ass. (0->23)
t33 = mrSges(5,1) + mrSges(6,1);
t32 = mrSges(5,2) - mrSges(6,3);
t31 = -m(1) - m(2);
t30 = m(3) + m(4);
t29 = -m(5) - m(6);
t10 = sin(pkin(7));
t11 = cos(pkin(7));
t8 = pkin(7) + qJ(4);
t2 = sin(t8);
t3 = cos(t8);
t28 = t10 * mrSges(4,1) + t11 * mrSges(4,2) + t33 * t2 + t32 * t3 - mrSges(2,2) + mrSges(3,3);
t27 = -m(4) * qJ(3) - mrSges(2,1) + mrSges(3,2) - mrSges(6,2) - mrSges(4,3) - mrSges(5,3);
t26 = pkin(3) * t10;
t9 = pkin(5) + r_base(3);
t13 = sin(qJ(1));
t25 = t13 * pkin(1) + r_base(2);
t24 = pkin(2) + t9;
t14 = cos(qJ(1));
t23 = t14 * pkin(1) + t13 * qJ(2) + r_base(1);
t22 = -qJ(2) - t26;
t18 = pkin(4) * t2 - qJ(5) * t3;
t12 = -pkin(6) - qJ(3);
t1 = (-m(1) * r_base(3) - m(4) * t24 - t11 * mrSges(4,1) + t10 * mrSges(4,2) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t9 + t29 * (t11 * pkin(3) + t24) + (-m(6) * pkin(4) - t33) * t3 + (-m(6) * qJ(5) + t32) * t2) * g(3) + (-mrSges(1,2) + t31 * r_base(2) - t30 * t25 + t29 * (-t13 * t12 + t25) + (-m(5) * t22 - m(6) * (-t18 + t22) + t30 * qJ(2) + t28) * t14 + t27 * t13) * g(2) + (-mrSges(1,1) + t31 * r_base(1) - t30 * t23 + t29 * (-t12 * t14 + t13 * t26 + t23) + t27 * t14 + (-m(6) * t18 - t28) * t13) * g(1);
U = t1;
