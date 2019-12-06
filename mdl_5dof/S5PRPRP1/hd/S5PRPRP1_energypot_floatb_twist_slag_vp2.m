% Calculate potential energy for
% S5PRPRP1
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRPRP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP1_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP1_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:14
% EndTime: 2019-12-05 15:28:14
% DurationCPUTime: 0.31s
% Computational Cost: add. (155->55), mult. (116->41), div. (0->0), fcn. (80->8), ass. (0->24)
t35 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t34 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t33 = -m(1) - m(2);
t32 = -m(3) - m(4);
t31 = -m(5) - m(6);
t15 = sin(pkin(8));
t17 = cos(pkin(8));
t13 = pkin(8) + qJ(4);
t5 = sin(t13);
t7 = cos(t13);
t30 = -m(4) * pkin(2) - t17 * mrSges(4,1) + t15 * mrSges(4,2) - t34 * t5 + t35 * t7 - mrSges(3,1);
t29 = m(4) * qJ(3) - mrSges(3,2) + mrSges(6,2) + mrSges(4,3) + mrSges(5,3);
t16 = sin(pkin(7));
t28 = t16 * pkin(1) + r_base(2);
t18 = cos(pkin(7));
t27 = t18 * pkin(1) + r_base(1);
t26 = qJ(1) + r_base(3);
t9 = pkin(5) + t26;
t19 = -pkin(6) - qJ(3);
t14 = pkin(7) + qJ(2);
t8 = cos(t14);
t6 = sin(t14);
t4 = t17 * pkin(3) + pkin(2);
t1 = (-m(1) * r_base(3) - m(2) * t26 - t15 * mrSges(4,1) - t17 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t32 * t9 + t31 * (t15 * pkin(3) + t9) + t34 * t7 + t35 * t5) * g(3) + (-t16 * mrSges(2,1) - t18 * mrSges(2,2) - mrSges(1,2) + t33 * r_base(2) + t32 * t28 + t31 * (t8 * t19 + t6 * t4 + t28) + t29 * t8 + t30 * t6) * g(2) + (-t18 * mrSges(2,1) + t16 * mrSges(2,2) - mrSges(1,1) + t33 * r_base(1) + t32 * t27 + t31 * (-t6 * t19 + t8 * t4 + t27) + t30 * t8 - t29 * t6) * g(1);
U = t1;
