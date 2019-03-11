% Calculate potential energy for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
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
% Datum: 2019-03-08 18:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR1_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:36:51
% EndTime: 2019-03-08 18:36:52
% DurationCPUTime: 0.28s
% Computational Cost: add. (150->51), mult. (134->41), div. (0->0), fcn. (102->10), ass. (0->24)
t13 = sin(qJ(5));
t16 = cos(qJ(5));
t36 = m(6) * pkin(4) + t16 * mrSges(6,1) - t13 * mrSges(6,2) + mrSges(5,1);
t35 = -m(6) * pkin(6) + mrSges(5,2) - mrSges(6,3);
t34 = -m(2) - m(3);
t33 = -m(5) - m(6);
t32 = -m(1) - m(4) + t34;
t14 = sin(qJ(2));
t17 = cos(qJ(2));
t12 = qJ(2) + qJ(3);
t9 = qJ(4) + t12;
t4 = sin(t9);
t5 = cos(t9);
t6 = t17 * pkin(2) + pkin(1);
t7 = sin(t12);
t8 = cos(t12);
t30 = -m(3) * pkin(1) - m(4) * t6 - t17 * mrSges(3,1) - t8 * mrSges(4,1) + t14 * mrSges(3,2) + t7 * mrSges(4,2) + t35 * t4 - t36 * t5 - mrSges(2,1);
t29 = t13 * mrSges(6,1) + t16 * mrSges(6,2) + mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t11 = pkin(5) + r_base(3);
t25 = -pkin(2) * t14 + t11;
t18 = cos(qJ(1));
t15 = sin(qJ(1));
t3 = pkin(3) * t8 + t6;
t1 = (-m(1) * r_base(3) - m(4) * t25 + mrSges(3,1) * t14 + t7 * mrSges(4,1) + mrSges(3,2) * t17 + t8 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) + t33 * (-pkin(3) * t7 + t25) + t35 * t5 + t36 * t4 + t34 * t11) * g(3) + (-mrSges(1,2) + t33 * (t15 * t3 + r_base(2)) + t32 * r_base(2) - t29 * t18 + t30 * t15) * g(2) + (-mrSges(1,1) + t33 * (t18 * t3 + r_base(1)) + t32 * r_base(1) + t30 * t18 + t29 * t15) * g(1);
U  = t1;
