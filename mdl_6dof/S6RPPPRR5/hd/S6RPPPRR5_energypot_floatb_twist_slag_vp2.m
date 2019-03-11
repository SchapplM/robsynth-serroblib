% Calculate potential energy for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPPRR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPPRR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR5_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:09
% EndTime: 2019-03-09 01:37:10
% DurationCPUTime: 0.43s
% Computational Cost: add. (157->64), mult. (211->50), div. (0->0), fcn. (205->8), ass. (0->28)
t43 = m(6) + m(7);
t15 = sin(qJ(6));
t17 = cos(qJ(6));
t42 = m(7) * pkin(5) + t17 * mrSges(7,1) - t15 * mrSges(7,2) + mrSges(6,1);
t41 = m(7) * pkin(8) - mrSges(6,2) + mrSges(7,3);
t40 = -m(1) - m(2);
t38 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t37 = mrSges(2,2) - mrSges(3,3) - mrSges(4,1);
t36 = t15 * mrSges(7,1) + t17 * mrSges(7,2) + t43 * pkin(7) - mrSges(5,2) + mrSges(6,3);
t16 = sin(qJ(5));
t18 = cos(qJ(5));
t35 = t41 * t16 + t42 * t18 + mrSges(5,1);
t34 = sin(qJ(1));
t33 = sin(pkin(9));
t13 = pkin(6) + r_base(3);
t32 = t34 * pkin(1) + r_base(2);
t30 = pkin(2) + t13;
t19 = cos(qJ(1));
t29 = t19 * pkin(1) + t34 * qJ(2) + r_base(1);
t28 = t19 * qJ(3) + t29;
t26 = t34 * pkin(3) + t28;
t25 = -t19 * qJ(2) + t32;
t7 = t34 * qJ(3);
t21 = t7 + (-pkin(3) - qJ(2)) * t19 + t32;
t14 = cos(pkin(9));
t4 = t14 * t34 + t19 * t33;
t3 = t19 * t14 - t34 * t33;
t1 = (-m(1) * r_base(3) - m(4) * t30 - mrSges(3,1) + mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(5,3) + (-m(5) - t43) * (qJ(4) + t30) + t41 * t18 - t42 * t16 + (-m(2) - m(3)) * t13) * g(3) + (-mrSges(1,2) - m(3) * t25 - m(4) * (t25 + t7) - m(5) * t21 + t40 * r_base(2) - t43 * (-t3 * pkin(4) + t21) + t36 * t4 + t38 * t34 + t35 * t3 - t37 * t19) * g(2) + (-m(3) * t29 - m(4) * t28 - m(5) * t26 - mrSges(1,1) + t40 * r_base(1) - t43 * (t4 * pkin(4) + t26) - t35 * t4 + t37 * t34 + t36 * t3 + t38 * t19) * g(1);
U  = t1;
