% Calculate potential energy for
% S6RPPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
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
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPPRR4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPPRR4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:14
% EndTime: 2019-03-09 01:35:14
% DurationCPUTime: 0.44s
% Computational Cost: add. (170->62), mult. (237->47), div. (0->0), fcn. (243->8), ass. (0->26)
t44 = m(6) + m(7);
t45 = m(5) + t44;
t15 = sin(qJ(6));
t17 = cos(qJ(6));
t43 = m(7) * pkin(5) + t17 * mrSges(7,1) - t15 * mrSges(7,2) + mrSges(6,1);
t42 = m(7) * pkin(8) - mrSges(6,2) + mrSges(7,3);
t41 = -m(1) - m(2);
t39 = -mrSges(2,1) - mrSges(3,1);
t38 = mrSges(2,2) - mrSges(3,3);
t36 = t15 * mrSges(7,1) + t17 * mrSges(7,2) + t44 * pkin(7) + mrSges(4,1) - mrSges(5,2) + mrSges(6,3);
t16 = sin(qJ(5));
t18 = cos(qJ(5));
t35 = t45 * qJ(4) + t43 * t16 - t42 * t18 - mrSges(4,2) + mrSges(5,3);
t34 = cos(qJ(1));
t33 = sin(qJ(1));
t31 = cos(pkin(9));
t30 = sin(pkin(9));
t14 = pkin(6) + r_base(3);
t29 = t34 * pkin(1) + t33 * qJ(2) + r_base(1);
t8 = -qJ(3) + t14;
t26 = t34 * pkin(2) + t29;
t24 = t33 * pkin(1) - t34 * qJ(2) + r_base(2);
t20 = t33 * pkin(2) + t24;
t4 = t34 * t30 - t33 * t31;
t3 = -t33 * t30 - t34 * t31;
t1 = (-m(1) * r_base(3) + mrSges(5,1) - mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) + (-m(4) - m(5)) * t8 - t44 * (-pkin(4) + t8) + t43 * t18 + t42 * t16 + (-m(2) - m(3)) * t14) * g(3) + (-m(3) * t24 - m(4) * t20 - mrSges(1,2) + t41 * r_base(2) - t38 * t34 + t39 * t33 - t45 * (-t4 * pkin(3) + t20) + t36 * t4 + t35 * t3) * g(2) + (-m(3) * t29 - m(4) * t26 - mrSges(1,1) + t41 * r_base(1) + t39 * t34 + t38 * t33 - t45 * (-t3 * pkin(3) + t26) - t35 * t4 + t36 * t3) * g(1);
U  = t1;
