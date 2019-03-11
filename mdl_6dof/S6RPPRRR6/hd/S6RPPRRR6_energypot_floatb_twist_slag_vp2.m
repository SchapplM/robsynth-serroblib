% Calculate potential energy for
% S6RPPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRR6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRR6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR6_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR6_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:30:31
% EndTime: 2019-03-09 02:30:32
% DurationCPUTime: 0.36s
% Computational Cost: add. (149->69), mult. (177->48), div. (0->0), fcn. (143->8), ass. (0->26)
t12 = sin(qJ(5));
t15 = cos(qJ(5));
t11 = qJ(5) + qJ(6);
t2 = sin(t11);
t3 = cos(t11);
t43 = -m(6) * pkin(4) - t15 * mrSges(6,1) + t12 * mrSges(6,2) - mrSges(5,1) - m(7) * (pkin(5) * t15 + pkin(4)) - t3 * mrSges(7,1) + t2 * mrSges(7,2);
t42 = mrSges(5,2) + m(7) * (-pkin(9) - pkin(8)) - mrSges(7,3) - m(6) * pkin(8) - mrSges(6,3);
t41 = -m(1) - m(2);
t40 = -m(6) - m(7);
t39 = m(5) - t40;
t13 = sin(qJ(4));
t16 = cos(qJ(4));
t36 = t43 * t13 - t42 * t16 - mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t35 = -t12 * mrSges(6,1) - t2 * mrSges(7,1) - t15 * mrSges(6,2) - t3 * mrSges(7,2) - mrSges(2,2) + mrSges(4,2) + mrSges(3,3) - mrSges(5,3);
t34 = pkin(5) * t12;
t10 = pkin(6) + r_base(3);
t14 = sin(qJ(1));
t33 = t14 * pkin(1) + r_base(2);
t32 = pkin(2) + t10;
t17 = cos(qJ(1));
t31 = t17 * pkin(1) + t14 * qJ(2) + r_base(1);
t24 = -t17 * qJ(2) + t33;
t4 = t14 * qJ(3);
t23 = t24 + t4;
t8 = t17 * pkin(7);
t1 = (-m(1) * r_base(3) - m(4) * t32 - mrSges(3,1) - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t10 - t39 * (pkin(3) + t32) + t43 * t16 + t42 * t13) * g(3) + (-mrSges(1,2) - m(3) * t24 - m(4) * t23 - m(5) * (t23 + t8) + t41 * r_base(2) + t40 * (t4 + t8 + t33) + (m(6) * qJ(2) - m(7) * (-qJ(2) + t34) + t35) * t17 + t36 * t14) * g(2) + (-m(3) * t31 - mrSges(1,1) + t41 * r_base(1) + (-m(4) - t39) * (t17 * qJ(3) + t31) + t36 * t17 + (m(7) * t34 + t39 * pkin(7) - t35) * t14) * g(1);
U  = t1;
