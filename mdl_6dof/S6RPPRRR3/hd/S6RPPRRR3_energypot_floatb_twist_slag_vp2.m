% Calculate potential energy for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:22:43
% EndTime: 2019-03-09 02:22:44
% DurationCPUTime: 0.42s
% Computational Cost: add. (214->66), mult. (175->49), div. (0->0), fcn. (141->10), ass. (0->25)
t15 = qJ(5) + qJ(6);
t10 = sin(t15);
t11 = cos(t15);
t16 = sin(qJ(5));
t19 = cos(qJ(5));
t44 = mrSges(5,1) + m(6) * pkin(4) + t19 * mrSges(6,1) - t16 * mrSges(6,2) + m(7) * (t19 * pkin(5) + pkin(4)) + t11 * mrSges(7,1) - t10 * mrSges(7,2);
t43 = mrSges(5,2) - m(6) * pkin(8) + m(7) * (-pkin(9) - pkin(8)) - mrSges(6,3) - mrSges(7,3);
t17 = sin(qJ(4));
t20 = cos(qJ(4));
t42 = t44 * t17 + t43 * t20 - mrSges(3,2) + mrSges(4,3);
t41 = -m(1) - m(2);
t39 = -m(5) - m(6) - m(7);
t37 = -t10 * mrSges(7,1) - t19 * mrSges(6,2) - t11 * mrSges(7,2) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t16;
t35 = pkin(6) + r_base(3);
t18 = sin(qJ(1));
t34 = t18 * pkin(1) + r_base(2);
t21 = cos(qJ(1));
t33 = t21 * pkin(1) + r_base(1);
t14 = qJ(1) + pkin(10);
t7 = sin(t14);
t32 = t7 * pkin(2) + t34;
t9 = qJ(2) + t35;
t8 = cos(t14);
t30 = t8 * pkin(2) + t7 * qJ(3) + t33;
t1 = (-m(1) * r_base(3) - m(2) * t35 - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + (-m(3) - m(4)) * t9 + t39 * (pkin(3) + t9) - t44 * t20 + t43 * t17) * g(3) + (-m(3) * t34 - m(4) * t32 - t18 * mrSges(2,1) - t21 * mrSges(2,2) - mrSges(1,2) + t41 * r_base(2) + t39 * (t7 * pkin(7) + t32) + ((m(4) - t39) * qJ(3) + t42) * t8 + t37 * t7) * g(2) + (-m(3) * t33 - m(4) * t30 - t21 * mrSges(2,1) + t18 * mrSges(2,2) - mrSges(1,1) + t41 * r_base(1) + t39 * (t8 * pkin(7) + t30) + t37 * t8 - t42 * t7) * g(1);
U  = t1;
