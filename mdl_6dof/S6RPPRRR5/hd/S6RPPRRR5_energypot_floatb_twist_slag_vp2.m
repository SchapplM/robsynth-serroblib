% Calculate potential energy for
% S6RPPRRR5
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
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR5_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR5_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:28:11
% EndTime: 2019-03-09 02:28:11
% DurationCPUTime: 0.38s
% Computational Cost: add. (151->69), mult. (157->49), div. (0->0), fcn. (119->8), ass. (0->27)
t14 = sin(qJ(6));
t17 = cos(qJ(6));
t44 = -m(7) * pkin(5) - t17 * mrSges(7,1) + t14 * mrSges(7,2) - mrSges(6,1);
t43 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t42 = -m(1) - m(2);
t41 = -m(6) - m(7);
t15 = sin(qJ(4));
t18 = cos(qJ(4));
t13 = qJ(4) + qJ(5);
t4 = sin(t13);
t5 = cos(t13);
t39 = -t15 * mrSges(5,1) - t18 * mrSges(5,2) + t44 * t4 - t43 * t5 - mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t38 = -t14 * mrSges(7,1) - t17 * mrSges(7,2) - mrSges(2,2) + mrSges(4,2) + mrSges(3,3) - mrSges(5,3) - mrSges(6,3);
t37 = pkin(4) * t15;
t12 = pkin(6) + r_base(3);
t16 = sin(qJ(1));
t35 = t16 * pkin(1) + r_base(2);
t34 = pkin(2) + t12;
t6 = t16 * qJ(3);
t33 = t6 + t35;
t19 = cos(qJ(1));
t32 = t19 * pkin(1) + t16 * qJ(2) + r_base(1);
t31 = pkin(3) + t34;
t29 = t19 * qJ(3) + t32;
t24 = -t19 * qJ(2) + t35;
t20 = -pkin(8) - pkin(7);
t1 = (-m(1) * r_base(3) - m(4) * t34 - m(5) * t31 - t18 * mrSges(5,1) + t15 * mrSges(5,2) - mrSges(3,1) - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) + t44 * t5 + t41 * (t18 * pkin(4) + t31) + t43 * t4 + (-m(2) - m(3)) * t12) * g(3) + (-mrSges(1,2) - m(3) * t24 - m(4) * (t24 + t6) - m(5) * t33 + t42 * r_base(2) + t41 * (t16 * t37 + t33) + (-m(5) * (pkin(7) - qJ(2)) + t41 * (-qJ(2) - t20) + t38) * t19 + t39 * t16) * g(2) + (-m(3) * t32 - mrSges(1,1) + t42 * r_base(1) + (-m(4) - m(5)) * t29 + t41 * (t16 * t20 + t19 * t37 + t29) + t39 * t19 + (m(5) * pkin(7) - t38) * t16) * g(1);
U  = t1;
