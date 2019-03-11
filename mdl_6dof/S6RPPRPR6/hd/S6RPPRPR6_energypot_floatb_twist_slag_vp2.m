% Calculate potential energy for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
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
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRPR6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRPR6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR6_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR6_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:50:42
% EndTime: 2019-03-09 01:50:42
% DurationCPUTime: 0.37s
% Computational Cost: add. (131->71), mult. (164->49), div. (0->0), fcn. (126->6), ass. (0->28)
t40 = m(6) + m(7);
t39 = -m(7) * pkin(8) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t12 = sin(qJ(6));
t15 = cos(qJ(6));
t38 = -t12 * mrSges(7,1) - t15 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t37 = -m(1) - m(2);
t13 = sin(qJ(4));
t16 = cos(qJ(4));
t35 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3) + (t40 * qJ(5) - t38) * t16 + t39 * t13;
t34 = -t15 * mrSges(7,1) + t12 * mrSges(7,2) - mrSges(6,1) - mrSges(2,2) + mrSges(4,2) + mrSges(3,3) - mrSges(5,3);
t33 = pkin(4) * t13;
t11 = pkin(6) + r_base(3);
t14 = sin(qJ(1));
t32 = t14 * pkin(1) + r_base(2);
t31 = pkin(2) + t11;
t17 = cos(qJ(1));
t30 = t17 * pkin(1) + t14 * qJ(2) + r_base(1);
t29 = pkin(3) + t31;
t28 = t17 * qJ(3) + t30;
t24 = -t17 * qJ(2) + t32;
t4 = t14 * qJ(3);
t22 = t24 + t4;
t21 = -t14 * pkin(7) + t28;
t9 = t17 * pkin(7);
t20 = t22 + t9;
t2 = t17 * t33;
t1 = t14 * t33;
t3 = (-m(1) * r_base(3) - m(4) * t31 - m(5) * t29 - mrSges(3,1) - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) - t40 * (t16 * pkin(4) + t13 * qJ(5) + t29) + (-m(2) - m(3)) * t11 + t39 * t16 + t38 * t13) * g(3) + (-mrSges(1,2) - m(3) * t24 - m(4) * t22 - m(5) * t20 - m(6) * (t1 + t20) - m(7) * (t1 + t4 + t9 + t32) + t37 * r_base(2) + (-m(7) * (pkin(5) - qJ(2)) + t34) * t17 + t35 * t14) * g(2) + (-mrSges(1,1) - m(3) * t30 - m(4) * t28 - m(5) * t21 - m(6) * (t2 + t21) - m(7) * (t2 + t28) + t37 * r_base(1) + t35 * t17 + (-m(7) * (-pkin(5) - pkin(7)) - t34) * t14) * g(1);
U  = t3;
