% Calculate potential energy for
% S4PPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PPRR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4PPRR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR1_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR1_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:15:27
% EndTime: 2019-03-08 18:15:27
% DurationCPUTime: 0.16s
% Computational Cost: add. (80->45), mult. (84->39), div. (0->0), fcn. (62->6), ass. (0->18)
t22 = m(3) + m(4) + m(5);
t21 = qJ(1) + r_base(3);
t20 = -pkin(4) + t21;
t19 = -m(1) - m(2) - t22;
t15 = sin(qJ(3));
t18 = m(5) * pkin(3) * t15 - mrSges(2,2) + mrSges(3,3);
t16 = cos(qJ(3));
t17 = -m(4) * pkin(2) - m(5) * (pkin(3) * t16 + pkin(2)) - mrSges(2,1) - mrSges(3,1);
t14 = cos(pkin(6));
t13 = sin(pkin(6));
t12 = qJ(3) + qJ(4);
t8 = cos(t12);
t7 = sin(t12);
t4 = t13 * t16 - t14 * t15;
t3 = -t13 * t15 - t14 * t16;
t2 = t13 * t8 - t14 * t7;
t1 = -t13 * t7 - t14 * t8;
t5 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - mrSges(3,2) - m(4) * t20 + mrSges(4,3) - m(5) * (-pkin(5) + t20) + mrSges(5,3) + (-m(2) - m(3)) * t21) * g(3) + (-t4 * mrSges(4,1) - t2 * mrSges(5,1) - t3 * mrSges(4,2) - t1 * mrSges(5,2) - mrSges(1,2) + t19 * r_base(2) + (t22 * qJ(2) + t18) * t14 + (-t22 * pkin(1) + t17) * t13) * g(2) + (t3 * mrSges(4,1) + t1 * mrSges(5,1) - t4 * mrSges(4,2) - t2 * mrSges(5,2) - t18 * t13 + t17 * t14 + t19 * r_base(1) - mrSges(1,1) - t22 * (t14 * pkin(1) + t13 * qJ(2))) * g(1);
U  = t5;
