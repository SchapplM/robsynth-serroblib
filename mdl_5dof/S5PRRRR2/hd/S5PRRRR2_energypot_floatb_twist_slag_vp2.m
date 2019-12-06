% Calculate potential energy for
% S5PRRRR2
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
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR2_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:40
% EndTime: 2019-12-05 17:04:41
% DurationCPUTime: 0.20s
% Computational Cost: add. (104->51), mult. (72->35), div. (0->0), fcn. (36->8), ass. (0->20)
t27 = -m(2) - m(3);
t26 = -m(5) - m(6);
t13 = sin(qJ(5));
t15 = cos(qJ(5));
t25 = -t15 * mrSges(6,1) + t13 * mrSges(6,2) - mrSges(5,1);
t24 = m(6) * pkin(6) - mrSges(5,2) + mrSges(6,3);
t12 = qJ(2) + qJ(3);
t23 = pkin(1) + r_base(1);
t14 = sin(qJ(2));
t22 = t14 * pkin(2) + r_base(2);
t11 = qJ(1) + r_base(3);
t16 = cos(qJ(2));
t20 = t16 * pkin(2) + t23;
t19 = pkin(4) + t11;
t8 = qJ(4) + t12;
t7 = cos(t12);
t6 = sin(t12);
t4 = cos(t8);
t3 = sin(t8);
t1 = (-m(1) * r_base(3) - m(4) * t19 - mrSges(6,1) * t13 - mrSges(6,2) * t15 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) + t26 * (pkin(5) + t19) + t27 * t11) * g(3) + (-m(4) * t22 - t14 * mrSges(3,1) - t6 * mrSges(4,1) - mrSges(3,2) * t16 - t7 * mrSges(4,2) - mrSges(1,2) - mrSges(2,2) + t26 * (pkin(3) * t6 + t22) + t24 * t4 + t25 * t3 + (-m(1) + t27) * r_base(2)) * g(2) + (-m(3) * t23 - m(4) * t20 - t16 * mrSges(3,1) - t7 * mrSges(4,1) + t14 * mrSges(3,2) + t6 * mrSges(4,2) - mrSges(1,1) - mrSges(2,1) + (-m(1) - m(2)) * r_base(1) + t25 * t4 + t26 * (pkin(3) * t7 + t20) - t24 * t3) * g(1);
U = t1;
