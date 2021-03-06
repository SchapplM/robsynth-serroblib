% Calculate potential energy for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPRP6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RPRP6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_energypot_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP6_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP6_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:56
% EndTime: 2019-12-31 16:45:57
% DurationCPUTime: 0.22s
% Computational Cost: add. (70->39), mult. (83->25), div. (0->0), fcn. (53->4), ass. (0->13)
t22 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1);
t21 = mrSges(4,2) + mrSges(5,2);
t20 = -m(1) - m(2);
t19 = -m(4) - m(5);
t18 = m(3) - t19;
t16 = -m(4) * pkin(5) + m(5) * (-qJ(4) - pkin(5)) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t6 = sin(qJ(3));
t8 = cos(qJ(3));
t15 = -t21 * t8 + t22 * t6 + mrSges(2,2) - mrSges(3,3);
t4 = pkin(4) + r_base(3);
t9 = cos(qJ(1));
t7 = sin(qJ(1));
t1 = (-m(1) * r_base(3) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + t19 * (pkin(2) + t4) + t22 * t8 + t21 * t6 + (-m(2) - m(3)) * t4) * g(3) + (-mrSges(1,2) + t20 * r_base(2) - t18 * (t7 * pkin(1) + r_base(2)) + (t18 * qJ(2) - t15) * t9 + t16 * t7) * g(2) + (-mrSges(1,1) + t20 * r_base(1) - t18 * (t9 * pkin(1) + t7 * qJ(2) + r_base(1)) + t16 * t9 + t15 * t7) * g(1);
U = t1;
