% Calculate potential energy for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RRRR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR1_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR1_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:06
% EndTime: 2019-12-31 17:22:06
% DurationCPUTime: 0.18s
% Computational Cost: add. (99->45), mult. (68->35), div. (0->0), fcn. (38->8), ass. (0->19)
t26 = -m(1) - m(2);
t25 = -m(4) - m(5);
t12 = sin(qJ(4));
t14 = cos(qJ(4));
t24 = -m(5) * pkin(3) - t14 * mrSges(5,1) + t12 * mrSges(5,2) - mrSges(4,1);
t23 = m(5) * pkin(7) - mrSges(4,2) + mrSges(5,3);
t11 = qJ(1) + qJ(2);
t22 = pkin(4) + r_base(3);
t13 = sin(qJ(1));
t21 = t13 * pkin(1) + r_base(2);
t15 = cos(qJ(1));
t20 = t15 * pkin(1) + r_base(1);
t19 = pkin(5) + t22;
t8 = qJ(3) + t11;
t7 = cos(t11);
t6 = sin(t11);
t4 = cos(t8);
t3 = sin(t8);
t1 = (-m(1) * r_base(3) - m(2) * t22 - m(3) * t19 - mrSges(5,1) * t12 - mrSges(5,2) * t14 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) + t25 * (pkin(6) + t19)) * g(3) + (-m(3) * t21 - t13 * mrSges(2,1) - t6 * mrSges(3,1) - mrSges(2,2) * t15 - t7 * mrSges(3,2) - mrSges(1,2) + t26 * r_base(2) + t25 * (pkin(2) * t6 + t21) + t23 * t4 + t24 * t3) * g(2) + (-m(3) * t20 - mrSges(2,1) * t15 - t7 * mrSges(3,1) + t13 * mrSges(2,2) + t6 * mrSges(3,2) - mrSges(1,1) + t26 * r_base(1) + t24 * t4 + t25 * (pkin(2) * t7 + t20) - t23 * t3) * g(1);
U = t1;
