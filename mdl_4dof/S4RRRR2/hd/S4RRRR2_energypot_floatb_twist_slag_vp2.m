% Calculate potential energy for
% S4RRRR2
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
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RRRR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR2_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR2_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:05
% EndTime: 2019-12-31 17:23:06
% DurationCPUTime: 0.19s
% Computational Cost: add. (99->44), mult. (79->34), div. (0->0), fcn. (49->8), ass. (0->16)
t24 = -m(1) - m(2);
t23 = -m(3) - m(4) - m(5);
t11 = sin(qJ(3));
t13 = cos(qJ(3));
t9 = qJ(3) + qJ(4);
t2 = sin(t9);
t4 = cos(t9);
t22 = -m(4) * pkin(2) - t13 * mrSges(4,1) + t11 * mrSges(4,2) - mrSges(3,1) - m(5) * (pkin(3) * t13 + pkin(2)) - t4 * mrSges(5,1) + t2 * mrSges(5,2);
t21 = m(4) * pkin(6) - m(5) * (-pkin(7) - pkin(6)) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t20 = pkin(4) + r_base(3);
t14 = cos(qJ(1));
t12 = sin(qJ(1));
t10 = qJ(1) + qJ(2);
t5 = cos(t10);
t3 = sin(t10);
t1 = (-m(1) * r_base(3) - m(2) * t20 - t2 * mrSges(5,1) - mrSges(4,2) * t13 - t4 * mrSges(5,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + (-m(5) * pkin(3) - mrSges(4,1)) * t11 + t23 * (pkin(5) + t20)) * g(3) + (-mrSges(2,1) * t12 - mrSges(2,2) * t14 - mrSges(1,2) + t24 * r_base(2) + t23 * (t12 * pkin(1) + r_base(2)) + t21 * t5 + t22 * t3) * g(2) + (-mrSges(2,1) * t14 + mrSges(2,2) * t12 - mrSges(1,1) + t24 * r_base(1) + t22 * t5 + t23 * (t14 * pkin(1) + r_base(1)) - t21 * t3) * g(1);
U = t1;
