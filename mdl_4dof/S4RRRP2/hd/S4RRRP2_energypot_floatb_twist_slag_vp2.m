% Calculate potential energy for
% S4RRRP2
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRP2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RRRP2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP2_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP2_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:50
% EndTime: 2019-12-31 17:12:51
% DurationCPUTime: 0.20s
% Computational Cost: add. (93->40), mult. (79->28), div. (0->0), fcn. (49->6), ass. (0->16)
t24 = -m(4) - m(5);
t23 = mrSges(4,2) + mrSges(5,2);
t22 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1);
t21 = -m(1) - m(2);
t20 = -m(3) + t24;
t11 = cos(qJ(3));
t9 = sin(qJ(3));
t19 = t24 * pkin(2) + t22 * t11 + t23 * t9 - mrSges(3,1);
t18 = m(4) * pkin(6) - m(5) * (-qJ(4) - pkin(6)) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t17 = pkin(4) + r_base(3);
t12 = cos(qJ(1));
t10 = sin(qJ(1));
t7 = qJ(1) + qJ(2);
t3 = cos(t7);
t2 = sin(t7);
t1 = (-m(1) * r_base(3) - m(2) * t17 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t22 * t9 + t20 * (pkin(5) + t17) - t23 * t11) * g(3) + (-t10 * mrSges(2,1) - mrSges(2,2) * t12 - mrSges(1,2) + t21 * r_base(2) + t20 * (t10 * pkin(1) + r_base(2)) + t18 * t3 + t19 * t2) * g(2) + (-mrSges(2,1) * t12 + t10 * mrSges(2,2) - mrSges(1,1) + t21 * r_base(1) + t19 * t3 + t20 * (t12 * pkin(1) + r_base(1)) - t18 * t2) * g(1);
U = t1;
