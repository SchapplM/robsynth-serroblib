% Calculate potential energy for
% S4RRRP4
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
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRP4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RRRP4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP4_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:07
% EndTime: 2019-12-31 17:15:07
% DurationCPUTime: 0.21s
% Computational Cost: add. (92->40), mult. (91->28), div. (0->0), fcn. (61->6), ass. (0->17)
t22 = -m(4) - m(5);
t25 = mrSges(4,2) + mrSges(5,2);
t24 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1);
t23 = -m(2) - m(3);
t10 = sin(qJ(2));
t12 = cos(qJ(2));
t9 = qJ(2) + qJ(3);
t3 = sin(t9);
t4 = cos(t9);
t21 = -m(3) * pkin(1) - t12 * mrSges(3,1) + t10 * mrSges(3,2) - mrSges(2,1) + t22 * (t12 * pkin(2) + pkin(1)) + t24 * t4 + t25 * t3;
t20 = -m(1) + t22 + t23;
t14 = -pkin(6) - pkin(5);
t19 = m(3) * pkin(5) - m(4) * t14 - m(5) * (-qJ(4) + t14) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t8 = pkin(4) + r_base(3);
t13 = cos(qJ(1));
t11 = sin(qJ(1));
t1 = (-m(1) * r_base(3) - mrSges(3,1) * t10 - mrSges(3,2) * t12 - mrSges(1,3) - mrSges(2,3) + t23 * t8 - t25 * t4 + t22 * (t10 * pkin(2) + t8) + t24 * t3) * g(3) + (t21 * t11 + t19 * t13 + t20 * r_base(2) - mrSges(1,2)) * g(2) + (-t19 * t11 + t21 * t13 + t20 * r_base(1) - mrSges(1,1)) * g(1);
U = t1;
