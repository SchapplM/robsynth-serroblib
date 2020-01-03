% Calculate potential energy for
% S4RPPR4
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPPR4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RPPR4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR4_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR4_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:47
% EndTime: 2019-12-31 16:38:47
% DurationCPUTime: 0.19s
% Computational Cost: add. (88->43), mult. (70->32), div. (0->0), fcn. (40->6), ass. (0->15)
t23 = -m(1) - m(2);
t22 = m(4) + m(5);
t10 = sin(qJ(4));
t12 = cos(qJ(4));
t21 = -t10 * mrSges(5,1) - t12 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t20 = -m(5) * pkin(5) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t19 = pkin(4) + r_base(3);
t11 = sin(qJ(1));
t18 = t11 * pkin(1) + r_base(2);
t13 = cos(qJ(1));
t17 = t13 * pkin(1) + r_base(1);
t9 = qJ(1) + pkin(6);
t5 = cos(t9);
t4 = sin(t9);
t1 = (-m(1) * r_base(3) - m(2) * t19 - m(5) * pkin(3) - t12 * mrSges(5,1) + t10 * mrSges(5,2) - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + (-m(3) - t22) * (qJ(2) + t19)) * g(3) + (-m(3) * t18 - t11 * mrSges(2,1) - mrSges(2,2) * t13 - mrSges(1,2) + t23 * r_base(2) - t22 * (t4 * pkin(2) + t18) + (t22 * qJ(3) - t21) * t5 + t20 * t4) * g(2) + (-m(3) * t17 - mrSges(2,1) * t13 + t11 * mrSges(2,2) - mrSges(1,1) + t23 * r_base(1) - t22 * (t5 * pkin(2) + t4 * qJ(3) + t17) + t20 * t5 + t21 * t4) * g(1);
U = t1;
