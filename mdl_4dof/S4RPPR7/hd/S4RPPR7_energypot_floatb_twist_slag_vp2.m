% Calculate potential energy for
% S4RPPR7
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPPR7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RPPR7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_energypot_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR7_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:33
% EndTime: 2019-12-31 16:41:34
% DurationCPUTime: 0.21s
% Computational Cost: add. (76->43), mult. (83->31), div. (0->0), fcn. (53->6), ass. (0->14)
t22 = -m(1) - m(2);
t21 = m(3) + m(4) + m(5);
t6 = pkin(6) + qJ(4);
t1 = sin(t6);
t2 = cos(t6);
t8 = sin(pkin(6));
t9 = cos(pkin(6));
t20 = t1 * mrSges(5,1) + t9 * mrSges(4,2) + t2 * mrSges(5,2) - mrSges(2,2) + mrSges(3,3) + (m(5) * pkin(3) + mrSges(4,1)) * t8;
t19 = -m(4) * qJ(3) + m(5) * (-pkin(5) - qJ(3)) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t7 = pkin(4) + r_base(3);
t16 = pkin(2) + t7;
t12 = cos(qJ(1));
t11 = sin(qJ(1));
t3 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - mrSges(3,1) - m(4) * t16 - t9 * mrSges(4,1) + t8 * mrSges(4,2) - m(5) * (pkin(3) * t9 + t16) - t2 * mrSges(5,1) + t1 * mrSges(5,2) + (-m(2) - m(3)) * t7) * g(3) + (-mrSges(1,2) + t22 * r_base(2) - t21 * (t11 * pkin(1) + r_base(2)) + (t21 * qJ(2) + t20) * t12 + t19 * t11) * g(2) + (-mrSges(1,1) + t22 * r_base(1) - t21 * (t12 * pkin(1) + t11 * qJ(2) + r_base(1)) + t19 * t12 - t20 * t11) * g(1);
U = t3;
