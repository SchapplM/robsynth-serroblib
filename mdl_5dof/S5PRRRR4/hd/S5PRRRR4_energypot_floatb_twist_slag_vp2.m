% Calculate potential energy for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRR4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR4_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR4_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:37
% EndTime: 2019-12-05 17:07:37
% DurationCPUTime: 0.26s
% Computational Cost: add. (154->56), mult. (97->43), div. (0->0), fcn. (61->10), ass. (0->22)
t32 = -m(1) - m(2);
t31 = -m(4) - m(5) - m(6);
t15 = qJ(4) + qJ(5);
t10 = sin(t15);
t11 = cos(t15);
t18 = sin(qJ(4));
t19 = cos(qJ(4));
t30 = -mrSges(4,1) - m(5) * pkin(3) - mrSges(5,1) * t19 + mrSges(5,2) * t18 - m(6) * (pkin(4) * t19 + pkin(3)) - t11 * mrSges(6,1) + t10 * mrSges(6,2);
t29 = m(5) * pkin(7) - m(6) * (-pkin(8) - pkin(7)) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t14 = pkin(9) + qJ(2);
t16 = sin(pkin(9));
t28 = t16 * pkin(1) + r_base(2);
t17 = cos(pkin(9));
t27 = t17 * pkin(1) + r_base(1);
t26 = qJ(1) + r_base(3);
t23 = pkin(5) + t26;
t9 = qJ(3) + t14;
t8 = cos(t14);
t7 = sin(t14);
t4 = cos(t9);
t3 = sin(t9);
t1 = (-m(1) * r_base(3) - m(2) * t26 - m(3) * t23 - t10 * mrSges(6,1) - mrSges(5,2) * t19 - t11 * mrSges(6,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) + (-m(6) * pkin(4) - mrSges(5,1)) * t18 + t31 * (pkin(6) + t23)) * g(3) + (-m(3) * t28 - mrSges(2,1) * t16 - t7 * mrSges(3,1) - mrSges(2,2) * t17 - t8 * mrSges(3,2) - mrSges(1,2) + t32 * r_base(2) + t31 * (pkin(2) * t7 + t28) + t29 * t4 + t30 * t3) * g(2) + (-m(3) * t27 - mrSges(2,1) * t17 - t8 * mrSges(3,1) + mrSges(2,2) * t16 + t7 * mrSges(3,2) - mrSges(1,1) + t32 * r_base(1) + t30 * t4 + t31 * (pkin(2) * t8 + t27) - t29 * t3) * g(1);
U = t1;
