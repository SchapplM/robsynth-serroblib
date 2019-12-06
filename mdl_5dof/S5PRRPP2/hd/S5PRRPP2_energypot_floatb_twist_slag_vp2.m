% Calculate potential energy for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPP2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRPP2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP2_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP2_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:08:37
% EndTime: 2019-12-05 16:08:38
% DurationCPUTime: 0.47s
% Computational Cost: add. (167->65), mult. (207->59), div. (0->0), fcn. (187->8), ass. (0->28)
t22 = sin(qJ(3));
t24 = cos(qJ(3));
t50 = -m(4) * pkin(2) - t24 * mrSges(4,1) + t22 * mrSges(4,2) - mrSges(3,1);
t49 = -m(4) * pkin(6) + mrSges(3,2) - mrSges(6,2) - mrSges(4,3) - mrSges(5,3);
t48 = -m(1) - m(2);
t47 = m(3) + m(4);
t46 = -m(5) - m(6);
t45 = -t22 * mrSges(4,1) - t24 * mrSges(4,2) + mrSges(2,2) - mrSges(3,3);
t43 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t42 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t23 = sin(qJ(2));
t25 = cos(qJ(2));
t41 = t49 * t23 + t50 * t25 - mrSges(2,1);
t40 = pkin(3) * t22;
t19 = sin(pkin(7));
t39 = t19 * t25;
t20 = cos(pkin(7));
t38 = t20 * t25;
t21 = -qJ(4) - pkin(6);
t37 = t21 * t23;
t34 = t19 * pkin(1) + r_base(2);
t17 = qJ(1) + r_base(3);
t33 = t20 * pkin(1) + t19 * pkin(5) + r_base(1);
t18 = qJ(3) + pkin(8);
t13 = cos(t18);
t12 = sin(t18);
t11 = t24 * pkin(3) + pkin(2);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + t46 * (t23 * t11 + t25 * t21 + t17) + (-m(2) - t47) * t17 - t49 * t25 + (t42 * t12 + t43 * t13 + t50) * t23) * g(3) + (-mrSges(1,2) + t48 * r_base(2) - t47 * t34 + t46 * (-t19 * t37 + t11 * t39 + (-pkin(5) - t40) * t20 + t34) + t43 * (-t20 * t12 + t13 * t39) + t42 * (t12 * t39 + t20 * t13) + (t47 * pkin(5) - t45) * t20 + t41 * t19) * g(2) + (-mrSges(1,1) + t48 * r_base(1) + t46 * (t11 * t38 + t19 * t40 - t20 * t37 + t33) + t43 * (t19 * t12 + t13 * t38) - t47 * t33 + t42 * (t12 * t38 - t19 * t13) + t45 * t19 + t41 * t20) * g(1);
U = t1;
