% Calculate potential energy for
% S5RRRPP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPP6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRPP6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP6_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP6_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:59:54
% EndTime: 2019-12-31 20:59:54
% DurationCPUTime: 0.49s
% Computational Cost: add. (167->67), mult. (207->62), div. (0->0), fcn. (187->8), ass. (0->31)
t20 = sin(qJ(3));
t23 = cos(qJ(3));
t51 = -m(4) * pkin(2) - t23 * mrSges(4,1) + t20 * mrSges(4,2) - mrSges(3,1);
t50 = -m(4) * pkin(7) + mrSges(3,2) - mrSges(4,3);
t49 = -m(1) - m(2);
t48 = m(3) + m(4);
t47 = -m(5) - m(6);
t46 = -mrSges(5,3) - mrSges(6,2);
t21 = sin(qJ(2));
t24 = cos(qJ(2));
t45 = t50 * t21 + t51 * t24 - mrSges(2,1);
t44 = -t20 * mrSges(4,1) - t23 * mrSges(4,2) + mrSges(2,2) - mrSges(3,3);
t42 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t41 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t40 = pkin(3) * t20;
t22 = sin(qJ(1));
t39 = t22 * t21;
t38 = t22 * t24;
t17 = qJ(3) + pkin(8);
t12 = sin(t17);
t25 = cos(qJ(1));
t37 = t25 * t12;
t13 = cos(t17);
t36 = t25 * t13;
t35 = t25 * t21;
t18 = pkin(5) + r_base(3);
t34 = t22 * pkin(1) + r_base(2);
t33 = t25 * pkin(1) + t22 * pkin(6) + r_base(1);
t19 = -qJ(4) - pkin(7);
t11 = t23 * pkin(3) + pkin(2);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + t47 * (t21 * t11 + t24 * t19 + t18) + (-m(2) - t48) * t18 + (-t46 - t50) * t24 + (t41 * t12 + t42 * t13 + t51) * t21) * g(3) + (-mrSges(1,2) + t49 * r_base(2) + t46 * t39 - t48 * t34 + t47 * (-t19 * t39 + t11 * t38 + (-pkin(6) - t40) * t25 + t34) + t42 * (t13 * t38 - t37) + t41 * (t12 * t38 + t36) + (t48 * pkin(6) - t44) * t25 + t45 * t22) * g(2) + (-mrSges(1,1) + t49 * r_base(1) + t47 * (t25 * t24 * t11 - t19 * t35 + t22 * t40 + t33) + t42 * (t22 * t12 + t24 * t36) + t46 * t35 - t48 * t33 + t41 * (-t22 * t13 + t24 * t37) + t45 * t25 + t44 * t22) * g(1);
U = t1;
