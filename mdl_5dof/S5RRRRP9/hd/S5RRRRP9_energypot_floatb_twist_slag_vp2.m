% Calculate potential energy for
% S5RRRRP9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP9_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRP9_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP9_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP9_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP9_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:03:42
% EndTime: 2019-12-31 22:03:42
% DurationCPUTime: 0.48s
% Computational Cost: add. (167->67), mult. (207->62), div. (0->0), fcn. (187->8), ass. (0->31)
t19 = sin(qJ(3));
t22 = cos(qJ(3));
t51 = -m(4) * pkin(2) - t22 * mrSges(4,1) + t19 * mrSges(4,2) - mrSges(3,1);
t50 = -m(4) * pkin(7) + mrSges(3,2) - mrSges(4,3);
t49 = -m(1) - m(2);
t48 = m(3) + m(4);
t47 = -m(5) - m(6);
t46 = -mrSges(5,3) - mrSges(6,2);
t20 = sin(qJ(2));
t23 = cos(qJ(2));
t45 = t50 * t20 + t51 * t23 - mrSges(2,1);
t44 = -t19 * mrSges(4,1) - t22 * mrSges(4,2) + mrSges(2,2) - mrSges(3,3);
t42 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t41 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t40 = pkin(3) * t19;
t21 = sin(qJ(1));
t39 = t21 * t20;
t38 = t21 * t23;
t18 = qJ(3) + qJ(4);
t12 = sin(t18);
t24 = cos(qJ(1));
t37 = t24 * t12;
t13 = cos(t18);
t36 = t24 * t13;
t35 = t24 * t20;
t17 = pkin(5) + r_base(3);
t34 = t21 * pkin(1) + r_base(2);
t33 = t24 * pkin(1) + t21 * pkin(6) + r_base(1);
t25 = -pkin(8) - pkin(7);
t10 = t22 * pkin(3) + pkin(2);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + t47 * (t20 * t10 + t23 * t25 + t17) + (-m(2) - t48) * t17 + (-t46 - t50) * t23 + (t41 * t12 + t42 * t13 + t51) * t20) * g(3) + (-mrSges(1,2) + t49 * r_base(2) + t46 * t39 - t48 * t34 + t47 * (-t25 * t39 + t10 * t38 + (-pkin(6) - t40) * t24 + t34) + t42 * (t13 * t38 - t37) + t41 * (t12 * t38 + t36) + (t48 * pkin(6) - t44) * t24 + t45 * t21) * g(2) + (-mrSges(1,1) + t49 * r_base(1) + t47 * (t24 * t23 * t10 + t21 * t40 - t25 * t35 + t33) + t42 * (t21 * t12 + t23 * t36) + t46 * t35 - t48 * t33 + t41 * (-t21 * t13 + t23 * t37) + t45 * t24 + t44 * t21) * g(1);
U = t1;
