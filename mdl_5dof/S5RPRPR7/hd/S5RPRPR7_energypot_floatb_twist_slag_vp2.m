% Calculate potential energy for
% S5RPRPR7
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRPR7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR7_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:18:50
% EndTime: 2019-12-31 18:18:50
% DurationCPUTime: 0.33s
% Computational Cost: add. (167->57), mult. (129->45), div. (0->0), fcn. (97->10), ass. (0->26)
t16 = sin(qJ(5));
t19 = cos(qJ(5));
t40 = -m(6) * pkin(4) - t19 * mrSges(6,1) + t16 * mrSges(6,2) - mrSges(5,1);
t39 = -m(6) * pkin(7) + mrSges(5,2) - mrSges(6,3);
t38 = -m(1) - m(2);
t37 = -m(3) - m(4);
t36 = m(5) + m(6);
t17 = sin(qJ(3));
t20 = cos(qJ(3));
t13 = qJ(3) + pkin(9);
t5 = sin(t13);
t7 = cos(t13);
t34 = -m(4) * pkin(2) - t20 * mrSges(4,1) + t17 * mrSges(4,2) + t39 * t5 + t40 * t7 - mrSges(3,1);
t33 = m(4) * pkin(6) + t16 * mrSges(6,1) + t19 * mrSges(6,2) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t32 = pkin(5) + r_base(3);
t18 = sin(qJ(1));
t31 = t18 * pkin(1) + r_base(2);
t21 = cos(qJ(1));
t30 = t21 * pkin(1) + r_base(1);
t9 = qJ(2) + t32;
t15 = -qJ(4) - pkin(6);
t14 = qJ(1) + pkin(8);
t8 = cos(t14);
t6 = sin(t14);
t4 = pkin(3) * t20 + pkin(2);
t1 = (-m(1) * r_base(3) - m(2) * t32 - mrSges(4,1) * t17 - mrSges(4,2) * t20 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t37 * t9 - t36 * (t17 * pkin(3) + t9) - t39 * t7 + t40 * t5) * g(3) + (-t18 * mrSges(2,1) - mrSges(2,2) * t21 - mrSges(1,2) + t38 * r_base(2) + t37 * t31 - t36 * (t8 * t15 + t6 * t4 + t31) + t33 * t8 + t34 * t6) * g(2) + (-mrSges(2,1) * t21 + t18 * mrSges(2,2) - mrSges(1,1) + t38 * r_base(1) + t37 * t30 - t36 * (t8 * t4 + t30) + t34 * t8 + (t36 * t15 - t33) * t6) * g(1);
U = t1;
