% Calculate potential energy for
% S5RRPRR15
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR15_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRPRR15_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR15_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR15_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:04
% EndTime: 2019-12-31 20:41:04
% DurationCPUTime: 0.50s
% Computational Cost: add. (134->72), mult. (184->65), div. (0->0), fcn. (156->8), ass. (0->31)
t49 = -mrSges(3,1) + mrSges(4,2) + m(6) * (-pkin(8) - pkin(7)) - mrSges(6,3);
t14 = sin(qJ(4));
t13 = qJ(4) + qJ(5);
t6 = sin(t13);
t7 = cos(t13);
t48 = -m(6) * pkin(4) * t14 - t6 * mrSges(6,1) - t7 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3);
t47 = -m(1) - m(2);
t16 = sin(qJ(1));
t15 = sin(qJ(2));
t32 = qJ(3) * t15;
t18 = cos(qJ(2));
t36 = t16 * t18;
t46 = pkin(2) * t36 + t16 * t32;
t45 = m(4) + m(5) + m(6);
t44 = -m(5) * pkin(7) - mrSges(5,3);
t41 = t48 * t15 + t49 * t18 - mrSges(2,1);
t40 = t7 * mrSges(6,1) - t6 * mrSges(6,2) + mrSges(4,1) - mrSges(2,2) + mrSges(3,3);
t38 = t16 * t14;
t17 = cos(qJ(4));
t37 = t16 * t17;
t19 = cos(qJ(1));
t35 = t19 * t14;
t34 = t19 * t17;
t33 = t19 * t18;
t12 = pkin(5) + r_base(3);
t31 = t16 * pkin(1) + r_base(2);
t29 = t19 * pkin(1) + t16 * pkin(6) + r_base(1);
t28 = t31 + t46;
t22 = -t19 * pkin(6) + t31;
t5 = t17 * pkin(4) + pkin(3);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t12 - t45 * (t15 * pkin(2) + t12) + (t14 * mrSges(5,1) + t17 * mrSges(5,2) + qJ(3) * t45 - t48) * t18 + (t44 + t49) * t15) * g(3) + (-mrSges(1,2) - m(3) * t22 - m(4) * (t22 + t46) - m(5) * (pkin(7) * t36 + t28) - (t15 * t38 - t34) * mrSges(5,1) - (t15 * t37 + t35) * mrSges(5,2) - mrSges(5,3) * t36 - m(6) * t28 + t47 * r_base(2) + (-m(5) * (-pkin(3) - pkin(6)) - m(6) * (-pkin(6) - t5) + t40) * t19 + t41 * t16) * g(2) + (-mrSges(1,1) - m(3) * t29 - (t15 * t35 + t37) * mrSges(5,1) - (t15 * t34 - t38) * mrSges(5,2) + t47 * r_base(1) + t44 * t33 - t45 * (pkin(2) * t33 + t19 * t32 + t29) + t41 * t19 + (-m(5) * pkin(3) - m(6) * t5 - t40) * t16) * g(1);
U = t1;
