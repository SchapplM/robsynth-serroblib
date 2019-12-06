% Calculate potential energy for
% S5PRPRR7
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRPRR7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR7_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:38
% EndTime: 2019-12-05 15:59:39
% DurationCPUTime: 0.49s
% Computational Cost: add. (134->65), mult. (184->53), div. (0->0), fcn. (156->8), ass. (0->27)
t50 = -mrSges(3,1) + mrSges(4,2) - mrSges(5,3) + m(6) * (-pkin(7) - pkin(6)) - mrSges(6,3);
t16 = sin(qJ(4));
t18 = cos(qJ(4));
t13 = qJ(4) + qJ(5);
t6 = sin(t13);
t7 = cos(t13);
t49 = -t6 * mrSges(6,1) - t18 * mrSges(5,2) - t7 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3) + (-m(6) * pkin(4) - mrSges(5,1)) * t16;
t46 = m(5) * pkin(6);
t45 = -m(1) - m(2);
t14 = sin(pkin(8));
t17 = sin(qJ(2));
t32 = qJ(3) * t17;
t19 = cos(qJ(2));
t37 = t14 * t19;
t44 = pkin(2) * t37 + t14 * t32;
t43 = m(4) + m(5) + m(6);
t40 = t18 * mrSges(5,1) + t7 * mrSges(6,1) - t16 * mrSges(5,2) - t6 * mrSges(6,2) + mrSges(4,1) - mrSges(2,2) + mrSges(3,3);
t39 = t49 * t17 + t50 * t19 - mrSges(2,1);
t15 = cos(pkin(8));
t36 = t15 * t19;
t31 = t14 * pkin(1) + r_base(2);
t12 = qJ(1) + r_base(3);
t30 = t15 * pkin(1) + t14 * pkin(5) + r_base(1);
t28 = t31 + t44;
t25 = -t15 * pkin(5) + t31;
t5 = t18 * pkin(4) + pkin(3);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t12 - t43 * (t17 * pkin(2) + t12) + (t43 * qJ(3) - t49) * t19 + (-t46 + t50) * t17) * g(3) + (-mrSges(1,2) - m(3) * t25 - m(4) * (t25 + t44) - m(5) * (pkin(6) * t37 + t28) - m(6) * t28 + t45 * r_base(2) + (-m(5) * (-pkin(3) - pkin(5)) - m(6) * (-pkin(5) - t5) + t40) * t15 + t39 * t14) * g(2) + (-t36 * t46 - m(3) * t30 - mrSges(1,1) + t45 * r_base(1) - t43 * (pkin(2) * t36 + t15 * t32 + t30) + (-m(5) * pkin(3) - m(6) * t5 - t40) * t14 + t39 * t15) * g(1);
U = t1;
