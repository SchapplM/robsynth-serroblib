% Calculate potential energy for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
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
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_energypot_floatb_twist_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR3_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:55:01
% EndTime: 2019-12-05 18:55:02
% DurationCPUTime: 0.45s
% Computational Cost: add. (144->61), mult. (156->58), div. (0->0), fcn. (128->10), ass. (0->31)
t17 = cos(qJ(4));
t12 = qJ(4) + qJ(5);
t6 = sin(t12);
t8 = cos(t12);
t46 = -mrSges(4,1) - m(5) * pkin(2) - m(6) * (pkin(3) * t17 + pkin(2)) - t8 * mrSges(6,1) + t6 * mrSges(6,2);
t45 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t15 = sin(qJ(2));
t18 = cos(qJ(2));
t13 = qJ(2) + qJ(3);
t7 = sin(t13);
t9 = cos(t13);
t44 = -mrSges(3,1) * t18 + mrSges(3,2) * t15 + t45 * t7 + t46 * t9 - mrSges(2,1);
t43 = -m(2) - m(3);
t42 = -m(1) + t43;
t40 = t6 * mrSges(6,1) + t8 * mrSges(6,2) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3);
t37 = pkin(5) * t7;
t36 = pkin(1) * t18;
t14 = sin(qJ(4));
t19 = cos(qJ(1));
t31 = t14 * t19;
t16 = sin(qJ(1));
t30 = t16 * t14;
t29 = t16 * t17;
t28 = t17 * t19;
t11 = pkin(4) + r_base(3);
t27 = t16 * t36 + r_base(2);
t26 = t19 * t36 + r_base(1);
t25 = t16 * t37 + t27;
t24 = t19 * t37 + t26;
t23 = t15 * pkin(1) + t11;
t1 = (-m(1) * r_base(3) - m(4) * t23 - mrSges(3,1) * t15 - mrSges(3,2) * t18 - mrSges(1,3) - mrSges(2,3) + (-m(5) - m(6)) * (-t9 * pkin(5) + t23) + t43 * t11 - t45 * t9 + (-t17 * mrSges(5,1) + t14 * mrSges(5,2) + t46) * t7) * g(3) + (-mrSges(1,2) - m(4) * t27 - m(5) * t25 - (t29 * t9 - t31) * mrSges(5,1) - (-t30 * t9 - t28) * mrSges(5,2) - m(6) * (-pkin(3) * t31 + t25) + t42 * r_base(2) + t40 * t19 + t44 * t16) * g(2) + (-mrSges(1,1) - m(4) * t26 - m(5) * t24 - (t28 * t9 + t30) * mrSges(5,1) - (-t31 * t9 + t29) * mrSges(5,2) - m(6) * (pkin(3) * t30 + t24) + t42 * r_base(1) - t40 * t16 + t44 * t19) * g(1);
U = t1;
