% Calculate potential energy for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_energypot_floatb_twist_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR2_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR2_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:52:52
% EndTime: 2019-12-05 18:52:52
% DurationCPUTime: 0.17s
% Computational Cost: add. (118->46), mult. (111->37), div. (0->0), fcn. (79->10), ass. (0->23)
t29 = -m(1) - m(2);
t28 = -m(3) - m(4);
t27 = -m(5) - m(6);
t18 = cos(qJ(3));
t26 = pkin(2) * t18;
t25 = mrSges(5,2) - mrSges(6,3);
t16 = sin(qJ(1));
t24 = t16 * pkin(1) + r_base(2);
t19 = cos(qJ(1));
t23 = t19 * pkin(1) + r_base(1);
t14 = sin(qJ(5));
t17 = cos(qJ(5));
t22 = -mrSges(6,1) * t17 + mrSges(6,2) * t14 - mrSges(5,1);
t21 = t14 * mrSges(6,1) + t17 * mrSges(6,2) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t15 = sin(qJ(3));
t12 = qJ(3) + qJ(4);
t6 = sin(t12);
t8 = cos(t12);
t20 = -mrSges(4,1) * t18 + mrSges(4,2) * t15 + t22 * t8 + t25 * t6 - mrSges(3,1);
t13 = qJ(1) + qJ(2);
t9 = cos(t13);
t7 = sin(t13);
t1 = (-mrSges(4,1) * t15 - mrSges(4,2) * t18 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - t25 * t8 + t27 * (pkin(2) * t15 + r_base(3)) + t22 * t6 + (t28 + t29) * r_base(3)) * g(3) + (-t16 * mrSges(2,1) - mrSges(2,2) * t19 - mrSges(1,2) + t29 * r_base(2) + t28 * t24 + t27 * (t7 * t26 + t24) + t21 * t9 + t20 * t7) * g(2) + (-mrSges(2,1) * t19 + t16 * mrSges(2,2) - mrSges(1,1) + t29 * r_base(1) + t28 * t23 + t27 * (t9 * t26 + t23) - t21 * t7 + t20 * t9) * g(1);
U = t1;
