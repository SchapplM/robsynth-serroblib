% Calculate potential energy for
% S5PRRRP3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRP3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRP3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP3_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:26
% EndTime: 2019-12-05 16:43:27
% DurationCPUTime: 0.29s
% Computational Cost: add. (147->53), mult. (109->37), div. (0->0), fcn. (73->8), ass. (0->22)
t30 = -m(5) - m(6);
t34 = mrSges(5,2) + mrSges(6,2);
t33 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t32 = -m(1) - m(2);
t31 = -m(3) - m(4);
t29 = t30 + t31;
t17 = sin(qJ(3));
t18 = cos(qJ(3));
t14 = qJ(3) + qJ(4);
t6 = sin(t14);
t7 = cos(t14);
t28 = -m(4) * pkin(2) - t18 * mrSges(4,1) + t17 * mrSges(4,2) - mrSges(3,1) + t30 * (t18 * pkin(3) + pkin(2)) + t33 * t7 + t34 * t6;
t19 = -pkin(7) - pkin(6);
t27 = m(4) * pkin(6) - m(5) * t19 - m(6) * (-qJ(5) + t19) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t24 = qJ(1) + r_base(3);
t5 = pkin(5) + t24;
t16 = cos(pkin(8));
t15 = sin(pkin(8));
t13 = pkin(8) + qJ(2);
t4 = cos(t13);
t3 = sin(t13);
t1 = (-m(1) * r_base(3) - m(2) * t24 - mrSges(4,1) * t17 - mrSges(4,2) * t18 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - t34 * t7 + t30 * (t17 * pkin(3) + t5) + t33 * t6 + t31 * t5) * g(3) + (-mrSges(2,1) * t15 - mrSges(2,2) * t16 - mrSges(1,2) + t32 * r_base(2) + t29 * (t15 * pkin(1) + r_base(2)) + t27 * t4 + t28 * t3) * g(2) + (-mrSges(2,1) * t16 + mrSges(2,2) * t15 - mrSges(1,1) + t32 * r_base(1) + t28 * t4 + t29 * (t16 * pkin(1) + r_base(1)) - t27 * t3) * g(1);
U = t1;
