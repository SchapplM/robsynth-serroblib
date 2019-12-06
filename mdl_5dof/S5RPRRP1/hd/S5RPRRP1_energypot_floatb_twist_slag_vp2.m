% Calculate potential energy for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPRRP1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_energypot_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP1_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:10
% EndTime: 2019-12-05 17:59:10
% DurationCPUTime: 0.32s
% Computational Cost: add. (114->50), mult. (115->33), div. (0->0), fcn. (79->6), ass. (0->18)
t31 = m(5) + m(6);
t30 = mrSges(5,2) + mrSges(6,2);
t29 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t28 = -m(1) - m(2);
t26 = m(3) + m(4) + t31;
t11 = sin(qJ(3));
t13 = cos(qJ(3));
t10 = qJ(3) + qJ(4);
t2 = sin(t10);
t3 = cos(t10);
t25 = t13 * mrSges(4,2) + t29 * t2 + t30 * t3 - mrSges(2,2) + mrSges(3,3) + (pkin(3) * t31 + mrSges(4,1)) * t11;
t15 = -pkin(7) - pkin(6);
t24 = -m(4) * pkin(6) + m(5) * t15 + m(6) * (-qJ(5) + t15) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t9 = pkin(5) + r_base(3);
t21 = pkin(2) + t9;
t14 = cos(qJ(1));
t12 = sin(qJ(1));
t1 = (-m(1) * r_base(3) - m(4) * t21 - t13 * mrSges(4,1) + t11 * mrSges(4,2) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t9 - t31 * (t13 * pkin(3) + t21) - t29 * t3 + t30 * t2) * g(3) + (-mrSges(1,2) + t28 * r_base(2) - t26 * (t12 * pkin(1) + r_base(2)) + (t26 * qJ(2) + t25) * t14 + t24 * t12) * g(2) + (-mrSges(1,1) + t28 * r_base(1) - t26 * (t14 * pkin(1) + t12 * qJ(2) + r_base(1)) + t24 * t14 - t25 * t12) * g(1);
U = t1;
