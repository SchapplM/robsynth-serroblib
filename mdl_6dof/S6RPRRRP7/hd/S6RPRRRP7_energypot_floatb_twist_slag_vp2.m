% Calculate potential energy for
% S6RPRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRP7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP7_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP7_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:18:32
% EndTime: 2019-03-09 06:18:33
% DurationCPUTime: 0.59s
% Computational Cost: add. (260->75), mult. (241->63), div. (0->0), fcn. (215->10), ass. (0->34)
t68 = m(5) * pkin(8) - mrSges(4,2) + mrSges(7,2) + mrSges(5,3) + mrSges(6,3);
t27 = sin(qJ(4));
t29 = cos(qJ(4));
t67 = -m(5) * pkin(3) - t29 * mrSges(5,1) + t27 * mrSges(5,2) - mrSges(4,1);
t61 = -m(6) - m(7);
t65 = m(3) * qJ(2) + mrSges(5,2) * t29 - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + (-pkin(4) * t61 + mrSges(5,1)) * t27;
t15 = pkin(4) * t29 + pkin(3);
t21 = pkin(10) + qJ(3);
t16 = sin(t21);
t17 = cos(t21);
t24 = sin(pkin(10));
t25 = cos(pkin(10));
t31 = -pkin(9) - pkin(8);
t64 = -m(3) * pkin(1) - t25 * mrSges(3,1) + t24 * mrSges(3,2) - mrSges(2,1) + (t61 * t15 + t67) * t17 + (-t61 * t31 - t68) * t16;
t63 = -m(2) - m(3);
t62 = -m(4) - m(5);
t59 = -m(1) + t63;
t55 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t54 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t23 = qJ(4) + qJ(5);
t18 = sin(t23);
t28 = sin(qJ(1));
t47 = t18 * t28;
t30 = cos(qJ(1));
t46 = t18 * t30;
t19 = cos(t23);
t45 = t19 * t30;
t44 = t28 * t19;
t22 = pkin(6) + r_base(3);
t13 = pkin(2) * t25 + pkin(1);
t43 = t30 * t13 + r_base(1);
t42 = t24 * pkin(2) + t22;
t26 = -pkin(7) - qJ(2);
t1 = (-m(1) * r_base(3) - t24 * mrSges(3,1) - t25 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t62 * t42 + t61 * (t16 * t15 + t17 * t31 + t42) + t63 * t22 + t68 * t17 + (t18 * t54 + t55 * t19 + t67) * t16) * g(3) + (-mrSges(1,2) + t55 * (t17 * t44 - t46) + t54 * (t17 * t47 + t45) + t59 * r_base(2) + (t62 + t61) * (t28 * t13 + t30 * t26 + r_base(2))) * g(2) + (-m(5) * t43 - mrSges(1,1) + t55 * (t17 * t45 + t47) + t54 * (t17 * t46 - t44) + t59 * r_base(1) + (-m(4) + t61) * (-t28 * t26 + t43)) * g(1) + (t64 * g(1) + t65 * g(2)) * t30 + (t64 * g(2) + (m(5) * t26 - t65) * g(1)) * t28;
U  = t1;
