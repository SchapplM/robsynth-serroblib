% Calculate potential energy for
% S6RPRRRP3
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
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRP3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:02:48
% EndTime: 2019-03-09 06:02:49
% DurationCPUTime: 0.52s
% Computational Cost: add. (271->76), mult. (227->64), div. (0->0), fcn. (201->10), ass. (0->33)
t24 = sin(qJ(4));
t27 = cos(qJ(4));
t60 = -m(5) * pkin(3) - t27 * mrSges(5,1) + t24 * mrSges(5,2) - mrSges(4,1);
t59 = -m(5) * pkin(8) + mrSges(4,2) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t13 = pkin(4) * t27 + pkin(3);
t25 = sin(qJ(3));
t28 = cos(qJ(3));
t30 = -pkin(9) - pkin(8);
t58 = t13 * t28 - t25 * t30;
t57 = -m(1) - m(2);
t56 = m(4) + m(5);
t55 = -m(6) - m(7);
t54 = -t24 * mrSges(5,1) - t27 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t52 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t51 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t50 = t59 * t25 + t60 * t28 - mrSges(3,1);
t49 = pkin(4) * t24;
t23 = qJ(4) + qJ(5);
t18 = sin(t23);
t47 = t18 * t28;
t19 = cos(t23);
t46 = t19 * t28;
t42 = pkin(6) + r_base(3);
t26 = sin(qJ(1));
t41 = t26 * pkin(1) + r_base(2);
t29 = cos(qJ(1));
t40 = t29 * pkin(1) + r_base(1);
t17 = qJ(2) + t42;
t22 = qJ(1) + pkin(10);
t15 = sin(t22);
t39 = t15 * pkin(2) + t41;
t16 = cos(t22);
t1 = (-m(1) * r_base(3) - m(2) * t42 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t55 * (t25 * t13 + t28 * t30 + t17) + (-m(3) - t56) * t17 - t59 * t28 + (t18 * t51 + t19 * t52 + t60) * t25) * g(3) + (-m(3) * t41 - t26 * mrSges(2,1) - t29 * mrSges(2,2) - mrSges(1,2) + t57 * r_base(2) - t56 * t39 + t55 * ((-pkin(7) - t49) * t16 + t39 + t58 * t15) + t52 * (t15 * t46 - t16 * t18) + t51 * (t15 * t47 + t16 * t19) + (t56 * pkin(7) - t54) * t16 + t50 * t15) * g(2) + (-m(3) * t40 - t29 * mrSges(2,1) + t26 * mrSges(2,2) + t57 * r_base(1) - mrSges(1,1) + (t55 - t56) * (t16 * pkin(2) + t15 * pkin(7) + t40) + (t52 * t18 - t51 * t19 + t55 * t49 + t54) * t15 + (t52 * t46 + t51 * t47 + t55 * t58 + t50) * t16) * g(1);
U  = t1;
