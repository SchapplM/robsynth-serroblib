% Calculate potential energy for
% S6RPRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPP4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:38:45
% EndTime: 2019-03-09 04:38:45
% DurationCPUTime: 0.49s
% Computational Cost: add. (260->80), mult. (241->73), div. (0->0), fcn. (215->10), ass. (0->40)
t28 = sin(qJ(4));
t30 = cos(qJ(4));
t62 = -m(5) * pkin(3) - t30 * mrSges(5,1) + t28 * mrSges(5,2) - mrSges(4,1);
t61 = -m(5) * pkin(8) + mrSges(4,2) - mrSges(5,3);
t60 = -m(2) - m(3);
t59 = -m(4) - m(5);
t58 = -m(6) - m(7);
t57 = -mrSges(6,3) - mrSges(7,2);
t56 = -m(1) + t60;
t21 = pkin(9) + qJ(3);
t16 = sin(t21);
t18 = cos(t21);
t24 = sin(pkin(9));
t25 = cos(pkin(9));
t54 = -m(3) * pkin(1) - t25 * mrSges(3,1) + t24 * mrSges(3,2) + t61 * t16 + t62 * t18 - mrSges(2,1);
t53 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t52 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t51 = m(3) * qJ(2) + mrSges(5,1) * t28 + mrSges(5,2) * t30 - mrSges(2,2) + mrSges(3,3) + mrSges(4,3);
t50 = pkin(4) * t28;
t29 = sin(qJ(1));
t49 = t16 * t29;
t31 = cos(qJ(1));
t48 = t16 * t31;
t47 = t18 * t31;
t22 = qJ(4) + pkin(10);
t19 = cos(t22);
t46 = t19 * t31;
t17 = sin(t22);
t45 = t29 * t17;
t44 = t29 * t19;
t23 = pkin(6) + r_base(3);
t13 = pkin(2) * t25 + pkin(1);
t43 = t31 * t13 + r_base(1);
t42 = t24 * pkin(2) + t23;
t27 = -pkin(7) - qJ(2);
t41 = t29 * t13 + t31 * t27 + r_base(2);
t38 = -t29 * t27 + t43;
t26 = -qJ(5) - pkin(8);
t15 = pkin(4) * t30 + pkin(3);
t1 = (-m(1) * r_base(3) - t24 * mrSges(3,1) - t25 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t59 * t42 + t58 * (t16 * t15 + t18 * t26 + t42) + t60 * t23 + (-t57 - t61) * t18 + (t52 * t17 + t53 * t19 + t62) * t16) * g(3) + (-mrSges(1,2) + t57 * t49 + t59 * t41 + t58 * (t29 * t18 * t15 - t26 * t49 - t31 * t50 + t41) + t53 * (-t17 * t31 + t18 * t44) + t52 * (t18 * t45 + t46) + t56 * r_base(2) + t51 * t31 + t54 * t29) * g(2) + (-m(4) * t38 - m(5) * t43 - mrSges(1,1) + t57 * t48 + t58 * (t15 * t47 - t26 * t48 + t29 * t50 + t38) + t53 * (t18 * t46 + t45) + t52 * (t17 * t47 - t44) + t56 * r_base(1) + t54 * t31 + (m(5) * t27 - t51) * t29) * g(1);
U  = t1;
