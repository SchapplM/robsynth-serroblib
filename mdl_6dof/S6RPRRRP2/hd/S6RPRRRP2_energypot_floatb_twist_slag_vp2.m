% Calculate potential energy for
% S6RPRRRP2
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
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRP2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:58:41
% EndTime: 2019-03-09 05:58:42
% DurationCPUTime: 0.49s
% Computational Cost: add. (256->74), mult. (212->64), div. (0->0), fcn. (182->10), ass. (0->31)
t25 = cos(qJ(4));
t10 = t25 * pkin(4) + pkin(3);
t21 = qJ(4) + qJ(5);
t15 = cos(t21);
t22 = sin(qJ(4));
t53 = -m(6) * t10 - m(7) * (pkin(5) * t15 + t10) - mrSges(4,1) - m(5) * pkin(3) - t25 * mrSges(5,1) + t22 * mrSges(5,2);
t28 = -pkin(9) - pkin(8);
t52 = m(6) * t28 + m(7) * (-qJ(6) + t28) + mrSges(4,2) - mrSges(6,3) - mrSges(7,3) - m(5) * pkin(8) - mrSges(5,3);
t51 = -m(1) - m(2);
t50 = -mrSges(6,1) - mrSges(7,1);
t49 = mrSges(6,2) + mrSges(7,2);
t48 = -m(4) - m(6) - m(7);
t47 = -m(5) + t48;
t23 = sin(qJ(3));
t26 = cos(qJ(3));
t45 = t52 * t23 + t53 * t26 - mrSges(3,1);
t14 = sin(t21);
t43 = pkin(4) * t22;
t44 = -m(6) * t43 - m(7) * (pkin(5) * t14 + t43) + mrSges(3,2) - mrSges(4,3) - t22 * mrSges(5,1) - t25 * mrSges(5,2);
t42 = t14 * t26;
t41 = t15 * t26;
t40 = pkin(6) + r_base(3);
t24 = sin(qJ(1));
t39 = t24 * pkin(1) + r_base(2);
t27 = cos(qJ(1));
t38 = t27 * pkin(1) + r_base(1);
t20 = qJ(1) + pkin(10);
t11 = sin(t20);
t37 = t11 * pkin(2) + t39;
t12 = cos(t20);
t1 = (-m(1) * r_base(3) - m(2) * t40 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + (-m(3) + t47) * (qJ(2) + t40) - t52 * t26 + (t49 * t14 + t50 * t15 + t53) * t23) * g(3) + (-m(3) * t39 - m(5) * t37 - t24 * mrSges(2,1) - t27 * mrSges(2,2) - mrSges(1,2) + t51 * r_base(2) + t50 * (t11 * t41 - t12 * t14) - t49 * (-t11 * t42 - t12 * t15) + t48 * (-t12 * pkin(7) + t37) + (m(5) * pkin(7) - t44) * t12 + t45 * t11) * g(2) + (-m(3) * t38 - t27 * mrSges(2,1) + t24 * mrSges(2,2) - mrSges(1,1) + t51 * r_base(1) + t50 * (t11 * t14 + t12 * t41) - t49 * (t11 * t15 - t12 * t42) + t47 * (t12 * pkin(2) + t11 * pkin(7) + t38) + t44 * t11 + t45 * t12) * g(1);
U  = t1;
