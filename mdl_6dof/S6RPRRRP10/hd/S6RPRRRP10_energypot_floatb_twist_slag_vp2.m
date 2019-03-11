% Calculate potential energy for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP10_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRP10_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP10_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP10_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP10_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:30:12
% EndTime: 2019-03-09 06:30:12
% DurationCPUTime: 0.57s
% Computational Cost: add. (192->75), mult. (229->63), div. (0->0), fcn. (203->8), ass. (0->32)
t21 = sin(qJ(4));
t24 = cos(qJ(4));
t56 = -m(5) * pkin(3) - t24 * mrSges(5,1) + t21 * mrSges(5,2) - mrSges(4,1);
t55 = -m(5) * pkin(8) + mrSges(4,2) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t54 = -m(1) - m(2);
t53 = -m(4) - m(5);
t52 = -m(6) - m(7);
t50 = -t21 * mrSges(5,1) - t24 * mrSges(5,2) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t49 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t48 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t22 = sin(qJ(3));
t25 = cos(qJ(3));
t47 = t56 * t22 - t55 * t25 + mrSges(2,2) - mrSges(3,3);
t46 = pkin(4) * t21;
t20 = qJ(4) + qJ(5);
t12 = sin(t20);
t26 = cos(qJ(1));
t45 = t12 * t26;
t23 = sin(qJ(1));
t44 = t22 * t23;
t27 = -pkin(9) - pkin(8);
t41 = t25 * t27;
t13 = cos(t20);
t40 = t26 * t13;
t19 = pkin(6) + r_base(3);
t39 = t23 * pkin(1) + r_base(2);
t38 = pkin(2) + t19;
t37 = t23 * pkin(7) + t39;
t36 = t26 * pkin(1) + t23 * qJ(2) + r_base(1);
t34 = t26 * pkin(7) + t36;
t11 = pkin(4) * t24 + pkin(3);
t1 = (-m(1) * r_base(3) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + t53 * t38 + t52 * (t25 * t11 - t22 * t27 + t38) + (-m(2) - m(3)) * t19 + (-t48 * t12 + t49 * t13 + t56) * t25 + t55 * t22) * g(3) + (-m(3) * t39 - mrSges(1,2) + t54 * r_base(2) + t52 * (t23 * t46 + t37) + t49 * (t12 * t23 - t22 * t40) + t53 * t37 + t48 * (t13 * t23 + t22 * t45) + t50 * t23 + (t52 * (-t11 * t22 - qJ(2) - t41) + (m(3) - t53) * qJ(2) - t47) * t26) * g(2) + (-m(3) * t36 - mrSges(1,1) + t54 * r_base(1) + t53 * t34 + t52 * (t11 * t44 + t23 * t41 + t26 * t46 + t34) + t49 * (t13 * t44 + t45) - t48 * (t12 * t44 - t40) + t50 * t26 + t47 * t23) * g(1);
U  = t1;
