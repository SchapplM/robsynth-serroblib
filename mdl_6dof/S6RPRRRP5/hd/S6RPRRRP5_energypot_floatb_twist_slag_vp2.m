% Calculate potential energy for
% S6RPRRRP5
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
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRRP5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP5_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP5_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:10:30
% EndTime: 2019-03-09 06:10:30
% DurationCPUTime: 0.44s
% Computational Cost: add. (260->77), mult. (213->69), div. (0->0), fcn. (183->10), ass. (0->37)
t58 = -m(2) - m(3);
t57 = -m(6) - m(7);
t56 = -mrSges(6,3) - mrSges(7,2);
t55 = -m(1) - m(4) + t58;
t25 = pkin(10) + qJ(3);
t21 = qJ(4) + t25;
t15 = sin(t21);
t16 = cos(t21);
t28 = cos(pkin(10));
t17 = t28 * pkin(2) + pkin(1);
t19 = sin(t25);
t20 = cos(t25);
t27 = sin(pkin(10));
t54 = -m(3) * pkin(1) - m(4) * t17 - t28 * mrSges(3,1) - t20 * mrSges(4,1) - t16 * mrSges(5,1) + t27 * mrSges(3,2) + t19 * mrSges(4,2) + t15 * mrSges(5,2) - mrSges(2,1);
t53 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t52 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t29 = -pkin(7) - qJ(2);
t51 = m(3) * qJ(2) - m(4) * t29 - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t50 = pkin(4) * t16;
t31 = sin(qJ(1));
t49 = t31 * t15;
t30 = sin(qJ(5));
t48 = t31 * t30;
t32 = cos(qJ(5));
t47 = t31 * t32;
t33 = cos(qJ(1));
t46 = t33 * t15;
t45 = t33 * t30;
t44 = t33 * t32;
t26 = pkin(6) + r_base(3);
t43 = t27 * pkin(2) + t26;
t24 = -pkin(8) + t29;
t7 = pkin(3) * t20 + t17;
t42 = t33 * t24 + t31 * t7 + r_base(2);
t41 = pkin(3) * t19 + t43;
t38 = -t31 * t24 + t33 * t7 + r_base(1);
t1 = (-m(1) * r_base(3) - m(4) * t43 - m(5) * t41 - t27 * mrSges(3,1) - t19 * mrSges(4,1) - t28 * mrSges(3,2) - t20 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) + t57 * (t15 * pkin(4) - t16 * pkin(9) + t41) + t58 * t26 + (-mrSges(5,2) - t56) * t16 + (t52 * t30 + t53 * t32 - mrSges(5,1)) * t15) * g(3) + (-m(5) * t42 - mrSges(1,2) + t56 * t49 + t57 * (pkin(9) * t49 + t31 * t50 + t42) + t53 * (t16 * t47 - t45) + t52 * (t16 * t48 + t44) + t55 * r_base(2) + t51 * t33 + t54 * t31) * g(2) + (-m(5) * t38 - mrSges(1,1) + t56 * t46 + t57 * (pkin(9) * t46 + t33 * t50 + t38) + t53 * (t16 * t44 + t48) + t52 * (t16 * t45 - t47) + t55 * r_base(1) + t54 * t33 - t51 * t31) * g(1);
U  = t1;
