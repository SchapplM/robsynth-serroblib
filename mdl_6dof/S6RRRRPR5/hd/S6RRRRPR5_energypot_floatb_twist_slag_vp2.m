% Calculate potential energy for
% S6RRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR5_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR5_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:11:55
% EndTime: 2019-03-09 22:11:56
% DurationCPUTime: 0.60s
% Computational Cost: add. (260->82), mult. (288->76), div. (0->0), fcn. (278->10), ass. (0->43)
t60 = -m(2) - m(3);
t59 = -m(6) - m(7);
t23 = qJ(2) + qJ(3);
t19 = sin(t23);
t25 = sin(qJ(4));
t29 = cos(qJ(4));
t58 = (pkin(4) * t29 + qJ(5) * t25) * t19;
t57 = -m(1) + t60;
t20 = cos(t23);
t26 = sin(qJ(2));
t30 = cos(qJ(2));
t56 = -m(3) * pkin(1) - t30 * mrSges(3,1) - t20 * mrSges(4,1) + t26 * mrSges(3,2) + t19 * mrSges(4,2) - mrSges(2,1);
t55 = mrSges(6,2) + mrSges(5,3) - mrSges(7,3);
t54 = m(3) * pkin(7) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3);
t53 = m(7) * pkin(10) - t55;
t24 = sin(qJ(6));
t28 = cos(qJ(6));
t52 = -t24 * mrSges(7,1) - t28 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t51 = -m(7) * pkin(5) - t28 * mrSges(7,1) + t24 * mrSges(7,2) - mrSges(5,1) - mrSges(6,1);
t50 = pkin(3) * t20;
t27 = sin(qJ(1));
t49 = t19 * t27;
t31 = cos(qJ(1));
t48 = t19 * t31;
t47 = t25 * t27;
t46 = t25 * t31;
t45 = t27 * t29;
t44 = t29 * t31;
t22 = pkin(6) + r_base(3);
t43 = t26 * pkin(2) + t22;
t17 = pkin(2) * t30 + pkin(1);
t32 = -pkin(8) - pkin(7);
t42 = t27 * t17 + t31 * t32 + r_base(2);
t41 = t19 * pkin(3) + t43;
t39 = t31 * t17 - t27 * t32 + r_base(1);
t38 = pkin(9) * t49 + t27 * t50 + t42;
t37 = -pkin(9) * t20 + t41;
t36 = pkin(9) * t48 + t31 * t50 + t39;
t6 = t20 * t44 + t47;
t5 = t20 * t46 - t45;
t4 = t20 * t45 - t46;
t3 = t20 * t47 + t44;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - t26 * mrSges(3,1) - t30 * mrSges(3,2) - m(4) * t43 - m(5) * t37 - m(6) * (t37 + t58) - m(7) * (t41 + t58) + t60 * t22 + (-mrSges(4,2) - m(7) * (-pkin(9) + pkin(10)) + t55) * t20 + (t52 * t25 + t51 * t29 - mrSges(4,1)) * t19) * g(3) + (-m(4) * t42 - m(5) * t38 - mrSges(1,2) + t59 * (t4 * pkin(4) + qJ(5) * t3 + t38) + t51 * t4 + t52 * t3 + t57 * r_base(2) + t53 * t49 + t54 * t31 + t56 * t27) * g(2) + (-m(4) * t39 - m(5) * t36 - mrSges(1,1) + t59 * (t6 * pkin(4) + t5 * qJ(5) + t36) + t51 * t6 + t52 * t5 + t57 * r_base(1) + t53 * t48 + t56 * t31 - t54 * t27) * g(1);
U  = t1;
