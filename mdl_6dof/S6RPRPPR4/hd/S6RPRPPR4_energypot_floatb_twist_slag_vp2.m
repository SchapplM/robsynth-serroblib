% Calculate potential energy for
% S6RPRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPPR4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPPR4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:46:40
% EndTime: 2019-03-09 02:46:41
% DurationCPUTime: 0.63s
% Computational Cost: add. (260->78), mult. (288->71), div. (0->0), fcn. (278->10), ass. (0->40)
t64 = -mrSges(4,2) + mrSges(6,2) + mrSges(5,3) - mrSges(7,3);
t22 = pkin(9) + qJ(3);
t19 = sin(t22);
t20 = cos(t22);
t63 = pkin(3) * t20 + qJ(4) * t19;
t25 = sin(pkin(9));
t27 = cos(pkin(9));
t62 = -m(3) * pkin(1) - t27 * mrSges(3,1) - t20 * mrSges(4,1) + t25 * mrSges(3,2) - mrSges(2,1) + (m(7) * pkin(8) - t64) * t19;
t61 = -m(2) - m(3);
t60 = -m(6) - m(7);
t24 = sin(pkin(10));
t26 = cos(pkin(10));
t59 = (pkin(4) * t26 + qJ(5) * t24) * t19;
t58 = -m(1) + t61;
t55 = m(3) * qJ(2) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3);
t29 = sin(qJ(6));
t31 = cos(qJ(6));
t53 = -t29 * mrSges(7,1) - t31 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t52 = -m(7) * pkin(5) - t31 * mrSges(7,1) + t29 * mrSges(7,2) - mrSges(5,1) - mrSges(6,1);
t32 = cos(qJ(1));
t48 = t24 * t32;
t47 = t26 * t32;
t30 = sin(qJ(1));
t46 = t30 * t24;
t45 = t30 * t26;
t23 = pkin(6) + r_base(3);
t43 = t25 * pkin(2) + t23;
t17 = pkin(2) * t27 + pkin(1);
t28 = -pkin(7) - qJ(2);
t42 = t30 * t17 + t32 * t28 + r_base(2);
t41 = t19 * pkin(3) + t43;
t39 = t32 * t17 - t30 * t28 + r_base(1);
t38 = t63 * t30 + t42;
t37 = -qJ(4) * t20 + t41;
t36 = t63 * t32 + t39;
t6 = t20 * t47 + t46;
t5 = t20 * t48 - t45;
t4 = t20 * t45 - t48;
t3 = t20 * t46 + t47;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - t25 * mrSges(3,1) - t27 * mrSges(3,2) - m(4) * t43 - m(5) * t37 - m(6) * (t37 + t59) - m(7) * (t41 + t59) + t61 * t23 + (-m(7) * (pkin(8) - qJ(4)) + t64) * t20 + (t53 * t24 + t52 * t26 - mrSges(4,1)) * t19) * g(3) + (-m(4) * t42 - m(5) * t38 - mrSges(1,2) + t60 * (t4 * pkin(4) + qJ(5) * t3 + t38) + t52 * t4 + t53 * t3 + t58 * r_base(2) + t55 * t32 + t62 * t30) * g(2) + (-m(4) * t39 - m(5) * t36 - mrSges(1,1) + t60 * (t6 * pkin(4) + t5 * qJ(5) + t36) + t52 * t6 + t53 * t5 + t58 * r_base(1) - t55 * t30 + t62 * t32) * g(1);
U  = t1;
