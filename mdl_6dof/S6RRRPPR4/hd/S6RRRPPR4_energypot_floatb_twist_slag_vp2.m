% Calculate potential energy for
% S6RRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPPR4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:31:39
% EndTime: 2019-03-09 15:31:40
% DurationCPUTime: 0.70s
% Computational Cost: add. (265->91), mult. (320->81), div. (0->0), fcn. (314->10), ass. (0->43)
t25 = sin(qJ(3));
t29 = cos(qJ(3));
t64 = -m(4) * pkin(2) - mrSges(4,1) * t29 + mrSges(4,2) * t25 - mrSges(3,1);
t63 = -m(4) * pkin(8) + mrSges(3,2) - mrSges(4,3) + mrSges(7,3);
t62 = -m(6) - m(7);
t61 = -m(1) - m(2);
t60 = -m(3) - m(4);
t21 = qJ(3) + pkin(10);
t16 = sin(t21);
t17 = cos(t21);
t31 = cos(qJ(1));
t27 = sin(qJ(1));
t30 = cos(qJ(2));
t48 = t27 * t30;
t3 = t16 * t48 + t17 * t31;
t4 = -t16 * t31 + t17 * t48;
t59 = t4 * pkin(4) + t3 * qJ(5);
t46 = t30 * t31;
t5 = t16 * t46 - t27 * t17;
t6 = t27 * t16 + t17 * t46;
t58 = t6 * pkin(4) + t5 * qJ(5);
t57 = -mrSges(5,3) - mrSges(6,2);
t56 = -t25 * mrSges(4,1) - t29 * mrSges(4,2) + mrSges(2,2) - mrSges(3,3);
t23 = -qJ(4) - pkin(8);
t26 = sin(qJ(2));
t54 = -mrSges(2,1) + t64 * t30 + (-m(7) * (-pkin(9) - t23) + t63) * t26;
t24 = sin(qJ(6));
t28 = cos(qJ(6));
t53 = -t24 * mrSges(7,1) - t28 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t52 = -m(7) * pkin(5) - t28 * mrSges(7,1) + t24 * mrSges(7,2) - mrSges(5,1) - mrSges(6,1);
t51 = pkin(3) * t25;
t50 = t26 * t27;
t49 = t26 * t31;
t22 = pkin(6) + r_base(3);
t44 = t27 * pkin(1) + r_base(2);
t43 = t31 * pkin(1) + t27 * pkin(7) + r_base(1);
t39 = -pkin(7) * t31 + t44;
t15 = pkin(3) * t29 + pkin(2);
t38 = t15 * t46 + t27 * t51 + t43;
t34 = -t23 * t49 + t38;
t7 = t15 * t48;
t33 = -t23 * t50 + t7 + (-pkin(7) - t51) * t31 + t44;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(2) + t60) * t22 + (-m(7) * pkin(9) - t57 - t63) * t30 + (-m(5) + t62) * (t26 * t15 + t30 * t23 + t22) + (t62 * (pkin(4) * t17 + qJ(5) * t16) + t53 * t16 + t52 * t17 + t64) * t26) * g(3) + (-mrSges(1,2) - m(3) * t39 - m(4) * t44 - m(5) * t33 - m(6) * (t33 + t59) - m(7) * (t39 + t7 + t59) + t61 * r_base(2) + t57 * t50 + t52 * t4 + t53 * t3 + (m(4) * pkin(7) + m(7) * t51 - t56) * t31 + t54 * t27) * g(2) + (-mrSges(1,1) - m(5) * t34 - m(6) * (t34 + t58) - m(7) * (t38 + t58) + t61 * r_base(1) + t52 * t6 + t53 * t5 + t57 * t49 + t60 * t43 + t56 * t27 + t54 * t31) * g(1);
U  = t1;
