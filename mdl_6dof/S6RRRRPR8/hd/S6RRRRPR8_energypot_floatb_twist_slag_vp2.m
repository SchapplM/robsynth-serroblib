% Calculate potential energy for
% S6RRRRPR8
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
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR8_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPR8_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR8_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR8_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR8_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:38:44
% EndTime: 2019-03-09 22:38:44
% DurationCPUTime: 0.74s
% Computational Cost: add. (265->89), mult. (320->78), div. (0->0), fcn. (314->10), ass. (0->41)
t68 = -m(4) * pkin(8) + mrSges(3,2) - mrSges(6,2) - mrSges(4,3) - mrSges(5,3) + mrSges(7,3);
t24 = sin(qJ(3));
t28 = cos(qJ(3));
t67 = -m(4) * pkin(2) - mrSges(4,1) * t28 + mrSges(4,2) * t24 - mrSges(3,1);
t65 = -m(6) - m(7);
t25 = sin(qJ(2));
t29 = cos(qJ(2));
t31 = -pkin(9) - pkin(8);
t64 = t67 * t29 - mrSges(2,1) + (-m(7) * (-pkin(10) - t31) + t68) * t25;
t63 = -m(1) - m(2);
t62 = -m(3) - m(4);
t22 = qJ(3) + qJ(4);
t16 = sin(t22);
t17 = cos(t22);
t30 = cos(qJ(1));
t26 = sin(qJ(1));
t48 = t26 * t29;
t3 = t16 * t48 + t17 * t30;
t4 = -t16 * t30 + t17 * t48;
t61 = t4 * pkin(4) + t3 * qJ(5);
t47 = t29 * t30;
t5 = t16 * t47 - t26 * t17;
t6 = t16 * t26 + t17 * t47;
t60 = t6 * pkin(4) + t5 * qJ(5);
t58 = -t24 * mrSges(4,1) - t28 * mrSges(4,2) + mrSges(2,2) - mrSges(3,3);
t23 = sin(qJ(6));
t27 = cos(qJ(6));
t54 = -t23 * mrSges(7,1) - t27 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t53 = -m(7) * pkin(5) - t27 * mrSges(7,1) + t23 * mrSges(7,2) - mrSges(5,1) - mrSges(6,1);
t52 = pkin(3) * t24;
t49 = t25 * t31;
t21 = pkin(6) + r_base(3);
t44 = t26 * pkin(1) + r_base(2);
t43 = t30 * pkin(1) + t26 * pkin(7) + r_base(1);
t39 = -t30 * pkin(7) + t44;
t14 = pkin(3) * t28 + pkin(2);
t38 = t14 * t47 + t26 * t52 + t43;
t34 = -t30 * t49 + t38;
t7 = t14 * t48;
t33 = -t26 * t49 + t7 + (-pkin(7) - t52) * t30 + t44;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(2) + t62) * t21 + (-m(7) * pkin(10) - t68) * t29 + (-m(5) + t65) * (t25 * t14 + t29 * t31 + t21) + (t65 * (pkin(4) * t17 + qJ(5) * t16) + t54 * t16 + t53 * t17 + t67) * t25) * g(3) + (-mrSges(1,2) - m(3) * t39 - m(4) * t44 - m(5) * t33 - m(6) * (t33 + t61) - m(7) * (t39 + t7 + t61) + t63 * r_base(2) + t53 * t4 + t54 * t3 + (m(4) * pkin(7) + m(7) * t52 - t58) * t30 + t64 * t26) * g(2) + (-mrSges(1,1) - m(5) * t34 - m(6) * (t34 + t60) - m(7) * (t38 + t60) + t63 * r_base(1) + t53 * t6 + t54 * t5 + t62 * t43 + t58 * t26 + t64 * t30) * g(1);
U  = t1;
