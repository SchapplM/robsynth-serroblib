% Calculate potential energy for
% S6RRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRR7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR7_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR7_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:55:35
% EndTime: 2019-03-09 13:55:35
% DurationCPUTime: 0.60s
% Computational Cost: add. (209->81), mult. (336->70), div. (0->0), fcn. (338->10), ass. (0->37)
t64 = -mrSges(3,1) - mrSges(4,1);
t63 = mrSges(3,2) - mrSges(4,3);
t27 = sin(qJ(4));
t55 = -m(6) * pkin(9) + m(7) * (-pkin(10) - pkin(9)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t62 = t27 * t55;
t28 = sin(qJ(2));
t32 = cos(qJ(2));
t31 = cos(qJ(4));
t5 = t27 * t28 + t31 * t32;
t50 = t28 * t31;
t25 = qJ(5) + qJ(6);
t16 = sin(t25);
t17 = cos(t25);
t26 = sin(qJ(5));
t30 = cos(qJ(5));
t53 = -m(6) * pkin(4) - m(7) * (pkin(5) * t30 + pkin(4)) - mrSges(6,1) * t30 - mrSges(7,1) * t17 + mrSges(6,2) * t26 + mrSges(7,2) * t16 - mrSges(5,1);
t61 = t63 * t28 + t64 * t32 + t5 * t53 - t50 * t55 - mrSges(2,1);
t60 = -m(1) - m(2);
t59 = -m(5) - m(6);
t29 = sin(qJ(1));
t47 = qJ(3) * t28;
t49 = t29 * t32;
t58 = pkin(2) * t49 + t29 * t47;
t33 = cos(qJ(1));
t57 = pkin(3) * t49 + t33 * pkin(8);
t54 = -t26 * mrSges(6,1) - t16 * mrSges(7,1) - t30 * mrSges(6,2) - t17 * mrSges(7,2) - mrSges(2,2) + mrSges(4,2) + mrSges(3,3) - mrSges(5,3);
t51 = pkin(5) * t26;
t48 = t32 * t33;
t24 = pkin(6) + r_base(3);
t46 = t29 * pkin(1) + r_base(2);
t45 = t33 * pkin(1) + t29 * pkin(7) + r_base(1);
t42 = -t33 * pkin(7) + t46;
t41 = pkin(2) * t48 + t33 * t47 + t45;
t40 = t28 * pkin(2) - t32 * qJ(3) + t24;
t39 = pkin(3) * t48 + t41;
t38 = t42 + t58;
t1 = (-m(1) * r_base(3) - m(4) * t40 - mrSges(1,3) - mrSges(2,3) - t63 * t32 + t64 * t28 + (-m(2) - m(3)) * t24 + t53 * (-t27 * t32 + t50) + (-m(7) + t59) * (t28 * pkin(3) + t40) + t55 * t5) * g(3) + (-mrSges(1,2) - m(3) * t42 - m(4) * t38 - m(7) * (t46 + t57 + t58) + t60 * r_base(2) + t59 * (t38 + t57) + t49 * t62 + (-m(7) * (-pkin(7) + t51) + t54) * t33 + t61 * t29) * g(2) + (-m(3) * t45 - m(4) * t41 - m(7) * t39 - mrSges(1,1) + t60 * r_base(1) + t59 * (-t29 * pkin(8) + t39) + t48 * t62 + (-m(7) * (-pkin(8) - t51) - t54) * t29 + t61 * t33) * g(1);
U  = t1;
