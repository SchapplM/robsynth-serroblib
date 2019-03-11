% Calculate potential energy for
% S6RRPRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPP2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRPP2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:49:48
% EndTime: 2019-03-09 09:49:49
% DurationCPUTime: 0.45s
% Computational Cost: add. (244->79), mult. (262->72), div. (0->0), fcn. (242->8), ass. (0->40)
t58 = -m(2) - m(3);
t57 = -m(6) - m(7);
t22 = qJ(2) + pkin(9);
t19 = sin(t22);
t25 = sin(qJ(4));
t28 = cos(qJ(4));
t56 = (pkin(4) * t28 + qJ(5) * t25) * t19;
t55 = -m(1) + t58;
t54 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t53 = -mrSges(5,3) - mrSges(6,2) + mrSges(7,3);
t52 = m(3) * pkin(7) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3);
t51 = -m(7) * pkin(5) - mrSges(5,1) - mrSges(6,1) - mrSges(7,1);
t20 = cos(t22);
t26 = sin(qJ(2));
t29 = cos(qJ(2));
t50 = -m(3) * pkin(1) - t29 * mrSges(3,1) - t20 * mrSges(4,1) + t26 * mrSges(3,2) - mrSges(2,1) + (m(7) * qJ(6) + mrSges(4,2)) * t19;
t49 = pkin(3) * t20;
t27 = sin(qJ(1));
t48 = t19 * t27;
t30 = cos(qJ(1));
t47 = t19 * t30;
t46 = t25 * t30;
t45 = t27 * t25;
t44 = t27 * t28;
t43 = t28 * t30;
t23 = pkin(6) + r_base(3);
t41 = t26 * pkin(2) + t23;
t18 = pkin(2) * t29 + pkin(1);
t24 = -qJ(3) - pkin(7);
t40 = t27 * t18 + t30 * t24 + r_base(2);
t39 = t19 * pkin(3) + t41;
t37 = t30 * t18 - t27 * t24 + r_base(1);
t36 = pkin(8) * t48 + t27 * t49 + t40;
t35 = -pkin(8) * t20 + t39;
t34 = pkin(8) * t47 + t30 * t49 + t37;
t6 = t20 * t43 + t45;
t5 = t20 * t46 - t44;
t4 = t20 * t44 - t46;
t3 = t20 * t45 + t43;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - t26 * mrSges(3,1) - t29 * mrSges(3,2) - m(4) * t41 - m(5) * t35 - m(6) * (t35 + t56) - m(7) * (t39 + t56) + t58 * t23 + (-mrSges(4,2) - m(7) * (-pkin(8) + qJ(6)) - t53) * t20 + (t54 * t25 + t51 * t28 - mrSges(4,1)) * t19) * g(3) + (-m(4) * t40 - m(5) * t36 - mrSges(1,2) + t57 * (t4 * pkin(4) + qJ(5) * t3 + t36) + t55 * r_base(2) + t53 * t48 + t51 * t4 + t52 * t30 + t54 * t3 + t50 * t27) * g(2) + (-m(4) * t37 - m(5) * t34 - mrSges(1,1) + t57 * (t6 * pkin(4) + t5 * qJ(5) + t34) + t55 * r_base(1) + t51 * t6 + t54 * t5 + t53 * t47 + t50 * t30 - t52 * t27) * g(1);
U  = t1;
