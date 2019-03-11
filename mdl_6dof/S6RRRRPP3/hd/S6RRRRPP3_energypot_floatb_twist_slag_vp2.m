% Calculate potential energy for
% S6RRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPP3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:53:38
% EndTime: 2019-03-09 20:53:38
% DurationCPUTime: 0.47s
% Computational Cost: add. (244->79), mult. (262->72), div. (0->0), fcn. (242->8), ass. (0->41)
t59 = -m(2) - m(3);
t58 = -m(6) - m(7);
t24 = qJ(2) + qJ(3);
t20 = sin(t24);
t25 = sin(qJ(4));
t28 = cos(qJ(4));
t57 = (pkin(4) * t28 + qJ(5) * t25) * t20;
t56 = -m(1) + t59;
t21 = cos(t24);
t26 = sin(qJ(2));
t29 = cos(qJ(2));
t55 = -m(3) * pkin(1) - t29 * mrSges(3,1) - t21 * mrSges(4,1) + t26 * mrSges(3,2) + mrSges(4,2) * t20 - mrSges(2,1);
t54 = mrSges(6,1) + mrSges(7,1) + mrSges(5,3);
t53 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t52 = m(3) * pkin(7) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3);
t51 = -m(7) * pkin(5) - t54;
t50 = -m(7) * qJ(6) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t49 = pkin(3) * t21;
t27 = sin(qJ(1));
t48 = t20 * t27;
t30 = cos(qJ(1));
t47 = t20 * t30;
t46 = t25 * t27;
t45 = t27 * t28;
t44 = t28 * t30;
t43 = t30 * t25;
t23 = pkin(6) + r_base(3);
t42 = t26 * pkin(2) + t23;
t18 = pkin(2) * t29 + pkin(1);
t31 = -pkin(8) - pkin(7);
t41 = t27 * t18 + t30 * t31 + r_base(2);
t40 = t20 * pkin(3) + t42;
t38 = t30 * t18 - t27 * t31 + r_base(1);
t37 = pkin(9) * t48 + t27 * t49 + t41;
t36 = -pkin(9) * t21 + t40;
t35 = pkin(9) * t47 + t30 * t49 + t38;
t6 = t21 * t44 + t46;
t5 = t21 * t43 - t45;
t4 = t21 * t45 - t43;
t3 = t21 * t46 + t44;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - t26 * mrSges(3,1) - t29 * mrSges(3,2) - m(4) * t42 - m(5) * t36 - m(6) * (t36 + t57) - m(7) * (t40 + t57) + t59 * t23 + (-mrSges(4,2) - m(7) * (-pkin(5) - pkin(9)) + t54) * t21 + (t53 * t25 + t50 * t28 - mrSges(4,1)) * t20) * g(3) + (-m(4) * t41 - m(5) * t37 - mrSges(1,2) + t58 * (t4 * pkin(4) + qJ(5) * t3 + t37) + t56 * r_base(2) + t51 * t48 + t50 * t4 + t52 * t30 + t53 * t3 + t55 * t27) * g(2) + (-m(4) * t38 - m(5) * t35 - mrSges(1,1) + t58 * (t6 * pkin(4) + t5 * qJ(5) + t35) + t56 * r_base(1) + t50 * t6 + t53 * t5 + t51 * t47 + t55 * t30 - t52 * t27) * g(1);
U  = t1;
