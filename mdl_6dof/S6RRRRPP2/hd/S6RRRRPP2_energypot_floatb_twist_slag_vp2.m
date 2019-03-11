% Calculate potential energy for
% S6RRRRPP2
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
% Datum: 2019-03-09 20:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPP2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:48:51
% EndTime: 2019-03-09 20:48:52
% DurationCPUTime: 0.48s
% Computational Cost: add. (244->79), mult. (262->72), div. (0->0), fcn. (242->8), ass. (0->40)
t58 = -m(2) - m(3);
t57 = -m(6) - m(7);
t23 = qJ(2) + qJ(3);
t19 = sin(t23);
t24 = sin(qJ(4));
t27 = cos(qJ(4));
t56 = (pkin(4) * t27 + qJ(5) * t24) * t19;
t55 = -m(1) + t58;
t54 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t53 = -mrSges(5,3) - mrSges(6,2) + mrSges(7,3);
t52 = m(3) * pkin(7) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3);
t51 = -m(7) * pkin(5) - mrSges(5,1) - mrSges(6,1) - mrSges(7,1);
t20 = cos(t23);
t25 = sin(qJ(2));
t28 = cos(qJ(2));
t50 = -m(3) * pkin(1) - t28 * mrSges(3,1) - t20 * mrSges(4,1) + t25 * mrSges(3,2) - mrSges(2,1) + (m(7) * qJ(6) + mrSges(4,2)) * t19;
t49 = pkin(3) * t20;
t26 = sin(qJ(1));
t48 = t19 * t26;
t29 = cos(qJ(1));
t47 = t19 * t29;
t46 = t24 * t26;
t45 = t24 * t29;
t44 = t26 * t27;
t43 = t27 * t29;
t22 = pkin(6) + r_base(3);
t41 = t25 * pkin(2) + t22;
t17 = pkin(2) * t28 + pkin(1);
t30 = -pkin(8) - pkin(7);
t40 = t26 * t17 + t29 * t30 + r_base(2);
t39 = t19 * pkin(3) + t41;
t37 = t29 * t17 - t26 * t30 + r_base(1);
t36 = pkin(9) * t48 + t26 * t49 + t40;
t35 = -pkin(9) * t20 + t39;
t34 = pkin(9) * t47 + t29 * t49 + t37;
t6 = t20 * t43 + t46;
t5 = t20 * t45 - t44;
t4 = t20 * t44 - t45;
t3 = t20 * t46 + t43;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - t25 * mrSges(3,1) - t28 * mrSges(3,2) - m(4) * t41 - m(5) * t35 - m(6) * (t35 + t56) - m(7) * (t39 + t56) + t58 * t22 + (-mrSges(4,2) - m(7) * (-pkin(9) + qJ(6)) - t53) * t20 + (t54 * t24 + t51 * t27 - mrSges(4,1)) * t19) * g(3) + (-m(4) * t40 - m(5) * t36 - mrSges(1,2) + t57 * (t4 * pkin(4) + qJ(5) * t3 + t36) + t55 * r_base(2) + t53 * t48 + t51 * t4 + t54 * t3 + t52 * t29 + t50 * t26) * g(2) + (-m(4) * t37 - m(5) * t34 - mrSges(1,1) + t57 * (t6 * pkin(4) + t5 * qJ(5) + t34) + t55 * r_base(1) + t51 * t6 + t54 * t5 + t53 * t47 + t50 * t29 - t52 * t26) * g(1);
U  = t1;
