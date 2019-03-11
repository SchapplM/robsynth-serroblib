% Calculate potential energy for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP8_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPP8_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP8_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP8_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP8_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:53:53
% EndTime: 2019-03-09 04:53:53
% DurationCPUTime: 0.53s
% Computational Cost: add. (164->77), mult. (250->67), div. (0->0), fcn. (230->6), ass. (0->36)
t57 = mrSges(6,1) - mrSges(4,2) + mrSges(5,3);
t54 = -m(6) - m(7);
t56 = -m(5) + t54;
t55 = -m(1) - m(2);
t53 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t52 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t51 = -m(7) * pkin(5) - mrSges(7,1);
t22 = sin(qJ(3));
t25 = cos(qJ(3));
t50 = -t22 * mrSges(4,1) + t57 * t25 + mrSges(2,2) - mrSges(3,3);
t49 = m(7) * qJ(6) + mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t48 = pkin(3) * t22;
t21 = sin(qJ(4));
t26 = cos(qJ(1));
t47 = t21 * t26;
t23 = sin(qJ(1));
t46 = t23 * t21;
t24 = cos(qJ(4));
t45 = t23 * t24;
t44 = t23 * t25;
t41 = t26 * t24;
t20 = pkin(6) + r_base(3);
t40 = pkin(8) * t44;
t39 = t23 * pkin(1) + r_base(2);
t38 = pkin(2) + t20;
t36 = t23 * pkin(7) + t39;
t35 = t26 * pkin(1) + t23 * qJ(2) + r_base(1);
t34 = t26 * t25 * pkin(8) + t36;
t33 = t26 * pkin(7) + t35;
t30 = t23 * t48 + t33;
t3 = t22 * t46 - t41;
t4 = t22 * t45 + t47;
t27 = t4 * pkin(4) + qJ(5) * t3 + t30;
t6 = t22 * t41 - t46;
t5 = t22 * t47 + t45;
t1 = (-m(1) * r_base(3) - m(4) * t38 - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t20 + (t51 - t57) * t22 + t56 * (t25 * pkin(3) + t22 * pkin(8) + t38) + (t54 * (pkin(4) * t24 + qJ(5) * t21) + t52 * t21 - t49 * t24 - mrSges(4,1)) * t25) * g(3) + (-m(3) * t39 - m(4) * t36 - m(5) * t34 - mrSges(1,2) + t55 * r_base(2) + t54 * (-t6 * pkin(4) - t5 * qJ(5) + t34) + t49 * t6 - t52 * t5 + t53 * t23 + (t51 * t25 + t56 * (-qJ(2) - t48) + (m(3) + m(4)) * qJ(2) - t50) * t26) * g(2) + (-mrSges(1,1) - m(3) * t35 - m(4) * t33 - m(5) * (t30 - t40) - m(6) * (t27 - t40) - m(7) * t27 - (m(7) * (-pkin(5) - pkin(8)) - mrSges(7,1)) * t44 + t55 * r_base(1) - t49 * t4 + t52 * t3 + t53 * t26 + t50 * t23) * g(1);
U  = t1;
