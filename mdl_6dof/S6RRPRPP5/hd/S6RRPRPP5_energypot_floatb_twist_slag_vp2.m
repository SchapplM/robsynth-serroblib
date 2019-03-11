% Calculate potential energy for
% S6RRPRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
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
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPP5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRPP5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP5_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP5_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:03:24
% EndTime: 2019-03-09 10:03:25
% DurationCPUTime: 0.59s
% Computational Cost: add. (176->78), mult. (281->70), div. (0->0), fcn. (261->6), ass. (0->38)
t60 = mrSges(3,2) - mrSges(4,3);
t59 = m(7) * qJ(6) - mrSges(3,1) + mrSges(4,2);
t57 = -m(1) - m(2);
t56 = -m(6) - m(7);
t24 = sin(qJ(1));
t23 = sin(qJ(2));
t42 = qJ(3) * t23;
t26 = cos(qJ(2));
t45 = t24 * t26;
t55 = pkin(2) * t45 + t24 * t42;
t54 = mrSges(5,1) + mrSges(6,1) + mrSges(7,1);
t53 = mrSges(2,2) - mrSges(3,3) - mrSges(4,1);
t52 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t51 = -mrSges(5,3) - mrSges(6,2) + mrSges(7,3);
t50 = -m(7) * pkin(5) - t54;
t49 = t60 * t23 + t59 * t26 - mrSges(2,1);
t22 = sin(qJ(4));
t48 = t22 * t24;
t27 = cos(qJ(1));
t47 = t22 * t27;
t25 = cos(qJ(4));
t46 = t24 * t25;
t44 = t25 * t27;
t43 = t26 * t27;
t21 = pkin(6) + r_base(3);
t40 = t24 * pkin(1) + r_base(2);
t39 = t23 * pkin(2) + t21;
t38 = t27 * pkin(1) + t24 * pkin(7) + r_base(1);
t37 = t23 * pkin(8) + t39;
t34 = -pkin(7) * t27 + t40;
t32 = pkin(2) * t43 + t27 * t42 + t38;
t31 = t24 * pkin(3) + pkin(8) * t43 + t32;
t30 = pkin(8) * t45 + (-pkin(3) - pkin(7)) * t27 + t40 + t55;
t6 = t23 * t48 - t44;
t5 = t23 * t46 + t47;
t4 = t23 * t47 + t46;
t3 = -t23 * t44 + t48;
t1 = (-m(1) * r_base(3) - m(4) * t39 - m(5) * t37 - mrSges(1,3) - mrSges(2,3) + t56 * (t26 * t25 * qJ(5) + t37) + (-m(2) - m(3)) * t21 + (t52 * t25 + (m(4) + m(5) - t56) * qJ(3) + (m(6) * pkin(4) - m(7) * (-pkin(4) - pkin(5)) + t54) * t22 - t60) * t26 + (t51 + t59) * t23) * g(3) + (-mrSges(1,2) - m(3) * t34 - m(4) * (t34 + t55) - m(5) * t30 + t57 * r_base(2) + t56 * (t6 * pkin(4) - t5 * qJ(5) + t30) + t50 * t6 - t52 * t5 + t51 * t45 - t53 * t27 + t49 * t24) * g(2) + (-m(3) * t38 - m(4) * t32 - m(5) * t31 - mrSges(1,1) + t57 * r_base(1) + t56 * (t4 * pkin(4) + t3 * qJ(5) + t31) + t51 * t43 + t50 * t4 + t52 * t3 + t49 * t27 + t53 * t24) * g(1);
U  = t1;
