% Calculate potential energy for
% S6RRPRPR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR11_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRPR11_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR11_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR11_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR11_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:11:35
% EndTime: 2019-03-09 11:11:35
% DurationCPUTime: 0.59s
% Computational Cost: add. (207->85), mult. (251->73), div. (0->0), fcn. (221->10), ass. (0->37)
t21 = -qJ(5) - pkin(8);
t62 = -mrSges(3,1) + mrSges(4,2) + m(7) * (-pkin(9) + t21) - mrSges(7,3) + m(6) * t21 - mrSges(6,3);
t19 = qJ(4) + pkin(10);
t10 = sin(t19);
t11 = cos(t19);
t22 = sin(qJ(4));
t49 = pkin(4) * t22;
t12 = qJ(6) + t19;
t7 = sin(t12);
t8 = cos(t12);
t61 = -m(6) * t49 - t10 * mrSges(6,1) - t11 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3) - m(7) * (pkin(5) * t10 + t49) - t7 * mrSges(7,1) - t8 * mrSges(7,2);
t60 = -m(1) - m(2);
t24 = sin(qJ(1));
t23 = sin(qJ(2));
t42 = qJ(3) * t23;
t26 = cos(qJ(2));
t45 = t24 * t26;
t59 = pkin(2) * t45 + t24 * t42;
t58 = -m(5) - m(6) - m(7);
t57 = m(4) - t58;
t56 = -m(5) * pkin(8) - mrSges(5,3);
t51 = t61 * t23 + t62 * t26 - mrSges(2,1);
t50 = t11 * mrSges(6,1) + t8 * mrSges(7,1) - t10 * mrSges(6,2) - t7 * mrSges(7,2) + mrSges(4,1) - mrSges(2,2) + mrSges(3,3);
t25 = cos(qJ(4));
t9 = t25 * pkin(4) + pkin(3);
t27 = cos(qJ(1));
t48 = t22 * t27;
t47 = t24 * t22;
t46 = t24 * t25;
t44 = t25 * t27;
t43 = t26 * t27;
t20 = pkin(6) + r_base(3);
t41 = t24 * pkin(1) + r_base(2);
t39 = t27 * pkin(1) + t24 * pkin(7) + r_base(1);
t31 = -pkin(7) * t27 + t41;
t1 = pkin(5) * t11 + t9;
t2 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t20 - t57 * (t23 * pkin(2) + t20) + (t22 * mrSges(5,1) + t25 * mrSges(5,2) + t57 * qJ(3) - t61) * t26 + (t56 + t62) * t23) * g(3) + (-mrSges(1,2) - m(3) * t31 - m(4) * (t31 + t59) - (t23 * t47 - t44) * mrSges(5,1) - (t23 * t46 + t48) * mrSges(5,2) + t60 * r_base(2) + t56 * t45 + t58 * (t41 + t59) + (-m(5) * (-pkin(3) - pkin(7)) - m(6) * (-pkin(7) - t9) - m(7) * (-pkin(7) - t1) + t50) * t27 + t51 * t24) * g(2) + (-mrSges(1,1) - m(3) * t39 - (t23 * t48 + t46) * mrSges(5,1) - (t23 * t44 - t47) * mrSges(5,2) + t60 * r_base(1) + t56 * t43 - t57 * (pkin(2) * t43 + t27 * t42 + t39) + t51 * t27 + (-m(5) * pkin(3) - m(6) * t9 - m(7) * t1 - t50) * t24) * g(1);
U  = t2;
