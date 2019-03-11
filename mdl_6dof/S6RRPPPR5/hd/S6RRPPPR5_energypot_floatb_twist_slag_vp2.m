% Calculate potential energy for
% S6RRPPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
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
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPPR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPPR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR5_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR5_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:21:15
% EndTime: 2019-03-09 08:21:15
% DurationCPUTime: 0.57s
% Computational Cost: add. (202->87), mult. (349->81), div. (0->0), fcn. (351->8), ass. (0->43)
t62 = -m(1) - m(2);
t23 = sin(pkin(9));
t26 = sin(qJ(2));
t24 = cos(pkin(9));
t53 = t24 * t26;
t61 = t26 * t23 * qJ(4) + pkin(3) * t53;
t29 = cos(qJ(2));
t60 = -t29 * mrSges(3,1) + t26 * mrSges(3,2) - mrSges(2,1);
t59 = mrSges(2,2) - mrSges(3,3);
t58 = mrSges(5,1) + mrSges(4,3) + mrSges(7,3) - mrSges(6,2);
t57 = -m(7) * pkin(8) - t58;
t25 = sin(qJ(6));
t28 = cos(qJ(6));
t56 = -t25 * mrSges(7,1) - t28 * mrSges(7,2) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t55 = -t28 * mrSges(7,1) + t25 * mrSges(7,2) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t54 = -m(7) * (pkin(5) + qJ(4)) + t55;
t27 = sin(qJ(1));
t52 = t26 * t27;
t30 = cos(qJ(1));
t51 = t26 * t30;
t50 = t27 * t29;
t49 = t29 * t30;
t3 = t23 * t50 + t24 * t30;
t48 = t3 * qJ(4);
t47 = t30 * t23;
t5 = -t27 * t24 + t29 * t47;
t46 = t5 * qJ(4);
t45 = -pkin(4) - qJ(3);
t43 = qJ(3) * t26;
t22 = pkin(6) + r_base(3);
t42 = t26 * pkin(2) + t22;
t41 = t30 * pkin(1) + t27 * pkin(7) + r_base(1);
t39 = t27 * pkin(1) - pkin(7) * t30 + r_base(2);
t38 = pkin(2) * t49 + t30 * t43 + t41;
t36 = -qJ(3) * t29 + t42;
t6 = t27 * t23 + t24 * t49;
t35 = t6 * pkin(3) + t38;
t34 = pkin(2) * t50 + t27 * t43 + t39;
t4 = t24 * t50 - t47;
t33 = t4 * pkin(3) + t34;
t32 = pkin(4) * t51 + t6 * qJ(5) + t35;
t31 = pkin(4) * t52 + t4 * qJ(5) + t33;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - m(4) * t36 - m(5) * (t36 + t61) + (-m(6) - m(7)) * (qJ(5) * t53 + t42 + t61) + (-m(2) - m(3)) * t22 + (-mrSges(3,2) - m(6) * t45 - m(7) * (-pkin(8) + t45) + t58) * t29 + (-mrSges(3,1) + t56 * t24 + (-m(7) * pkin(5) + t55) * t23) * t26) * g(3) + (-mrSges(1,2) - m(3) * t39 - m(4) * t34 - m(5) * (t33 + t48) - m(6) * (t31 + t48) - m(7) * t31 + t62 * r_base(2) - t59 * t30 + t60 * t27 + t56 * t4 + t54 * t3 + t57 * t52) * g(2) + (-mrSges(1,1) - m(3) * t41 - m(4) * t38 - m(5) * (t35 + t46) - m(6) * (t32 + t46) - m(7) * t32 + t62 * r_base(1) + t60 * t30 + t59 * t27 + t56 * t6 + t54 * t5 + t57 * t51) * g(1);
U  = t1;
