% Calculate potential energy for
% S6RPRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRPRR4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:44:19
% EndTime: 2019-03-09 03:44:20
% DurationCPUTime: 0.56s
% Computational Cost: add. (232->77), mult. (204->62), div. (0->0), fcn. (170->10), ass. (0->33)
t58 = -mrSges(4,1) + mrSges(5,2) - mrSges(6,3) + m(7) * (-pkin(9) - pkin(8)) - mrSges(7,3);
t18 = qJ(5) + qJ(6);
t12 = sin(t18);
t13 = cos(t18);
t19 = sin(qJ(5));
t22 = cos(qJ(5));
t57 = -t12 * mrSges(7,1) - t22 * mrSges(6,2) - t13 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t19;
t54 = m(6) * pkin(8);
t53 = -m(1) - m(2);
t20 = sin(qJ(3));
t40 = qJ(4) * t20;
t23 = cos(qJ(3));
t17 = qJ(1) + pkin(10);
t9 = sin(t17);
t45 = t23 * t9;
t52 = pkin(3) * t45 + t40 * t9;
t51 = m(5) + m(6) + m(7);
t48 = t22 * mrSges(6,1) + t13 * mrSges(7,1) - t19 * mrSges(6,2) - t12 * mrSges(7,2) + mrSges(5,1) - mrSges(3,2) + mrSges(4,3);
t47 = t20 * t57 + t23 * t58 - mrSges(3,1);
t10 = cos(t17);
t44 = t10 * t23;
t39 = pkin(6) + r_base(3);
t21 = sin(qJ(1));
t38 = pkin(1) * t21 + r_base(2);
t24 = cos(qJ(1));
t37 = pkin(1) * t24 + r_base(1);
t36 = pkin(2) * t9 + t38;
t11 = qJ(2) + t39;
t35 = pkin(2) * t10 + pkin(7) * t9 + t37;
t33 = t36 + t52;
t27 = -t10 * pkin(7) + t36;
t8 = pkin(5) * t22 + pkin(4);
t1 = (-m(1) * r_base(3) - m(2) * t39 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + (-m(3) - m(4)) * t11 - t51 * (pkin(3) * t20 + t11) + (qJ(4) * t51 - t57) * t23 + (-t54 + t58) * t20) * g(3) + (-mrSges(1,2) - t21 * mrSges(2,1) - t24 * mrSges(2,2) - m(3) * t38 - m(4) * t27 - m(5) * (t27 + t52) - m(6) * (pkin(8) * t45 + t33) - m(7) * t33 + t53 * r_base(2) + (-m(6) * (-pkin(4) - pkin(7)) - m(7) * (-pkin(7) - t8) + t48) * t10 + t47 * t9) * g(2) + (-t44 * t54 - m(3) * t37 - m(4) * t35 - t24 * mrSges(2,1) + t21 * mrSges(2,2) - mrSges(1,1) + t53 * r_base(1) - t51 * (pkin(3) * t44 + t10 * t40 + t35) + (-m(6) * pkin(4) - m(7) * t8 - t48) * t9 + t47 * t10) * g(1);
U  = t1;
