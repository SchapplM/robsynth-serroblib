% Calculate potential energy for
% S6RPRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPP3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:35:07
% EndTime: 2019-03-09 04:35:07
% DurationCPUTime: 0.45s
% Computational Cost: add. (253->79), mult. (248->74), div. (0->0), fcn. (228->8), ass. (0->38)
t56 = -m(1) - m(2);
t55 = -m(6) - m(7);
t26 = sin(qJ(3));
t29 = cos(qJ(3));
t54 = -t29 * mrSges(4,1) + mrSges(4,2) * t26 - mrSges(3,1);
t53 = mrSges(3,2) - mrSges(4,3);
t25 = sin(qJ(4));
t28 = cos(qJ(4));
t52 = (pkin(4) * t28 + qJ(5) * t25) * t26;
t51 = mrSges(6,1) + mrSges(7,1) + mrSges(5,3);
t50 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t49 = -m(7) * pkin(5) - t51;
t48 = -m(7) * qJ(6) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t47 = pkin(3) * t29;
t24 = qJ(1) + pkin(9);
t18 = sin(t24);
t46 = t18 * t26;
t19 = cos(t24);
t45 = t19 * t26;
t44 = t25 * t29;
t43 = t28 * t29;
t42 = pkin(6) + r_base(3);
t27 = sin(qJ(1));
t41 = t27 * pkin(1) + r_base(2);
t30 = cos(qJ(1));
t40 = t30 * pkin(1) + r_base(1);
t20 = qJ(2) + t42;
t39 = t26 * pkin(3) + t20;
t38 = t19 * pkin(2) + t18 * pkin(7) + t40;
t36 = t18 * pkin(2) - pkin(7) * t19 + t41;
t35 = pkin(8) * t45 + t19 * t47 + t38;
t34 = -pkin(8) * t29 + t39;
t33 = pkin(8) * t46 + t18 * t47 + t36;
t6 = t18 * t25 + t19 * t43;
t5 = -t18 * t28 + t19 * t44;
t4 = t18 * t43 - t19 * t25;
t3 = t18 * t44 + t19 * t28;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t42 - mrSges(2,3) - mrSges(3,3) - m(5) * t34 - m(6) * (t34 + t52) - m(7) * (t39 + t52) + (-m(3) - m(4)) * t20 + (-mrSges(4,2) - m(7) * (-pkin(5) - pkin(8)) + t51) * t29 + (t50 * t25 + t48 * t28 - mrSges(4,1)) * t26) * g(3) + (-m(3) * t41 - m(4) * t36 - m(5) * t33 - t27 * mrSges(2,1) - t30 * mrSges(2,2) - mrSges(1,2) + t56 * r_base(2) + t55 * (t4 * pkin(4) + qJ(5) * t3 + t33) - t53 * t19 + t54 * t18 + t49 * t46 + t48 * t4 + t50 * t3) * g(2) + (-m(3) * t40 - m(4) * t38 - m(5) * t35 - t30 * mrSges(2,1) + t27 * mrSges(2,2) - mrSges(1,1) + t56 * r_base(1) + t55 * (t6 * pkin(4) + qJ(5) * t5 + t35) + t54 * t19 + t53 * t18 + t48 * t6 + t50 * t5 + t49 * t45) * g(1);
U  = t1;
