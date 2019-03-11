% Calculate potential energy for
% S6RPRRPP2
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
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPRRPP2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:31:23
% EndTime: 2019-03-09 04:31:23
% DurationCPUTime: 0.47s
% Computational Cost: add. (253->79), mult. (248->74), div. (0->0), fcn. (228->8), ass. (0->37)
t55 = -m(1) - m(2);
t54 = -m(6) - m(7);
t53 = mrSges(3,2) - mrSges(4,3);
t24 = sin(qJ(4));
t25 = sin(qJ(3));
t27 = cos(qJ(4));
t52 = (pkin(4) * t27 + qJ(5) * t24) * t25;
t51 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t50 = -mrSges(5,3) - mrSges(6,2) + mrSges(7,3);
t28 = cos(qJ(3));
t49 = -t28 * mrSges(4,1) - mrSges(3,1) + (m(7) * qJ(6) + mrSges(4,2)) * t25;
t48 = -m(7) * pkin(5) - mrSges(5,1) - mrSges(6,1) - mrSges(7,1);
t47 = pkin(3) * t28;
t23 = qJ(1) + pkin(9);
t17 = sin(t23);
t46 = t17 * t25;
t18 = cos(t23);
t45 = t18 * t25;
t44 = t24 * t28;
t43 = t27 * t28;
t41 = pkin(6) + r_base(3);
t26 = sin(qJ(1));
t40 = t26 * pkin(1) + r_base(2);
t29 = cos(qJ(1));
t39 = t29 * pkin(1) + r_base(1);
t19 = qJ(2) + t41;
t38 = t25 * pkin(3) + t19;
t37 = t18 * pkin(2) + t17 * pkin(7) + t39;
t35 = t17 * pkin(2) - pkin(7) * t18 + t40;
t34 = pkin(8) * t45 + t18 * t47 + t37;
t33 = -pkin(8) * t28 + t38;
t32 = pkin(8) * t46 + t17 * t47 + t35;
t6 = t17 * t24 + t18 * t43;
t5 = -t17 * t27 + t18 * t44;
t4 = t17 * t43 - t18 * t24;
t3 = t17 * t44 + t18 * t27;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t41 - mrSges(2,3) - mrSges(3,3) - m(5) * t33 - m(6) * (t33 + t52) - m(7) * (t38 + t52) + (-m(3) - m(4)) * t19 + (-mrSges(4,2) - m(7) * (-pkin(8) + qJ(6)) - t50) * t28 + (t51 * t24 + t48 * t27 - mrSges(4,1)) * t25) * g(3) + (-m(3) * t40 - m(4) * t35 - m(5) * t32 - t26 * mrSges(2,1) - t29 * mrSges(2,2) - mrSges(1,2) + t55 * r_base(2) + t54 * (t4 * pkin(4) + qJ(5) * t3 + t32) - t53 * t18 + t49 * t17 + t50 * t46 + t48 * t4 + t51 * t3) * g(2) + (-m(3) * t39 - m(4) * t37 - m(5) * t34 - t29 * mrSges(2,1) + t26 * mrSges(2,2) - mrSges(1,1) + t55 * r_base(1) + t54 * (t6 * pkin(4) + qJ(5) * t5 + t34) + t49 * t18 + t53 * t17 + t48 * t6 + t51 * t5 + t50 * t45) * g(1);
U  = t1;
