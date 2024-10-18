% Calculate potential energy for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR14_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRR14_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR14_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR14_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR14_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:42:17
% EndTime: 2024-09-27 18:42:17
% DurationCPUTime: 0.16s
% Computational Cost: add. (260->75), mult. (190->66), div. (0->0), fcn. (154->22), ass. (0->34)
t38 = sin(qJ(3));
t40 = cos(qJ(3));
t63 = (-m(5) * pkin(3) - mrSges(4,1)) * t38 - m(6) * (pkin(4) * sin(qJ(4)) * t40 + (cos(qJ(4)) * pkin(4) + pkin(3)) * t38) - t40 * mrSges(4,2);
t42 = pkin(8) + pkin(9);
t60 = m(5) * t42 + m(6) * (pkin(10) + t42) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3) + m(4) * pkin(8);
t57 = -m(1) - m(2);
t56 = -m(3) - m(4) - m(5) - m(6);
t21 = t40 * pkin(3) + pkin(2);
t33 = qJ(3) + qJ(4);
t26 = cos(t33);
t28 = qJ(5) + t33;
t55 = -m(4) * pkin(2) - m(5) * t21 - m(6) * (pkin(4) * t26 + t21) - t40 * mrSges(4,1) - t26 * mrSges(5,1) - cos(t28) * mrSges(6,1) + t38 * mrSges(4,2) + sin(t33) * mrSges(5,2) + sin(t28) * mrSges(6,2) - mrSges(3,1);
t22 = pkin(5) + t33;
t10 = sin(t22) / 0.2e1;
t23 = pkin(5) - t33;
t11 = cos(t23) / 0.2e1;
t17 = -qJ(5) + t23;
t12 = sin(t17);
t16 = qJ(5) + t22;
t13 = cos(t16);
t14 = sin(t23);
t15 = cos(t22);
t35 = sin(pkin(5));
t36 = cos(pkin(5));
t8 = sin(t16) / 0.2e1;
t9 = cos(t17) / 0.2e1;
t54 = -(t10 - t14 / 0.2e1) * mrSges(5,1) - (t8 - t12 / 0.2e1) * mrSges(6,1) - (t11 + t15 / 0.2e1) * mrSges(5,2) - (t9 + t13 / 0.2e1) * mrSges(6,2) - mrSges(3,2) + t60 * t35 + t63 * t36;
t47 = pkin(6) + r_base(3);
t41 = cos(qJ(1));
t39 = sin(qJ(1));
t34 = qJ(1) + qJ(2);
t27 = cos(t34);
t25 = sin(t34);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t47 - mrSges(2,3) - mrSges(3,3) - (t11 - t15 / 0.2e1) * mrSges(5,1) - (t10 + t14 / 0.2e1) * mrSges(5,2) - (t9 - t13 / 0.2e1) * mrSges(6,1) - (t8 + t12 / 0.2e1) * mrSges(6,2) + t63 * t35 + t56 * (pkin(7) + t47) - t60 * t36) * g(3) + (-t39 * mrSges(2,1) - t41 * mrSges(2,2) - mrSges(1,2) + t57 * r_base(2) + t56 * (t39 * pkin(1) + r_base(2)) + t55 * t25 + t54 * t27) * g(2) + (-t41 * mrSges(2,1) + t39 * mrSges(2,2) - mrSges(1,1) + t57 * r_base(1) + t56 * (t41 * pkin(1) + r_base(1)) + t55 * t27 - t54 * t25) * g(1);
U = t1;
