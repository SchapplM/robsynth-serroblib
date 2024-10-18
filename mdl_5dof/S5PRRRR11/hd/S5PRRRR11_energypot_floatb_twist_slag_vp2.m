% Calculate potential energy for
% S5PRRRR11
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
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR11_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRR11_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR11_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR11_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:07
% EndTime: 2024-09-27 21:45:07
% DurationCPUTime: 0.09s
% Computational Cost: add. (260->75), mult. (190->66), div. (0->0), fcn. (154->22), ass. (0->34)
t40 = sin(qJ(3));
t41 = cos(qJ(3));
t63 = (-m(5) * pkin(3) - mrSges(4,1)) * t40 - m(6) * (pkin(4) * sin(qJ(4)) * t41 + (cos(qJ(4)) * pkin(4) + pkin(3)) * t40) - t41 * mrSges(4,2);
t42 = pkin(7) + pkin(8);
t60 = m(5) * t42 + m(6) * (pkin(9) + t42) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3) + m(4) * pkin(7);
t57 = -m(1) - m(2);
t56 = -m(3) - m(4) - m(5) - m(6);
t21 = t41 * pkin(3) + pkin(2);
t34 = qJ(3) + qJ(4);
t27 = cos(t34);
t30 = qJ(5) + t34;
t55 = -m(4) * pkin(2) - m(5) * t21 - m(6) * (pkin(4) * t27 + t21) - t41 * mrSges(4,1) - t27 * mrSges(5,1) - cos(t30) * mrSges(6,1) + t40 * mrSges(4,2) + sin(t34) * mrSges(5,2) + sin(t30) * mrSges(6,2) - mrSges(3,1);
t24 = pkin(5) + t34;
t10 = sin(t24) / 0.2e1;
t25 = pkin(5) - t34;
t11 = cos(t25) / 0.2e1;
t17 = -qJ(5) + t25;
t12 = sin(t17);
t16 = qJ(5) + t24;
t13 = cos(t16);
t14 = sin(t25);
t15 = cos(t24);
t36 = sin(pkin(5));
t38 = cos(pkin(5));
t8 = sin(t16) / 0.2e1;
t9 = cos(t17) / 0.2e1;
t54 = -(t10 - t14 / 0.2e1) * mrSges(5,1) - (t8 - t12 / 0.2e1) * mrSges(6,1) - (t11 + t15 / 0.2e1) * mrSges(5,2) - (t9 + t13 / 0.2e1) * mrSges(6,2) - mrSges(3,2) + t60 * t36 + t63 * t38;
t45 = qJ(1) + r_base(3);
t37 = cos(pkin(10));
t35 = sin(pkin(10));
t32 = pkin(10) + qJ(2);
t23 = cos(t32);
t22 = sin(t32);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t45 - mrSges(2,3) - mrSges(3,3) - (t11 - t15 / 0.2e1) * mrSges(5,1) - (t10 + t14 / 0.2e1) * mrSges(5,2) - (t9 - t13 / 0.2e1) * mrSges(6,1) - (t8 + t12 / 0.2e1) * mrSges(6,2) + t63 * t36 + t56 * (pkin(6) + t45) - t60 * t38) * g(3) + (-t35 * mrSges(2,1) - t37 * mrSges(2,2) - mrSges(1,2) + t57 * r_base(2) + t56 * (t35 * pkin(1) + r_base(2)) + t55 * t22 + t54 * t23) * g(2) + (-t37 * mrSges(2,1) + t35 * mrSges(2,2) - mrSges(1,1) + t57 * r_base(1) + t56 * (t37 * pkin(1) + r_base(1)) + t55 * t23 - t54 * t22) * g(1);
U = t1;
