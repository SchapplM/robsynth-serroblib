% Calculate potential energy for
% S5PRRPR5
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% m_mdh [6x1]
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
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRPR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR5_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR5_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:25:12
% EndTime: 2019-12-05 16:25:13
% DurationCPUTime: 0.70s
% Computational Cost: add. (245->98), mult. (383->114), div. (0->0), fcn. (422->12), ass. (0->43)
t64 = -m(1) - m(2);
t63 = -m(5) - m(6);
t62 = m(4) * pkin(7) - mrSges(3,2) + mrSges(4,3);
t61 = -m(6) * pkin(8) + mrSges(5,2) - mrSges(6,3);
t35 = sin(qJ(5));
t38 = cos(qJ(5));
t60 = -m(6) * pkin(4) - t38 * mrSges(6,1) + t35 * mrSges(6,2) - mrSges(5,1);
t59 = -t35 * mrSges(6,1) - t38 * mrSges(6,2) - mrSges(5,3) - t62;
t36 = sin(qJ(3));
t58 = pkin(3) * t36;
t30 = sin(pkin(9));
t31 = sin(pkin(5));
t57 = t30 * t31;
t32 = cos(pkin(9));
t56 = t31 * t32;
t55 = t31 * t36;
t37 = sin(qJ(2));
t54 = t31 * t37;
t39 = cos(qJ(3));
t53 = t31 * t39;
t40 = cos(qJ(2));
t52 = t31 * t40;
t33 = cos(pkin(5));
t51 = t33 * t37;
t50 = t33 * t40;
t49 = t30 * pkin(1) + r_base(2);
t48 = t30 * t55;
t47 = qJ(1) + r_base(3);
t46 = t32 * pkin(1) + pkin(6) * t57 + r_base(1);
t45 = t33 * pkin(6) + t47;
t44 = -pkin(6) * t56 + t49;
t23 = pkin(3) * t39 + pkin(2);
t34 = -qJ(4) - pkin(7);
t42 = t23 * t54 + t33 * t58 + t34 * t52 + t45;
t29 = qJ(3) + pkin(10);
t25 = cos(t29);
t24 = sin(t29);
t14 = -t30 * t51 + t32 * t40;
t13 = t30 * t50 + t32 * t37;
t12 = t30 * t40 + t32 * t51;
t11 = t30 * t37 - t32 * t50;
t8 = t24 * t33 + t25 * t54;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t47 - mrSges(2,3) - m(5) * t42 - t8 * mrSges(5,1) + mrSges(5,3) * t52 - m(6) * (pkin(4) * t8 + t42) - (-t35 * t52 + t8 * t38) * mrSges(6,1) - (-t8 * t35 - t38 * t52) * mrSges(6,2) + t61 * (t24 * t54 - t33 * t25) + (-m(3) - m(4)) * t45 + (-t36 * mrSges(4,1) - t39 * mrSges(4,2) - mrSges(3,3)) * t33 + (t62 * t40 + (-m(4) * pkin(2) - t39 * mrSges(4,1) + t36 * mrSges(4,2) - mrSges(3,1)) * t37) * t31) * g(3) + (-mrSges(1,2) - t30 * mrSges(2,1) - t32 * mrSges(2,2) - m(3) * t44 - t12 * mrSges(3,1) + mrSges(3,3) * t56 - m(4) * (pkin(2) * t12 + t44) - (t12 * t39 - t32 * t55) * mrSges(4,1) - (-t12 * t36 - t32 * t53) * mrSges(4,2) + t64 * r_base(2) + t63 * (t12 * t23 - t11 * t34 + (-pkin(6) - t58) * t56 + t49) + t61 * (t12 * t24 + t25 * t56) + t60 * (t12 * t25 - t24 * t56) + t59 * t11) * g(2) + (-mrSges(1,1) - t32 * mrSges(2,1) + t30 * mrSges(2,2) - m(3) * t46 - t14 * mrSges(3,1) - mrSges(3,3) * t57 - m(4) * (pkin(2) * t14 + t46) - (t14 * t39 + t48) * mrSges(4,1) - (-t14 * t36 + t30 * t53) * mrSges(4,2) + t64 * r_base(1) + t63 * (pkin(3) * t48 - t13 * t34 + t14 * t23 + t46) + t61 * (t14 * t24 - t25 * t57) + t60 * (t14 * t25 + t24 * t57) + t59 * t13) * g(1);
U = t1;
