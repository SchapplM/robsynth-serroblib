% Calculate potential energy for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPP3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRPP3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP3_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP3_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:11:54
% EndTime: 2019-12-05 16:11:55
% DurationCPUTime: 0.51s
% Computational Cost: add. (163->74), mult. (291->81), div. (0->0), fcn. (301->8), ass. (0->35)
t60 = mrSges(3,2) - mrSges(4,3);
t59 = -m(1) - m(2);
t58 = -m(5) - m(6);
t57 = mrSges(2,2) - mrSges(3,3);
t56 = -mrSges(5,3) - mrSges(6,2);
t32 = sin(qJ(2));
t34 = cos(qJ(2));
t55 = -t34 * mrSges(3,1) + t32 * t60 - mrSges(2,1);
t54 = mrSges(4,2) + t56;
t53 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t52 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t51 = pkin(2) * t34;
t28 = sin(pkin(7));
t50 = t28 * t32;
t30 = cos(pkin(7));
t49 = t30 * t32;
t31 = sin(qJ(3));
t48 = t31 * t34;
t46 = t32 * t31;
t33 = cos(qJ(3));
t45 = t32 * t33;
t44 = t33 * t34;
t26 = qJ(1) + r_base(3);
t43 = t30 * pkin(1) + t28 * pkin(5) + r_base(1);
t41 = t28 * pkin(1) - t30 * pkin(5) + r_base(2);
t40 = pkin(6) * t49 + t30 * t51 + t43;
t39 = t32 * pkin(2) - t34 * pkin(6) + t26;
t38 = pkin(6) * t50 + t28 * t51 + t41;
t29 = cos(pkin(8));
t27 = sin(pkin(8));
t12 = t28 * t31 + t30 * t44;
t11 = -t28 * t33 + t30 * t48;
t8 = t28 * t44 - t30 * t31;
t7 = t28 * t48 + t30 * t33;
t1 = (-m(1) * r_base(3) - m(4) * t39 - mrSges(1,3) - mrSges(2,3) + t58 * (pkin(3) * t45 + qJ(4) * t46 + t39) + t52 * (t27 * t45 + t34 * t29) + t56 * t46 - t60 * t34 + (-t33 * mrSges(4,1) + t31 * mrSges(4,2) - mrSges(3,1)) * t32 + (-m(2) - m(3)) * t26 + t53 * (-t34 * t27 + t29 * t45)) * g(3) + (-m(3) * t41 - m(4) * t38 - t8 * mrSges(4,1) - mrSges(1,2) + t59 * r_base(2) + t58 * (t8 * pkin(3) + t7 * qJ(4) + t38) - t57 * t30 + t53 * (t27 * t50 + t8 * t29) + t52 * (t8 * t27 - t29 * t50) + t54 * t7 + t55 * t28) * g(2) + (-m(3) * t43 - m(4) * t40 - t12 * mrSges(4,1) - mrSges(1,1) + t59 * r_base(1) + t58 * (t12 * pkin(3) + t11 * qJ(4) + t40) + t53 * (t12 * t29 + t27 * t49) + t52 * (t12 * t27 - t29 * t49) + t57 * t28 + t55 * t30 + t54 * t11) * g(1);
U = t1;
