% Calculate potential energy for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPPRR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PPPRR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR2_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:22
% EndTime: 2019-12-05 14:59:22
% DurationCPUTime: 0.65s
% Computational Cost: add. (176->82), mult. (325->93), div. (0->0), fcn. (348->10), ass. (0->37)
t60 = mrSges(3,2) - mrSges(4,3);
t59 = -m(1) - m(2);
t58 = -m(5) - m(6);
t57 = mrSges(2,2) - mrSges(3,3);
t28 = sin(pkin(8));
t31 = cos(pkin(8));
t56 = -t31 * mrSges(3,1) + t28 * t60 - mrSges(2,1);
t55 = -m(6) * pkin(6) + mrSges(5,2) - mrSges(6,3);
t33 = sin(qJ(5));
t35 = cos(qJ(5));
t54 = -t33 * mrSges(6,1) - t35 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t53 = -m(6) * pkin(4) - mrSges(6,1) * t35 + mrSges(6,2) * t33 - mrSges(5,1);
t27 = sin(pkin(9));
t51 = t27 * t28;
t34 = sin(qJ(4));
t50 = t28 * t34;
t36 = cos(qJ(4));
t49 = t28 * t36;
t29 = sin(pkin(7));
t48 = t29 * t31;
t32 = cos(pkin(7));
t47 = t31 * t32;
t46 = qJ(3) * t28;
t26 = qJ(1) + r_base(3);
t45 = t32 * pkin(1) + t29 * qJ(2) + r_base(1);
t43 = t29 * pkin(1) - qJ(2) * t32 + r_base(2);
t42 = pkin(2) * t47 + t32 * t46 + t45;
t41 = t28 * pkin(2) - qJ(3) * t31 + t26;
t40 = pkin(2) * t48 + t29 * t46 + t43;
t30 = cos(pkin(9));
t39 = t28 * t30 * pkin(3) + pkin(5) * t51 + t41;
t12 = t30 * t49 - t31 * t34;
t10 = t27 * t29 + t30 * t47;
t9 = t27 * t47 - t29 * t30;
t8 = -t27 * t32 + t30 * t48;
t7 = t27 * t48 + t30 * t32;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - m(4) * t41 - m(5) * t39 - t12 * mrSges(5,1) - mrSges(5,3) * t51 - m(6) * (pkin(4) * t12 + t39) - (t12 * t35 + t33 * t51) * mrSges(6,1) - (-t12 * t33 + t35 * t51) * mrSges(6,2) - t60 * t31 + (-t30 * mrSges(4,1) + mrSges(4,2) * t27 - mrSges(3,1)) * t28 + (-m(2) - m(3)) * t26 + t55 * (t30 * t50 + t31 * t36)) * g(3) + (-m(3) * t43 - m(4) * t40 - t8 * mrSges(4,1) - mrSges(1,2) + t59 * r_base(2) + t58 * (t8 * pkin(3) + pkin(5) * t7 + t40) + t53 * (t29 * t50 + t36 * t8) + t54 * t7 - t57 * t32 + t55 * (-t29 * t49 + t34 * t8) + t56 * t29) * g(2) + (-m(3) * t45 - m(4) * t42 - t10 * mrSges(4,1) - mrSges(1,1) + t59 * r_base(1) + t58 * (t10 * pkin(3) + pkin(5) * t9 + t42) + t53 * (t10 * t36 + t32 * t50) + t54 * t9 + t55 * (t10 * t34 - t32 * t49) + t57 * t29 + t56 * t32) * g(1);
U = t1;
