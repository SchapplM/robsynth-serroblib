% Calculate potential energy for
% S5PRPRR6
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRPRR6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR6_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:21
% EndTime: 2019-12-05 15:56:22
% DurationCPUTime: 0.71s
% Computational Cost: add. (245->98), mult. (383->112), div. (0->0), fcn. (422->12), ass. (0->41)
t62 = -m(1) - m(2);
t61 = -m(5) - m(6);
t60 = -m(6) * pkin(8) + mrSges(5,2) - mrSges(6,3);
t59 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3);
t37 = sin(qJ(5));
t39 = cos(qJ(5));
t58 = -m(6) * pkin(4) - t39 * mrSges(6,1) + t37 * mrSges(6,2) - mrSges(5,1);
t57 = -t37 * mrSges(6,1) - t39 * mrSges(6,2) - mrSges(5,3) - t59;
t30 = sin(pkin(10));
t56 = pkin(3) * t30;
t31 = sin(pkin(9));
t32 = sin(pkin(5));
t55 = t31 * t32;
t34 = cos(pkin(9));
t54 = t32 * t34;
t38 = sin(qJ(2));
t53 = t32 * t38;
t40 = cos(qJ(2));
t52 = t32 * t40;
t35 = cos(pkin(5));
t51 = t35 * t38;
t50 = t35 * t40;
t49 = t31 * pkin(1) + r_base(2);
t48 = t30 * t55;
t47 = qJ(1) + r_base(3);
t46 = t34 * pkin(1) + pkin(6) * t55 + r_base(1);
t45 = t35 * pkin(6) + t47;
t44 = -pkin(6) * t54 + t49;
t33 = cos(pkin(10));
t23 = pkin(3) * t33 + pkin(2);
t36 = -pkin(7) - qJ(3);
t42 = t23 * t53 + t35 * t56 + t36 * t52 + t45;
t29 = pkin(10) + qJ(4);
t25 = cos(t29);
t24 = sin(t29);
t14 = -t31 * t51 + t34 * t40;
t13 = t31 * t50 + t34 * t38;
t12 = t31 * t40 + t34 * t51;
t11 = t31 * t38 - t34 * t50;
t8 = t24 * t35 + t25 * t53;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t47 - mrSges(2,3) - m(5) * t42 - t8 * mrSges(5,1) + mrSges(5,3) * t52 - m(6) * (pkin(4) * t8 + t42) - (-t37 * t52 + t8 * t39) * mrSges(6,1) - (-t8 * t37 - t39 * t52) * mrSges(6,2) + t60 * (t24 * t53 - t35 * t25) + (-m(3) - m(4)) * t45 + (-t30 * mrSges(4,1) - t33 * mrSges(4,2) - mrSges(3,3)) * t35 + (t59 * t40 + (-m(4) * pkin(2) - t33 * mrSges(4,1) + t30 * mrSges(4,2) - mrSges(3,1)) * t38) * t32) * g(3) + (-mrSges(1,2) - t31 * mrSges(2,1) - t34 * mrSges(2,2) - m(3) * t44 - t12 * mrSges(3,1) + mrSges(3,3) * t54 - m(4) * (pkin(2) * t12 + t44) - (t12 * t33 - t30 * t54) * mrSges(4,1) - (-t12 * t30 - t33 * t54) * mrSges(4,2) + t62 * r_base(2) + t61 * (t12 * t23 - t11 * t36 + (-pkin(6) - t56) * t54 + t49) + t60 * (t12 * t24 + t25 * t54) + t58 * (t12 * t25 - t24 * t54) + t57 * t11) * g(2) + (-mrSges(1,1) - t34 * mrSges(2,1) + t31 * mrSges(2,2) - m(3) * t46 - t14 * mrSges(3,1) - mrSges(3,3) * t55 - m(4) * (pkin(2) * t14 + t46) - (t14 * t33 + t48) * mrSges(4,1) - (-t14 * t30 + t33 * t55) * mrSges(4,2) + t62 * r_base(1) + t61 * (pkin(3) * t48 - t13 * t36 + t14 * t23 + t46) + t60 * (t14 * t24 - t25 * t55) + t58 * (t14 * t25 + t24 * t55) + t57 * t13) * g(1);
U = t1;
