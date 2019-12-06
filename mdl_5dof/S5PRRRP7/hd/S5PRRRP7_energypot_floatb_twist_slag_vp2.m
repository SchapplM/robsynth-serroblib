% Calculate potential energy for
% S5PRRRP7
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRP7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRP7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP7_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP7_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:54:01
% EndTime: 2019-12-05 16:54:02
% DurationCPUTime: 0.55s
% Computational Cost: add. (225->92), mult. (447->109), div. (0->0), fcn. (511->10), ass. (0->43)
t62 = -m(1) - m(2);
t61 = -mrSges(5,1) - mrSges(6,1);
t60 = -mrSges(5,2) - mrSges(6,2);
t33 = sin(qJ(4));
t48 = pkin(4) * t33 + pkin(7);
t59 = -m(6) * t48 + mrSges(3,2) - mrSges(4,3);
t58 = -m(5) * pkin(8) + m(6) * (-qJ(5) - pkin(8)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t28 = sin(pkin(9));
t29 = sin(pkin(5));
t57 = t28 * t29;
t30 = cos(pkin(9));
t56 = t29 * t30;
t34 = sin(qJ(3));
t55 = t29 * t34;
t35 = sin(qJ(2));
t54 = t29 * t35;
t37 = cos(qJ(3));
t53 = t29 * t37;
t38 = cos(qJ(2));
t52 = t29 * t38;
t31 = cos(pkin(5));
t51 = t31 * t35;
t50 = t31 * t38;
t49 = qJ(1) + r_base(3);
t47 = t30 * pkin(1) + pkin(6) * t57 + r_base(1);
t46 = t31 * pkin(6) + t49;
t16 = -t28 * t51 + t30 * t38;
t45 = t16 * pkin(2) + t47;
t44 = pkin(2) * t54 + t46;
t43 = t28 * pkin(1) - pkin(6) * t56 + r_base(2);
t14 = t28 * t38 + t30 * t51;
t42 = t14 * pkin(2) + t43;
t15 = t28 * t50 + t30 * t35;
t41 = pkin(7) * t15 + t45;
t40 = -pkin(7) * t52 + t44;
t13 = t28 * t35 - t30 * t50;
t39 = pkin(7) * t13 + t42;
t36 = cos(qJ(4));
t24 = pkin(4) * t36 + pkin(3);
t18 = t31 * t34 + t35 * t53;
t8 = t16 * t37 + t28 * t55;
t6 = t14 * t37 - t30 * t55;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t49 - mrSges(2,3) - m(3) * t46 - t31 * mrSges(3,3) - (t35 * mrSges(3,1) + t38 * mrSges(3,2)) * t29 - m(4) * t40 - t18 * mrSges(4,1) + mrSges(4,3) * t52 - m(5) * (t18 * pkin(3) + t40) - m(6) * (t18 * t24 - t48 * t52 + t44) + t60 * (-t18 * t33 - t36 * t52) + t61 * (t18 * t36 - t33 * t52) + t58 * (-t31 * t37 + t34 * t54)) * g(3) + (-mrSges(1,2) - t28 * mrSges(2,1) - t30 * mrSges(2,2) - m(3) * t43 - t14 * mrSges(3,1) + mrSges(3,3) * t56 - m(4) * t39 - t6 * mrSges(4,1) - m(5) * (pkin(3) * t6 + t39) - m(6) * (t24 * t6 + t42) + t62 * r_base(2) + t61 * (t13 * t33 + t36 * t6) + t59 * t13 + t60 * (t13 * t36 - t33 * t6) + t58 * (t14 * t34 + t30 * t53)) * g(2) + (-mrSges(1,1) - t30 * mrSges(2,1) + t28 * mrSges(2,2) - m(3) * t47 - t16 * mrSges(3,1) - mrSges(3,3) * t57 - m(4) * t41 - t8 * mrSges(4,1) - m(5) * (pkin(3) * t8 + t41) - m(6) * (t8 * t24 + t45) + t62 * r_base(1) + t61 * (t15 * t33 + t36 * t8) + t60 * (t15 * t36 - t33 * t8) + t59 * t15 + t58 * (t16 * t34 - t28 * t53)) * g(1);
U = t1;
