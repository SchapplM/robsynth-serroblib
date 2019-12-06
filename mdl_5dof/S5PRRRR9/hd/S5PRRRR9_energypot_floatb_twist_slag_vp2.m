% Calculate potential energy for
% S5PRRRR9
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR9_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRR9_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR9_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR9_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR9_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:18:57
% EndTime: 2019-12-05 17:18:57
% DurationCPUTime: 0.63s
% Computational Cost: add. (237->92), mult. (447->107), div. (0->0), fcn. (511->12), ass. (0->42)
t59 = -m(1) - m(2);
t58 = -m(4) - m(5);
t57 = -m(5) * pkin(8) + m(6) * (-pkin(9) - pkin(8)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t24 = qJ(4) + qJ(5);
t19 = sin(t24);
t20 = cos(t24);
t29 = sin(qJ(4));
t32 = cos(qJ(4));
t45 = pkin(4) * t29 + pkin(7);
t56 = -m(6) * t45 - t29 * mrSges(5,1) - t19 * mrSges(6,1) - t32 * mrSges(5,2) - t20 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3);
t18 = pkin(4) * t32 + pkin(3);
t55 = -m(5) * pkin(3) - m(6) * t18 - t32 * mrSges(5,1) - t20 * mrSges(6,1) + t29 * mrSges(5,2) + t19 * mrSges(6,2) - mrSges(4,1);
t25 = sin(pkin(10));
t26 = sin(pkin(5));
t54 = t25 * t26;
t27 = cos(pkin(10));
t53 = t26 * t27;
t30 = sin(qJ(3));
t52 = t26 * t30;
t31 = sin(qJ(2));
t51 = t26 * t31;
t33 = cos(qJ(3));
t50 = t26 * t33;
t34 = cos(qJ(2));
t49 = t26 * t34;
t28 = cos(pkin(5));
t48 = t28 * t31;
t47 = t28 * t34;
t46 = qJ(1) + r_base(3);
t44 = t27 * pkin(1) + pkin(6) * t54 + r_base(1);
t43 = t28 * pkin(6) + t46;
t10 = -t25 * t48 + t27 * t34;
t42 = t10 * pkin(2) + t44;
t41 = pkin(2) * t51 + t43;
t40 = t25 * pkin(1) - pkin(6) * t53 + r_base(2);
t8 = t25 * t34 + t27 * t48;
t38 = t8 * pkin(2) + t40;
t37 = -pkin(7) * t49 + t41;
t12 = t28 * t30 + t31 * t50;
t9 = t25 * t47 + t27 * t31;
t7 = t25 * t31 - t27 * t47;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t46 - mrSges(2,3) - m(3) * t43 - t28 * mrSges(3,3) - (t31 * mrSges(3,1) + t34 * mrSges(3,2)) * t26 - m(4) * t37 - t12 * mrSges(4,1) + mrSges(4,3) * t49 - m(5) * (pkin(3) * t12 + t37) - (t12 * t32 - t29 * t49) * mrSges(5,1) - (-t12 * t29 - t32 * t49) * mrSges(5,2) - m(6) * (t12 * t18 - t45 * t49 + t41) - (t12 * t20 - t19 * t49) * mrSges(6,1) - (-t12 * t19 - t20 * t49) * mrSges(6,2) + t57 * (-t28 * t33 + t30 * t51)) * g(3) + (-m(3) * t40 - m(6) * t38 - t25 * mrSges(2,1) - t8 * mrSges(3,1) - t27 * mrSges(2,2) + mrSges(3,3) * t53 - mrSges(1,2) + t59 * r_base(2) + t58 * (t7 * pkin(7) + t38) + t55 * (-t27 * t52 + t33 * t8) + t56 * t7 + t57 * (t27 * t50 + t30 * t8)) * g(2) + (-m(3) * t44 - m(6) * t42 - t27 * mrSges(2,1) - t10 * mrSges(3,1) + t25 * mrSges(2,2) - mrSges(3,3) * t54 - mrSges(1,1) + t59 * r_base(1) + t58 * (t9 * pkin(7) + t42) + t55 * (t10 * t33 + t25 * t52) + t56 * t9 + t57 * (t10 * t30 - t25 * t50)) * g(1);
U = t1;
