% Calculate potential energy for
% S5PRRPR7
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
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPR7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRPR7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR7_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR7_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:34:53
% EndTime: 2019-12-05 16:34:54
% DurationCPUTime: 0.56s
% Computational Cost: add. (261->87), mult. (548->106), div. (0->0), fcn. (651->12), ass. (0->44)
t67 = -m(1) - m(2);
t66 = -m(5) - m(6);
t65 = mrSges(3,2) - mrSges(4,3);
t64 = -m(6) * pkin(8) + mrSges(5,2) - mrSges(6,3);
t38 = sin(qJ(5));
t41 = cos(qJ(5));
t63 = -t38 * mrSges(6,1) - t41 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t62 = -m(6) * pkin(4) - t41 * mrSges(6,1) + t38 * mrSges(6,2) - mrSges(5,1);
t33 = sin(pkin(9));
t34 = sin(pkin(5));
t61 = t33 * t34;
t36 = cos(pkin(9));
t60 = t34 * t36;
t39 = sin(qJ(3));
t59 = t34 * t39;
t40 = sin(qJ(2));
t58 = t34 * t40;
t42 = cos(qJ(3));
t57 = t34 * t42;
t43 = cos(qJ(2));
t56 = t34 * t43;
t37 = cos(pkin(5));
t55 = t37 * t40;
t54 = t37 * t43;
t53 = qJ(1) + r_base(3);
t52 = t36 * pkin(1) + pkin(6) * t61 + r_base(1);
t51 = t37 * pkin(6) + t53;
t50 = t33 * pkin(1) - pkin(6) * t60 + r_base(2);
t20 = t33 * t54 + t36 * t40;
t21 = -t33 * t55 + t36 * t43;
t49 = t21 * pkin(2) + pkin(7) * t20 + t52;
t48 = pkin(2) * t58 - pkin(7) * t56 + t51;
t18 = t33 * t40 - t36 * t54;
t19 = t33 * t43 + t36 * t55;
t47 = t19 * pkin(2) + pkin(7) * t18 + t50;
t35 = cos(pkin(10));
t32 = sin(pkin(10));
t23 = t37 * t39 + t40 * t57;
t22 = -t37 * t42 + t39 * t58;
t12 = t21 * t42 + t33 * t59;
t11 = t21 * t39 - t33 * t57;
t10 = t19 * t42 - t36 * t59;
t9 = t19 * t39 + t36 * t57;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t53 - mrSges(2,3) - m(3) * t51 - t37 * mrSges(3,3) - (t40 * mrSges(3,1) + t43 * mrSges(3,2)) * t34 - m(4) * t48 - t23 * mrSges(4,1) + mrSges(4,3) * t56 + t66 * (t23 * pkin(3) + t22 * qJ(4) + t48) + t64 * (t23 * t32 + t35 * t56) + t62 * (t23 * t35 - t32 * t56) + t63 * t22) * g(3) + (-m(3) * t50 - m(4) * t47 - t33 * mrSges(2,1) - t19 * mrSges(3,1) - t10 * mrSges(4,1) - t36 * mrSges(2,2) + mrSges(3,3) * t60 - mrSges(1,2) + t67 * r_base(2) + t66 * (t10 * pkin(3) + qJ(4) * t9 + t47) + t62 * (t10 * t35 + t18 * t32) + t63 * t9 + t65 * t18 + t64 * (t10 * t32 - t18 * t35)) * g(2) + (-m(3) * t52 - m(4) * t49 - t36 * mrSges(2,1) - t21 * mrSges(3,1) - t12 * mrSges(4,1) + t33 * mrSges(2,2) - mrSges(3,3) * t61 - mrSges(1,1) + t67 * r_base(1) + t66 * (t12 * pkin(3) + qJ(4) * t11 + t49) + t64 * (t12 * t32 - t20 * t35) + t65 * t20 + t62 * (t12 * t35 + t20 * t32) + t63 * t11) * g(1);
U = t1;
