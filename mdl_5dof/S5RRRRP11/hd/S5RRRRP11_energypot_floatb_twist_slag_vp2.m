% Calculate potential energy for
% S5RRRRP11
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP11_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRP11_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP11_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP11_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP11_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:14:20
% EndTime: 2019-12-31 22:14:21
% DurationCPUTime: 0.50s
% Computational Cost: add. (240->84), mult. (493->99), div. (0->0), fcn. (575->10), ass. (0->43)
t66 = -m(1) - m(2);
t65 = -m(5) - m(6);
t64 = mrSges(3,2) - mrSges(4,3);
t63 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2);
t62 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t61 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t32 = sin(pkin(5));
t36 = sin(qJ(2));
t60 = t32 * t36;
t37 = sin(qJ(1));
t59 = t32 * t37;
t39 = cos(qJ(3));
t58 = t32 * t39;
t40 = cos(qJ(2));
t57 = t32 * t40;
t41 = cos(qJ(1));
t56 = t32 * t41;
t55 = t36 * t37;
t54 = t36 * t41;
t53 = t37 * t40;
t52 = t40 * t41;
t51 = pkin(6) + r_base(3);
t33 = cos(pkin(5));
t50 = t33 * pkin(7) + t51;
t49 = t41 * pkin(1) + pkin(7) * t59 + r_base(1);
t48 = t37 * pkin(1) - pkin(7) * t56 + r_base(2);
t22 = t33 * t53 + t54;
t23 = -t33 * t55 + t52;
t47 = t23 * pkin(2) + pkin(8) * t22 + t49;
t46 = pkin(2) * t60 - pkin(8) * t57 + t50;
t20 = -t33 * t52 + t55;
t21 = t33 * t54 + t53;
t45 = t21 * pkin(2) + t20 * pkin(8) + t48;
t38 = cos(qJ(4));
t35 = sin(qJ(3));
t34 = sin(qJ(4));
t19 = t33 * t35 + t36 * t58;
t18 = -t33 * t39 + t35 * t60;
t12 = t23 * t39 + t35 * t59;
t11 = t23 * t35 - t37 * t58;
t10 = t21 * t39 - t35 * t56;
t9 = t21 * t35 + t39 * t56;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t51 - mrSges(2,3) - m(3) * t50 - t33 * mrSges(3,3) - (t36 * mrSges(3,1) + t40 * mrSges(3,2)) * t32 - m(4) * t46 - t19 * mrSges(4,1) + mrSges(4,3) * t57 + t65 * (t19 * pkin(3) + pkin(9) * t18 + t46) + t62 * (t19 * t38 - t34 * t57) + t61 * (t19 * t34 + t38 * t57) + t63 * t18) * g(3) + (-m(3) * t48 - m(4) * t45 - t37 * mrSges(2,1) - t21 * mrSges(3,1) - t10 * mrSges(4,1) - t41 * mrSges(2,2) + mrSges(3,3) * t56 - mrSges(1,2) + t66 * r_base(2) + t65 * (t10 * pkin(3) + t9 * pkin(9) + t45) + t64 * t20 + t62 * (t10 * t38 + t20 * t34) + t61 * (t10 * t34 - t20 * t38) + t63 * t9) * g(2) + (-m(3) * t49 - m(4) * t47 - t41 * mrSges(2,1) - t23 * mrSges(3,1) - t12 * mrSges(4,1) + t37 * mrSges(2,2) - mrSges(3,3) * t59 - mrSges(1,1) + t66 * r_base(1) + t65 * (t12 * pkin(3) + pkin(9) * t11 + t47) + t62 * (t12 * t38 + t22 * t34) + t61 * (t12 * t34 - t22 * t38) + t64 * t22 + t63 * t11) * g(1);
U = t1;
