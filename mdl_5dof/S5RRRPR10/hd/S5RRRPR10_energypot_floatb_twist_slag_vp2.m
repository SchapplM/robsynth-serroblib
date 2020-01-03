% Calculate potential energy for
% S5RRRPR10
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR10_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRPR10_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR10_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR10_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR10_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:26:32
% EndTime: 2019-12-31 21:26:32
% DurationCPUTime: 0.70s
% Computational Cost: add. (245->98), mult. (383->110), div. (0->0), fcn. (422->12), ass. (0->43)
t64 = -m(1) - m(2);
t63 = -m(5) - m(6);
t62 = m(4) * pkin(8) - mrSges(3,2) + mrSges(4,3);
t61 = -m(6) * pkin(9) + mrSges(5,2) - mrSges(6,3);
t33 = sin(qJ(5));
t37 = cos(qJ(5));
t60 = -m(6) * pkin(4) - t37 * mrSges(6,1) + t33 * mrSges(6,2) - mrSges(5,1);
t59 = -t33 * mrSges(6,1) - t37 * mrSges(6,2) - mrSges(5,3) - t62;
t34 = sin(qJ(3));
t58 = pkin(3) * t34;
t30 = sin(pkin(5));
t35 = sin(qJ(2));
t57 = t30 * t35;
t36 = sin(qJ(1));
t56 = t30 * t36;
t39 = cos(qJ(2));
t55 = t30 * t39;
t40 = cos(qJ(1));
t54 = t30 * t40;
t53 = t35 * t36;
t52 = t35 * t40;
t51 = t36 * t39;
t50 = t39 * t40;
t49 = pkin(6) + r_base(3);
t48 = t36 * pkin(1) + r_base(2);
t47 = t34 * t56;
t31 = cos(pkin(5));
t46 = t31 * pkin(7) + t49;
t45 = t40 * pkin(1) + pkin(7) * t56 + r_base(1);
t44 = -pkin(7) * t54 + t48;
t38 = cos(qJ(3));
t23 = pkin(3) * t38 + pkin(2);
t32 = -qJ(4) - pkin(8);
t43 = t23 * t57 + t31 * t58 + t32 * t55 + t46;
t29 = qJ(3) + pkin(10);
t25 = cos(t29);
t24 = sin(t29);
t14 = -t31 * t53 + t50;
t13 = t31 * t51 + t52;
t12 = t31 * t52 + t51;
t11 = -t31 * t50 + t53;
t8 = t24 * t31 + t25 * t57;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t49 - mrSges(2,3) - m(5) * t43 - t8 * mrSges(5,1) + mrSges(5,3) * t55 - m(6) * (pkin(4) * t8 + t43) - (-t33 * t55 + t37 * t8) * mrSges(6,1) - (-t33 * t8 - t37 * t55) * mrSges(6,2) + t61 * (t24 * t57 - t31 * t25) + (-m(3) - m(4)) * t46 + (-t34 * mrSges(4,1) - t38 * mrSges(4,2) - mrSges(3,3)) * t31 + (t62 * t39 + (-m(4) * pkin(2) - t38 * mrSges(4,1) + t34 * mrSges(4,2) - mrSges(3,1)) * t35) * t30) * g(3) + (-mrSges(1,2) - t36 * mrSges(2,1) - t40 * mrSges(2,2) - m(3) * t44 - t12 * mrSges(3,1) + mrSges(3,3) * t54 - m(4) * (t12 * pkin(2) + t44) - (t12 * t38 - t34 * t54) * mrSges(4,1) - (-t12 * t34 - t38 * t54) * mrSges(4,2) + t64 * r_base(2) + t63 * (t12 * t23 - t11 * t32 + (-pkin(7) - t58) * t54 + t48) + t61 * (t12 * t24 + t25 * t54) + t60 * (t12 * t25 - t24 * t54) + t59 * t11) * g(2) + (-mrSges(1,1) - t40 * mrSges(2,1) + t36 * mrSges(2,2) - m(3) * t45 - t14 * mrSges(3,1) - mrSges(3,3) * t56 - m(4) * (pkin(2) * t14 + t45) - (t14 * t38 + t47) * mrSges(4,1) - (-t14 * t34 + t38 * t56) * mrSges(4,2) + t64 * r_base(1) + t63 * (pkin(3) * t47 - t13 * t32 + t14 * t23 + t45) + t61 * (t14 * t24 - t25 * t56) + t60 * (t14 * t25 + t24 * t56) + t59 * t13) * g(1);
U = t1;
