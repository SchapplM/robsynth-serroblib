% Calculate potential energy for
% S6RRRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP9_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRP9_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP9_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_energypot_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP9_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP9_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:21:46
% EndTime: 2019-03-09 17:21:46
% DurationCPUTime: 0.42s
% Computational Cost: add. (224->111), mult. (392->127), div. (0->0), fcn. (430->8), ass. (0->42)
t65 = rSges(7,1) + pkin(5);
t64 = rSges(7,3) + qJ(6);
t33 = sin(qJ(2));
t34 = sin(qJ(1));
t38 = cos(qJ(1));
t63 = (g(1) * t38 + g(2) * t34) * t33;
t32 = sin(qJ(3));
t58 = t32 * t33;
t57 = t33 * t34;
t36 = cos(qJ(3));
t56 = t33 * t36;
t55 = t33 * t38;
t37 = cos(qJ(2));
t54 = t34 * t37;
t53 = t38 * t32;
t52 = t38 * t36;
t51 = rSges(5,3) + qJ(4);
t50 = pkin(6) + r_base(3);
t49 = t34 * pkin(1) + r_base(2);
t48 = t33 * pkin(2) + t50;
t47 = t38 * pkin(1) + t34 * pkin(7) + r_base(1);
t46 = pkin(3) * t56 + qJ(4) * t58 + t48;
t45 = t38 * t37 * pkin(2) + pkin(8) * t55 + t47;
t16 = t34 * t32 + t37 * t52;
t44 = t16 * pkin(3) + t45;
t43 = pkin(2) * t54 - t38 * pkin(7) + pkin(8) * t57 + t49;
t42 = pkin(4) * t56 + t37 * pkin(9) + t46;
t14 = t36 * t54 - t53;
t41 = t14 * pkin(3) + t43;
t15 = -t34 * t36 + t37 * t53;
t40 = t16 * pkin(4) + t15 * qJ(4) + t44;
t13 = t32 * t54 + t52;
t39 = t14 * pkin(4) + t13 * qJ(4) + t41;
t35 = cos(qJ(5));
t31 = sin(qJ(5));
t8 = (t31 * t32 + t35 * t36) * t33;
t7 = t31 * t56 - t35 * t58;
t4 = t15 * t31 + t16 * t35;
t3 = -t15 * t35 + t16 * t31;
t2 = t13 * t31 + t14 * t35;
t1 = -t13 * t35 + t14 * t31;
t5 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t38 * rSges(2,1) - t34 * rSges(2,2) + r_base(1)) + g(2) * (t34 * rSges(2,1) + t38 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t50)) - m(3) * (g(1) * (t34 * rSges(3,3) + t47) + g(2) * (rSges(3,1) * t54 - rSges(3,2) * t57 + t49) + g(3) * (t33 * rSges(3,1) + t37 * rSges(3,2) + t50) + (g(1) * (rSges(3,1) * t37 - rSges(3,2) * t33) + g(2) * (-rSges(3,3) - pkin(7))) * t38) - m(4) * (g(1) * (t16 * rSges(4,1) - t15 * rSges(4,2) + rSges(4,3) * t55 + t45) + g(2) * (t14 * rSges(4,1) - t13 * rSges(4,2) + rSges(4,3) * t57 + t43) + g(3) * ((-rSges(4,3) - pkin(8)) * t37 + (rSges(4,1) * t36 - rSges(4,2) * t32) * t33 + t48)) - m(5) * (g(1) * (t16 * rSges(5,1) + rSges(5,2) * t55 + t51 * t15 + t44) + g(2) * (t14 * rSges(5,1) + rSges(5,2) * t57 + t51 * t13 + t41) + g(3) * ((-rSges(5,2) - pkin(8)) * t37 + (rSges(5,1) * t36 + rSges(5,3) * t32) * t33 + t46)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t40) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t39) + g(3) * (t8 * rSges(6,1) - t7 * rSges(6,2) + (rSges(6,3) - pkin(8)) * t37 + t42) + (-rSges(6,3) - pkin(9)) * t63) - m(7) * (g(1) * (t64 * t3 + t65 * t4 + t40) + g(2) * (t64 * t1 + t65 * t2 + t39) + g(3) * (t65 * t8 + t64 * t7 + (rSges(7,2) - pkin(8)) * t37 + t42) + (-rSges(7,2) - pkin(9)) * t63);
U  = t5;
