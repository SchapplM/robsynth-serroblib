% Calculate joint inertia matrix for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:14:19
% EndTime: 2019-12-05 18:14:20
% DurationCPUTime: 0.38s
% Computational Cost: add. (1399->99), mult. (760->132), div. (0->0), fcn. (614->10), ass. (0->60)
t53 = qJ(1) + pkin(9);
t52 = qJ(3) + t53;
t49 = qJ(4) + t52;
t45 = sin(t49);
t83 = t45 ^ 2;
t46 = cos(t49);
t82 = t46 ^ 2;
t54 = sin(qJ(5));
t86 = Icges(6,5) * t54;
t56 = cos(qJ(5));
t85 = Icges(6,6) * t56;
t35 = t85 + t86;
t84 = t45 * t46;
t38 = t54 * rSges(6,1) + t56 * rSges(6,2);
t79 = m(6) * t38;
t47 = sin(t52);
t78 = pkin(3) * t47;
t48 = cos(t52);
t77 = pkin(3) * t48;
t55 = sin(qJ(1));
t76 = t55 * pkin(1);
t57 = cos(qJ(1));
t75 = t57 * pkin(1);
t74 = rSges(6,1) * t56;
t73 = rSges(6,2) * t54;
t72 = t46 * rSges(6,3) + t45 * t73;
t69 = Icges(6,2) * t56 ^ 2 + Icges(5,3) + (Icges(6,1) * t54 + 0.2e1 * Icges(6,4) * t56) * t54;
t68 = t35 * t82 + (t86 / 0.2e1 + t85 / 0.2e1 + t35 / 0.2e1) * t83;
t67 = -pkin(4) - t74;
t29 = -t48 * rSges(4,1) + t47 * rSges(4,2);
t25 = -t46 * rSges(5,1) + t45 * rSges(5,2);
t66 = Icges(4,3) + t69;
t50 = sin(t53);
t65 = -pkin(2) * t50 - t76;
t51 = cos(t53);
t64 = -pkin(2) * t51 - t75;
t28 = -t47 * rSges(4,1) - t48 * rSges(4,2);
t24 = -t45 * rSges(5,1) - t46 * rSges(5,2);
t58 = Icges(6,5) * t56 - Icges(6,6) * t54;
t23 = t25 - t77;
t22 = t24 - t78;
t10 = t46 * pkin(8) + t67 * t45 + t72;
t33 = t46 * t73;
t11 = t33 + t67 * t46 + (-rSges(6,3) - pkin(8)) * t45;
t8 = t10 - t78;
t9 = t11 - t77;
t40 = -t57 * rSges(2,1) + t55 * rSges(2,2);
t39 = -t55 * rSges(2,1) - t57 * rSges(2,2);
t27 = -t51 * rSges(3,1) + t50 * rSges(3,2) - t75;
t26 = -t50 * rSges(3,1) - t51 * rSges(3,2) - t76;
t21 = t29 + t64;
t20 = t28 + t65;
t15 = Icges(6,3) * t45 + t58 * t46;
t14 = Icges(6,3) * t46 - t58 * t45;
t13 = t23 + t64;
t12 = t22 + t65;
t5 = t64 + t9;
t4 = t65 + t8;
t3 = t46 * (t45 * rSges(6,3) + t46 * t74 - t33) - t45 * (-t45 * t74 + t72);
t1 = [Icges(2,3) + Icges(3,3) + m(2) * (t39 ^ 2 + t40 ^ 2) + m(3) * (t26 ^ 2 + t27 ^ 2) + m(4) * (t20 ^ 2 + t21 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(6) * (t4 ^ 2 + t5 ^ 2) + t66; 0; m(3) + m(4) + m(5) + m(6); m(4) * (t28 * t20 + t29 * t21) + m(5) * (t22 * t12 + t23 * t13) + m(6) * (t8 * t4 + t9 * t5) + t66; 0; m(4) * (t28 ^ 2 + t29 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2) + m(6) * (t8 ^ 2 + t9 ^ 2) + t66; m(5) * (t24 * t12 + t25 * t13) + m(6) * (t10 * t4 + t11 * t5) + t69; 0; m(5) * (t24 * t22 + t25 * t23) + m(6) * (t10 * t8 + t11 * t9) + t69; m(5) * (t24 ^ 2 + t25 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2) + t69; (-t4 * t46 + t45 * t5) * t79 + t68; m(6) * t3; (t45 * t9 - t46 * t8) * t79 + t68; (-t10 * t46 + t11 * t45) * t79 + t68; m(6) * (t3 ^ 2 + (t82 + t83) * t38 ^ 2) + t46 * (t82 * t14 + t15 * t84) + t45 * (t14 * t84 + t83 * t15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
