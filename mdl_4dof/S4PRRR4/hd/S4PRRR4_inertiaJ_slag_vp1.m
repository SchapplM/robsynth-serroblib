% Calculate joint inertia matrix for
% S4PRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR4_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR4_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR4_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:30
% EndTime: 2019-12-31 16:32:31
% DurationCPUTime: 0.43s
% Computational Cost: add. (1069->111), mult. (876->167), div. (0->0), fcn. (778->6), ass. (0->61)
t55 = pkin(7) + qJ(2);
t51 = sin(t55);
t52 = cos(t55);
t91 = t51 * t52;
t56 = qJ(3) + qJ(4);
t53 = sin(t56);
t54 = cos(t56);
t72 = rSges(5,1) * t54 - rSges(5,2) * t53;
t49 = t51 ^ 2;
t50 = t52 ^ 2;
t90 = t51 / 0.2e1;
t89 = -t52 / 0.2e1;
t88 = t51 * pkin(5);
t58 = cos(qJ(3));
t87 = rSges(4,1) * t58;
t57 = sin(qJ(3));
t85 = rSges(4,2) * t57;
t83 = t52 * rSges(4,3);
t60 = t51 * rSges(5,3) + t72 * t52;
t8 = t51 * (-t52 * rSges(5,3) + t72 * t51) + t52 * t60;
t82 = t51 * rSges(4,3) + t52 * t87;
t81 = t49 + t50;
t80 = Icges(4,4) * t57;
t79 = Icges(4,4) * t58;
t78 = Icges(5,4) * t53;
t77 = Icges(5,4) * t54;
t33 = Icges(5,5) * t53 + Icges(5,6) * t54;
t63 = -Icges(5,2) * t53 + t77;
t65 = Icges(5,1) * t54 - t78;
t34 = Icges(5,2) * t54 + t78;
t35 = Icges(5,1) * t53 + t77;
t67 = -t34 * t53 + t35 * t54;
t76 = (t54 * (Icges(5,6) * t51 + t63 * t52) + t53 * (Icges(5,5) * t51 + t65 * t52) + t51 * t33 + t67 * t52) * t90 + (t54 * (-Icges(5,6) * t52 + t63 * t51) + t53 * (-Icges(5,5) * t52 + t65 * t51) - t52 * t33 + t67 * t51) * t89;
t61 = Icges(5,5) * t54 - Icges(5,6) * t53;
t16 = -Icges(5,3) * t52 + t61 * t51;
t17 = Icges(5,3) * t51 + t61 * t52;
t75 = -t52 * (t50 * t16 - t17 * t91) + t51 * (-t16 * t91 + t49 * t17);
t36 = t53 * rSges(5,1) + t54 * rSges(5,2);
t74 = -pkin(3) * t57 - t36;
t73 = -t85 + t87;
t66 = Icges(4,1) * t58 - t80;
t64 = -Icges(4,2) * t57 + t79;
t62 = Icges(4,5) * t58 - Icges(4,6) * t57;
t59 = -pkin(6) - pkin(5);
t48 = t58 * pkin(3) + pkin(2);
t47 = t52 * pkin(5);
t43 = t57 * rSges(4,1) + t58 * rSges(4,2);
t37 = t52 * t48;
t31 = t52 * rSges(3,1) - t51 * rSges(3,2);
t30 = -t51 * rSges(3,1) - t52 * rSges(3,2);
t25 = Icges(4,3) * t51 + t62 * t52;
t24 = -Icges(4,3) * t52 + t62 * t51;
t23 = t74 * t52;
t22 = t74 * t51;
t13 = t88 + (pkin(2) - t85) * t52 + t82;
t12 = t83 + t47 + (-pkin(2) - t73) * t51;
t11 = -t51 * t59 + t37 + t60;
t10 = (rSges(5,3) - t59) * t52 + (-t48 - t72) * t51;
t9 = t52 * (-t52 * t85 + t82) + (t73 * t51 - t83) * t51;
t3 = (t47 + (-pkin(2) + t48) * t51) * t51 + (-t52 * pkin(2) + t37 - t88) * t52 + t8;
t1 = [m(2) + m(3) + m(4) + m(5); 0; t54 * t34 + t53 * t35 + t58 * (Icges(4,2) * t58 + t80) + t57 * (Icges(4,1) * t57 + t79) + Icges(3,3) + m(3) * (t30 ^ 2 + t31 ^ 2) + m(4) * (t12 ^ 2 + t13 ^ 2) + m(5) * (t10 ^ 2 + t11 ^ 2); m(4) * t9 + m(5) * t3; (t58 * (Icges(4,6) * t51 + t64 * t52) + t57 * (Icges(4,5) * t51 + t66 * t52)) * t90 + (t58 * (-Icges(4,6) * t52 + t64 * t51) + t57 * (-Icges(4,5) * t52 + t66 * t51)) * t89 + m(5) * (t23 * t10 + t22 * t11) + m(4) * (-t12 * t52 - t13 * t51) * t43 + (t49 / 0.2e1 + t50 / 0.2e1) * (Icges(4,5) * t57 + Icges(4,6) * t58) + t76; m(4) * (t81 * t43 ^ 2 + t9 ^ 2) + t51 * (-t24 * t91 + t49 * t25) - t52 * (t50 * t24 - t25 * t91) + m(5) * (t22 ^ 2 + t23 ^ 2 + t3 ^ 2) + t75; m(5) * t8; m(5) * (-t10 * t52 - t11 * t51) * t36 + t76; m(5) * (t8 * t3 + (-t22 * t51 - t23 * t52) * t36) + t75; m(5) * (t81 * t36 ^ 2 + t8 ^ 2) + t75;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
