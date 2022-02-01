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
% m [6x1]
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
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:34:27
% EndTime: 2022-01-23 09:34:28
% DurationCPUTime: 0.38s
% Computational Cost: add. (1399->98), mult. (760->131), div. (0->0), fcn. (614->10), ass. (0->59)
t59 = qJ(1) + pkin(9);
t57 = qJ(3) + t59;
t54 = qJ(4) + t57;
t49 = sin(t54);
t87 = t49 ^ 2;
t50 = cos(t54);
t86 = t50 ^ 2;
t60 = sin(qJ(5));
t90 = Icges(6,5) * t60;
t62 = cos(qJ(5));
t89 = Icges(6,6) * t62;
t35 = t89 + t90;
t88 = t49 * t50;
t38 = t60 * rSges(6,1) + t62 * rSges(6,2);
t83 = m(6) * t38;
t52 = sin(t57);
t82 = pkin(3) * t52;
t61 = sin(qJ(1));
t81 = t61 * pkin(1);
t80 = rSges(6,1) * t62;
t79 = rSges(6,2) * t60;
t78 = t50 * rSges(6,3) + t49 * t79;
t56 = cos(t59);
t63 = cos(qJ(1));
t58 = t63 * pkin(1);
t77 = pkin(2) * t56 + t58;
t74 = Icges(6,2) * t62 ^ 2 + Icges(5,3) + (Icges(6,1) * t60 + 0.2e1 * Icges(6,4) * t62) * t60;
t73 = t35 * t86 + (t90 / 0.2e1 + t89 / 0.2e1 + t35 / 0.2e1) * t87;
t53 = cos(t57);
t29 = t53 * rSges(4,1) - t52 * rSges(4,2);
t25 = t50 * rSges(5,1) - t49 * rSges(5,2);
t72 = Icges(4,3) + t74;
t48 = pkin(3) * t53;
t23 = t25 + t48;
t55 = sin(t59);
t71 = -pkin(2) * t55 - t81;
t28 = -t52 * rSges(4,1) - t53 * rSges(4,2);
t24 = -t49 * rSges(5,1) - t50 * rSges(5,2);
t65 = Icges(6,5) * t62 - Icges(6,6) * t60;
t64 = t49 * rSges(6,3) + (-t79 + t80) * t50;
t11 = t50 * pkin(4) + t49 * pkin(8) + t64;
t22 = t24 - t82;
t9 = t11 + t48;
t10 = t50 * pkin(8) + (-pkin(4) - t80) * t49 + t78;
t8 = t10 - t82;
t40 = t63 * rSges(2,1) - t61 * rSges(2,2);
t39 = -t61 * rSges(2,1) - t63 * rSges(2,2);
t27 = t56 * rSges(3,1) - t55 * rSges(3,2) + t58;
t26 = -t55 * rSges(3,1) - t56 * rSges(3,2) - t81;
t21 = t29 + t77;
t20 = t28 + t71;
t15 = Icges(6,3) * t49 + t65 * t50;
t14 = -Icges(6,3) * t50 + t65 * t49;
t13 = t23 + t77;
t12 = t22 + t71;
t5 = t9 + t77;
t4 = t71 + t8;
t3 = t49 * (t49 * t80 - t78) + t50 * t64;
t1 = [Icges(2,3) + Icges(3,3) + m(6) * (t4 ^ 2 + t5 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(4) * (t20 ^ 2 + t21 ^ 2) + m(3) * (t26 ^ 2 + t27 ^ 2) + m(2) * (t39 ^ 2 + t40 ^ 2) + t72; 0; m(3) + m(4) + m(5) + m(6); m(6) * (t8 * t4 + t9 * t5) + m(5) * (t22 * t12 + t23 * t13) + m(4) * (t28 * t20 + t29 * t21) + t72; 0; m(4) * (t28 ^ 2 + t29 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2) + m(6) * (t8 ^ 2 + t9 ^ 2) + t72; m(6) * (t10 * t4 + t11 * t5) + m(5) * (t24 * t12 + t25 * t13) + t74; 0; m(5) * (t24 * t22 + t25 * t23) + m(6) * (t10 * t8 + t11 * t9) + t74; m(5) * (t24 ^ 2 + t25 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2) + t74; (-t4 * t50 - t49 * t5) * t83 + t73; m(6) * t3; (-t49 * t9 - t50 * t8) * t83 + t73; (-t10 * t50 - t11 * t49) * t83 + t73; m(6) * (t3 ^ 2 + (t86 + t87) * t38 ^ 2) + t49 * (-t14 * t88 + t87 * t15) - t50 * (t86 * t14 - t15 * t88);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
