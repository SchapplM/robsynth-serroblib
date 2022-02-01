% Calculate joint inertia matrix for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:32
% EndTime: 2022-01-20 11:30:32
% DurationCPUTime: 0.41s
% Computational Cost: add. (1570->110), mult. (862->146), div. (0->0), fcn. (686->10), ass. (0->63)
t64 = qJ(1) + qJ(2);
t62 = qJ(3) + t64;
t57 = pkin(9) + t62;
t52 = sin(t57);
t91 = t52 ^ 2;
t53 = cos(t57);
t90 = t53 ^ 2;
t65 = sin(qJ(5));
t94 = Icges(6,5) * t65;
t67 = cos(qJ(5));
t93 = Icges(6,6) * t67;
t39 = t93 + t94;
t92 = t52 * t53;
t42 = t65 * rSges(6,1) + t67 * rSges(6,2);
t87 = m(6) * t42;
t60 = sin(t64);
t86 = pkin(2) * t60;
t58 = sin(t62);
t85 = pkin(3) * t58;
t66 = sin(qJ(1));
t84 = t66 * pkin(1);
t83 = rSges(6,1) * t67;
t82 = rSges(6,2) * t65;
t81 = t53 * rSges(6,3) + t52 * t82;
t78 = t39 * t90 + (t94 / 0.2e1 + t93 / 0.2e1 + t39 / 0.2e1) * t91;
t61 = cos(t64);
t33 = t61 * rSges(3,1) - t60 * rSges(3,2);
t59 = cos(t62);
t31 = t59 * rSges(4,1) - t58 * rSges(4,2);
t77 = Icges(6,2) * t67 ^ 2 + Icges(4,3) + Icges(5,3) + (Icges(6,1) * t65 + 0.2e1 * Icges(6,4) * t67) * t65;
t56 = pkin(2) * t61;
t27 = t31 + t56;
t54 = pkin(3) * t59;
t23 = t53 * rSges(5,1) - t52 * rSges(5,2) + t54;
t32 = -t60 * rSges(3,1) - t61 * rSges(3,2);
t30 = -t58 * rSges(4,1) - t59 * rSges(4,2);
t73 = Icges(3,3) + t77;
t70 = Icges(6,5) * t67 - Icges(6,6) * t65;
t69 = t52 * rSges(6,3) + (-t82 + t83) * t53;
t21 = t23 + t56;
t26 = t30 - t86;
t22 = -t52 * rSges(5,1) - t53 * rSges(5,2) - t85;
t11 = t53 * pkin(4) + t52 * pkin(8) + t54 + t69;
t9 = t11 + t56;
t20 = t22 - t86;
t10 = -t85 + t53 * pkin(8) + (-pkin(4) - t83) * t52 + t81;
t8 = t10 - t86;
t68 = cos(qJ(1));
t63 = t68 * pkin(1);
t44 = t68 * rSges(2,1) - t66 * rSges(2,2);
t43 = -t66 * rSges(2,1) - t68 * rSges(2,2);
t29 = t33 + t63;
t28 = t32 - t84;
t25 = t27 + t63;
t24 = t26 - t84;
t15 = Icges(6,3) * t52 + t70 * t53;
t14 = -Icges(6,3) * t53 + t70 * t52;
t13 = t21 + t63;
t12 = t20 - t84;
t7 = t63 + t9;
t6 = t8 - t84;
t3 = t52 * (t52 * t83 - t81) + t53 * t69;
t1 = [Icges(2,3) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(6) * (t6 ^ 2 + t7 ^ 2) + m(4) * (t24 ^ 2 + t25 ^ 2) + m(3) * (t28 ^ 2 + t29 ^ 2) + m(2) * (t43 ^ 2 + t44 ^ 2) + t73; m(5) * (t20 * t12 + t21 * t13) + m(6) * (t8 * t6 + t9 * t7) + m(4) * (t26 * t24 + t27 * t25) + m(3) * (t32 * t28 + t33 * t29) + t73; m(6) * (t8 ^ 2 + t9 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t26 ^ 2 + t27 ^ 2) + m(3) * (t32 ^ 2 + t33 ^ 2) + t73; m(5) * (t22 * t12 + t23 * t13) + m(6) * (t10 * t6 + t11 * t7) + m(4) * (t30 * t24 + t31 * t25) + t77; m(6) * (t10 * t8 + t11 * t9) + m(5) * (t22 * t20 + t23 * t21) + m(4) * (t30 * t26 + t31 * t27) + t77; m(4) * (t30 ^ 2 + t31 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2) + t77; 0; 0; 0; m(5) + m(6); (-t52 * t7 - t53 * t6) * t87 + t78; (-t52 * t9 - t53 * t8) * t87 + t78; (-t10 * t53 - t11 * t52) * t87 + t78; m(6) * t3; m(6) * (t3 ^ 2 + (t90 + t91) * t42 ^ 2) + t52 * (-t14 * t92 + t91 * t15) - t53 * (t90 * t14 - t15 * t92);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
