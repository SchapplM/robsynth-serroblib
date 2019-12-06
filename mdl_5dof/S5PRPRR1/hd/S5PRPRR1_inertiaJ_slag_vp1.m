% Calculate joint inertia matrix for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:42:38
% EndTime: 2019-12-05 15:42:41
% DurationCPUTime: 0.54s
% Computational Cost: add. (1587->131), mult. (1024->199), div. (0->0), fcn. (892->8), ass. (0->67)
t64 = pkin(8) + qJ(2);
t57 = sin(t64);
t59 = cos(t64);
t99 = t57 * t59;
t63 = pkin(9) + qJ(4);
t60 = qJ(5) + t63;
t51 = sin(t60);
t52 = cos(t60);
t81 = rSges(6,1) * t52 - rSges(6,2) * t51;
t54 = t57 ^ 2;
t55 = t59 ^ 2;
t98 = t57 / 0.2e1;
t97 = -t59 / 0.2e1;
t66 = cos(pkin(9));
t53 = t66 * pkin(3) + pkin(2);
t58 = cos(t63);
t96 = rSges(5,1) * t58;
t56 = sin(t63);
t94 = rSges(5,2) * t56;
t67 = -pkin(6) - qJ(3);
t69 = t57 * rSges(6,3) + t81 * t59;
t8 = t57 * (-t59 * rSges(6,3) + t81 * t57) + t59 * t69;
t92 = t57 * rSges(5,3) + t59 * t96;
t91 = t54 + t55;
t90 = Icges(5,4) * t56;
t89 = Icges(5,4) * t58;
t88 = Icges(6,4) * t51;
t87 = Icges(6,4) * t52;
t86 = rSges(4,3) + qJ(3);
t34 = Icges(6,5) * t51 + Icges(6,6) * t52;
t72 = -Icges(6,2) * t51 + t87;
t74 = Icges(6,1) * t52 - t88;
t35 = Icges(6,2) * t52 + t88;
t36 = Icges(6,1) * t51 + t87;
t76 = -t35 * t51 + t36 * t52;
t85 = (t52 * (Icges(6,6) * t57 + t72 * t59) + t51 * (Icges(6,5) * t57 + t74 * t59) + t57 * t34 + t76 * t59) * t98 + (t52 * (-Icges(6,6) * t59 + t72 * t57) + t51 * (-Icges(6,5) * t59 + t74 * t57) - t59 * t34 + t76 * t57) * t97;
t70 = Icges(6,5) * t52 - Icges(6,6) * t51;
t20 = -Icges(6,3) * t59 + t70 * t57;
t21 = Icges(6,3) * t57 + t70 * t59;
t84 = -t59 * (t55 * t20 - t21 * t99) + t57 * (-t20 * t99 + t54 * t21);
t37 = t51 * rSges(6,1) + t52 * rSges(6,2);
t83 = -pkin(4) * t56 - t37;
t82 = -t94 + t96;
t75 = Icges(5,1) * t58 - t90;
t73 = -Icges(5,2) * t56 + t89;
t71 = Icges(5,5) * t58 - Icges(5,6) * t56;
t65 = sin(pkin(9));
t68 = rSges(4,1) * t66 - rSges(4,2) * t65 + pkin(2);
t62 = -pkin(7) + t67;
t46 = pkin(4) * t58 + t53;
t44 = t59 * rSges(3,1) - t57 * rSges(3,2);
t43 = -t57 * rSges(3,1) - t59 * rSges(3,2);
t42 = t56 * rSges(5,1) + t58 * rSges(5,2);
t32 = t59 * t46;
t27 = Icges(5,3) * t57 + t71 * t59;
t26 = -Icges(5,3) * t59 + t71 * t57;
t19 = t83 * t59;
t18 = t83 * t57;
t15 = t86 * t57 + t68 * t59;
t14 = -t68 * t57 + t86 * t59;
t13 = -t57 * t67 + (t53 - t94) * t59 + t92;
t12 = (rSges(5,3) - t67) * t59 + (-t53 - t82) * t57;
t11 = -t57 * t62 + t32 + t69;
t10 = (rSges(6,3) - t62) * t59 + (-t46 - t81) * t57;
t9 = t59 * (-t59 * t94 + t92) + (-t59 * rSges(5,3) + t82 * t57) * t57;
t3 = t59 * (-t59 * t53 + t32) + (t46 - t53) * t54 + t8;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; Icges(4,2) * t66 ^ 2 + t52 * t35 + t51 * t36 + t58 * (Icges(5,2) * t58 + t90) + t56 * (Icges(5,1) * t56 + t89) + Icges(3,3) + (Icges(4,1) * t65 + 0.2e1 * Icges(4,4) * t66) * t65 + m(6) * (t10 ^ 2 + t11 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(4) * (t14 ^ 2 + t15 ^ 2) + m(3) * (t43 ^ 2 + t44 ^ 2); 0; m(6) * (t57 * t10 - t59 * t11) + m(5) * (t57 * t12 - t59 * t13) + m(4) * (t57 * t14 - t59 * t15); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t91; m(5) * t9 + m(6) * t3; (t58 * (Icges(5,6) * t57 + t73 * t59) + t56 * (Icges(5,5) * t57 + t75 * t59)) * t98 + (t58 * (-Icges(5,6) * t59 + t73 * t57) + t56 * (-Icges(5,5) * t59 + t75 * t57)) * t97 + m(6) * (t19 * t10 + t18 * t11) + m(5) * (-t12 * t59 - t13 * t57) * t42 + (t55 / 0.2e1 + t54 / 0.2e1) * (Icges(5,5) * t56 + Icges(5,6) * t58) + t85; m(6) * (-t18 * t59 + t19 * t57); m(5) * (t91 * t42 ^ 2 + t9 ^ 2) + t57 * (-t26 * t99 + t54 * t27) - t59 * (t55 * t26 - t27 * t99) + m(6) * (t18 ^ 2 + t19 ^ 2 + t3 ^ 2) + t84; m(6) * t8; m(6) * (-t10 * t59 - t11 * t57) * t37 + t85; 0; m(6) * (t8 * t3 + (-t18 * t57 - t19 * t59) * t37) + t84; m(6) * (t91 * t37 ^ 2 + t8 ^ 2) + t84;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
