% Calculate joint inertia matrix for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:41:38
% EndTime: 2019-12-05 17:41:40
% DurationCPUTime: 0.54s
% Computational Cost: add. (1611->147), mult. (1054->210), div. (0->0), fcn. (916->10), ass. (0->80)
t66 = qJ(1) + pkin(8);
t59 = sin(t66);
t61 = cos(t66);
t108 = t61 * t59;
t56 = t59 ^ 2;
t57 = t61 ^ 2;
t107 = t59 / 0.2e1;
t106 = t61 / 0.2e1;
t65 = pkin(9) + qJ(4);
t62 = qJ(5) + t65;
t53 = sin(t62);
t54 = cos(t62);
t73 = Icges(6,5) * t54 - Icges(6,6) * t53;
t19 = Icges(6,3) * t61 - t73 * t59;
t20 = Icges(6,3) * t59 + t73 * t61;
t105 = t59 * (t19 * t108 + t56 * t20) + t61 * (t20 * t108 + t57 * t19);
t70 = sin(qJ(1));
t104 = t70 * pkin(1);
t71 = cos(qJ(1));
t103 = t71 * pkin(1);
t68 = cos(pkin(9));
t55 = t68 * pkin(3) + pkin(2);
t60 = cos(t65);
t102 = rSges(5,1) * t60;
t101 = rSges(6,1) * t54;
t58 = sin(t65);
t100 = rSges(5,2) * t58;
t99 = rSges(6,2) * t53;
t98 = t59 * rSges(5,3);
t69 = -pkin(6) - qJ(3);
t97 = t61 * t69;
t96 = t61 * rSges(6,3) + t59 * t99;
t45 = pkin(4) * t60 + t55;
t95 = t45 - t55;
t94 = t61 * rSges(5,3) + t59 * t100;
t93 = t57 + t56;
t92 = Icges(5,4) * t58;
t91 = Icges(5,4) * t60;
t90 = Icges(6,4) * t53;
t89 = Icges(6,4) * t54;
t88 = rSges(4,3) + qJ(3);
t35 = Icges(6,5) * t53 + Icges(6,6) * t54;
t75 = -Icges(6,2) * t53 + t89;
t77 = Icges(6,1) * t54 - t90;
t36 = Icges(6,2) * t54 + t90;
t37 = Icges(6,1) * t53 + t89;
t79 = t36 * t53 - t37 * t54;
t87 = (t54 * (Icges(6,6) * t59 + t75 * t61) + t53 * (Icges(6,5) * t59 + t77 * t61) + t59 * t35 - t79 * t61) * t107 + (t54 * (Icges(6,6) * t61 - t75 * t59) + t53 * (Icges(6,5) * t61 - t77 * t59) + t61 * t35 + t79 * t59) * t106;
t38 = t53 * rSges(6,1) + t54 * rSges(6,2);
t86 = pkin(4) * t58 + t38;
t85 = -t100 + t102;
t84 = -t99 + t101;
t78 = Icges(5,1) * t60 - t92;
t76 = -Icges(5,2) * t58 + t91;
t74 = Icges(5,5) * t60 - Icges(5,6) * t58;
t67 = sin(pkin(9));
t72 = -rSges(4,1) * t68 + rSges(4,2) * t67 - pkin(2);
t64 = -pkin(7) + t69;
t50 = t59 * t69;
t49 = -t71 * rSges(2,1) + t70 * rSges(2,2);
t48 = -t70 * rSges(2,1) - t71 * rSges(2,2);
t43 = t58 * rSges(5,1) + t60 * rSges(5,2);
t33 = -t61 * rSges(3,1) + t59 * rSges(3,2) - t103;
t32 = -t59 * rSges(3,1) - t61 * rSges(3,2) - t104;
t27 = Icges(5,3) * t59 + t74 * t61;
t26 = Icges(5,3) * t61 - t74 * t59;
t25 = -t59 * t101 + t96;
t18 = t86 * t61;
t17 = t86 * t59;
t16 = t61 * (t59 * rSges(6,3) + t84 * t61);
t15 = -t88 * t59 + t72 * t61 - t103;
t14 = t72 * t59 + t88 * t61 - t104;
t13 = -t98 - t103 + t50 + (-t55 - t85) * t61;
t12 = -t104 - t97 + (-t55 - t102) * t59 + t94;
t11 = -t103 + (-rSges(6,3) + t64) * t59 + (-t45 - t84) * t61;
t10 = -t104 - t61 * t64 + (-t45 - t101) * t59 + t96;
t9 = -t59 * (-t59 * t102 + t94) + (t85 * t61 + t98) * t61;
t8 = -t59 * t25 + t16;
t3 = t16 + (t95 * t61 + t50) * t61 + (t95 * t59 - t25 - t97) * t59;
t1 = [Icges(4,2) * t68 ^ 2 + t54 * t36 + t53 * t37 + t60 * (Icges(5,2) * t60 + t92) + t58 * (Icges(5,1) * t58 + t91) + Icges(2,3) + Icges(3,3) + (Icges(4,1) * t67 + 0.2e1 * Icges(4,4) * t68) * t67 + m(2) * (t48 ^ 2 + t49 ^ 2) + m(3) * (t32 ^ 2 + t33 ^ 2) + m(4) * (t14 ^ 2 + t15 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2); 0; m(3) + m(4) + m(5) + m(6); m(4) * (t59 * t14 + t61 * t15) + m(5) * (t59 * t12 + t61 * t13) + m(6) * (t59 * t10 + t61 * t11); 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t93; (t60 * (Icges(5,6) * t61 - t76 * t59) + t58 * (Icges(5,5) * t61 - t78 * t59)) * t106 + (t60 * (Icges(5,6) * t59 + t76 * t61) + t58 * (Icges(5,5) * t59 + t78 * t61)) * t107 + m(6) * (-t18 * t10 + t17 * t11) + m(5) * (-t12 * t61 + t13 * t59) * t43 + (t57 / 0.2e1 + t56 / 0.2e1) * (Icges(5,5) * t58 + Icges(5,6) * t60) + t87; m(5) * t9 + m(6) * t3; m(6) * (t17 * t61 - t18 * t59); m(5) * (t93 * t43 ^ 2 + t9 ^ 2) + t61 * (t108 * t27 + t57 * t26) + t59 * (t108 * t26 + t56 * t27) + m(6) * (t17 ^ 2 + t18 ^ 2 + t3 ^ 2) + t105; m(6) * (-t10 * t61 + t11 * t59) * t38 + t87; m(6) * t8; 0; m(6) * (t8 * t3 + (t17 * t59 + t18 * t61) * t38) + t105; m(6) * (t93 * t38 ^ 2 + t8 ^ 2) + t105;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
