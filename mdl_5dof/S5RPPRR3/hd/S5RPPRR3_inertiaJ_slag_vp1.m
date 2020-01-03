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
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:27:28
% EndTime: 2020-01-03 11:27:32
% DurationCPUTime: 0.62s
% Computational Cost: add. (1611->142), mult. (1054->208), div. (0->0), fcn. (916->10), ass. (0->77)
t68 = qJ(1) + pkin(8);
t59 = sin(t68);
t61 = cos(t68);
t106 = t61 * t59;
t67 = pkin(9) + qJ(4);
t62 = qJ(5) + t67;
t53 = sin(t62);
t54 = cos(t62);
t87 = -rSges(6,1) * t54 + rSges(6,2) * t53;
t56 = t59 ^ 2;
t57 = t61 ^ 2;
t105 = -t59 / 0.2e1;
t104 = -t61 / 0.2e1;
t70 = cos(pkin(9));
t55 = t70 * pkin(3) + pkin(2);
t60 = cos(t67);
t103 = rSges(5,1) * t60;
t58 = sin(t67);
t101 = rSges(5,2) * t58;
t71 = -pkin(6) - qJ(3);
t99 = rSges(5,3) - t71;
t46 = pkin(4) * t60 + t55;
t66 = -pkin(7) + t71;
t98 = t59 * t46 + t61 * t66;
t97 = t57 + t56;
t96 = Icges(5,4) * t58;
t95 = Icges(5,4) * t60;
t94 = Icges(6,4) * t53;
t93 = Icges(6,4) * t54;
t92 = rSges(4,3) + qJ(3);
t36 = Icges(6,5) * t53 + Icges(6,6) * t54;
t78 = -Icges(6,2) * t53 + t93;
t80 = Icges(6,1) * t54 - t94;
t37 = Icges(6,2) * t54 + t94;
t38 = Icges(6,1) * t53 + t93;
t82 = t37 * t53 - t38 * t54;
t91 = (t54 * (-Icges(6,6) * t59 - t78 * t61) + t53 * (-Icges(6,5) * t59 - t80 * t61) - t59 * t36 + t82 * t61) * t105 + (t54 * (-Icges(6,6) * t61 + t78 * t59) + t53 * (-Icges(6,5) * t61 + t80 * t59) - t61 * t36 - t82 * t59) * t104;
t39 = t53 * rSges(6,1) + t54 * rSges(6,2);
t90 = pkin(4) * t58 + t39;
t76 = Icges(6,5) * t54 - Icges(6,6) * t53;
t19 = -Icges(6,3) * t61 + t76 * t59;
t20 = -Icges(6,3) * t59 - t76 * t61;
t89 = -t61 * (t20 * t106 + t57 * t19) - t59 * (t19 * t106 + t56 * t20);
t88 = -t101 + t103;
t81 = Icges(5,1) * t60 - t96;
t79 = -Icges(5,2) * t58 + t95;
t77 = Icges(5,5) * t60 - Icges(5,6) * t58;
t69 = sin(pkin(9));
t75 = rSges(4,1) * t70 - rSges(4,2) * t69 + pkin(2);
t74 = -t61 * rSges(6,3) - t87 * t59;
t73 = cos(qJ(1));
t72 = sin(qJ(1));
t65 = t73 * pkin(1);
t64 = t72 * pkin(1);
t51 = t73 * rSges(2,1) - t72 * rSges(2,2);
t50 = t72 * rSges(2,1) + t73 * rSges(2,2);
t49 = t59 * t103;
t48 = t61 * t55;
t44 = t58 * rSges(5,1) + t60 * rSges(5,2);
t33 = t61 * rSges(3,1) - t59 * rSges(3,2) + t65;
t32 = t59 * rSges(3,1) + t61 * rSges(3,2) + t64;
t27 = -Icges(5,3) * t59 - t77 * t61;
t26 = -Icges(5,3) * t61 + t77 * t59;
t25 = -t59 * rSges(6,3) + t87 * t61;
t18 = t90 * t61;
t17 = t90 * t59;
t16 = t59 * t74;
t15 = t92 * t59 + t75 * t61 + t65;
t14 = t75 * t59 - t92 * t61 + t64;
t13 = t99 * t59 + t88 * t61 + t48 + t65;
t12 = t49 + t64 - t99 * t61 + (t55 - t101) * t59;
t11 = t65 + (rSges(6,3) - t66) * t59 + (t46 - t87) * t61;
t10 = t64 + t74 + t98;
t9 = t59 * (-t59 * t101 + t49) + t88 * t57;
t8 = -t61 * t25 + t16;
t3 = t59 * (-t59 * t55 + t98) + t16 + (t61 * t46 - t59 * t66 - t25 - t48) * t61;
t1 = [Icges(4,2) * t70 ^ 2 + t54 * t37 + t53 * t38 + t60 * (Icges(5,2) * t60 + t96) + t58 * (Icges(5,1) * t58 + t95) + Icges(2,3) + Icges(3,3) + (Icges(4,1) * t69 + 0.2e1 * Icges(4,4) * t70) * t69 + m(2) * (t50 ^ 2 + t51 ^ 2) + m(3) * (t32 ^ 2 + t33 ^ 2) + m(4) * (t14 ^ 2 + t15 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2); 0; m(3) + m(4) + m(5) + m(6); m(4) * (-t59 * t14 - t61 * t15) + m(5) * (-t59 * t12 - t61 * t13) + m(6) * (-t59 * t10 - t61 * t11); 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t97; (t60 * (-Icges(5,6) * t61 + t79 * t59) + t58 * (-Icges(5,5) * t61 + t81 * t59)) * t104 + (t60 * (-Icges(5,6) * t59 - t79 * t61) + t58 * (-Icges(5,5) * t59 - t81 * t61)) * t105 + m(6) * (t18 * t10 - t17 * t11) + m(5) * (t12 * t61 - t13 * t59) * t44 + (t57 / 0.2e1 + t56 / 0.2e1) * (Icges(5,5) * t58 + Icges(5,6) * t60) + t91; m(5) * t9 + m(6) * t3; m(6) * (t17 * t61 - t18 * t59); m(5) * (t97 * t44 ^ 2 + t9 ^ 2) - t61 * (t27 * t106 + t57 * t26) - t59 * (t26 * t106 + t56 * t27) + m(6) * (t17 ^ 2 + t18 ^ 2 + t3 ^ 2) + t89; m(6) * (t10 * t61 - t11 * t59) * t39 + t91; m(6) * t8; 0; m(6) * (t8 * t3 + (t17 * t59 + t18 * t61) * t39) + t89; m(6) * (t97 * t39 ^ 2 + t8 ^ 2) + t89;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
