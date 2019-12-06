% Calculate joint inertia matrix for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:37:59
% EndTime: 2019-12-05 17:38:01
% DurationCPUTime: 0.52s
% Computational Cost: add. (749->137), mult. (1114->204), div. (0->0), fcn. (942->6), ass. (0->73)
t66 = sin(qJ(1));
t68 = cos(qJ(1));
t107 = t66 * t68;
t64 = qJ(4) + qJ(5);
t56 = sin(t64);
t57 = cos(t64);
t84 = rSges(6,1) * t56 + rSges(6,2) * t57;
t62 = t66 ^ 2;
t63 = t68 ^ 2;
t106 = -t66 / 0.2e1;
t105 = t68 / 0.2e1;
t51 = t62 + t63;
t40 = m(5) * t51;
t39 = m(6) * t51;
t65 = sin(qJ(4));
t104 = pkin(4) * t65;
t103 = -rSges(5,3) - pkin(6);
t67 = cos(qJ(4));
t101 = rSges(5,2) * t67;
t99 = t65 * t68;
t98 = -pkin(1) - qJ(3);
t97 = rSges(5,1) * t99 + t68 * t101;
t96 = t68 * pkin(1) + t66 * qJ(2);
t95 = Icges(5,4) * t65;
t94 = Icges(5,4) * t67;
t93 = Icges(6,4) * t56;
t92 = Icges(6,4) * t57;
t91 = m(4) * t51 + t39 + t40;
t23 = -t66 * rSges(6,3) + t84 * t68;
t69 = -pkin(7) - pkin(6);
t90 = -pkin(4) * t99 - t66 * t69 - t23;
t89 = t68 * qJ(3) + t96;
t35 = Icges(6,5) * t57 - Icges(6,6) * t56;
t74 = Icges(6,2) * t57 + t93;
t76 = Icges(6,1) * t56 + t92;
t36 = -Icges(6,2) * t56 + t92;
t37 = Icges(6,1) * t57 - t93;
t78 = t36 * t57 + t37 * t56;
t88 = (-t56 * (-Icges(6,6) * t66 + t74 * t68) + t57 * (-Icges(6,5) * t66 + t76 * t68) - t66 * t35 + t78 * t68) * t106 + (-t56 * (Icges(6,6) * t68 + t74 * t66) + t57 * (Icges(6,5) * t68 + t76 * t66) + t68 * t35 + t78 * t66) * t105;
t72 = Icges(6,5) * t56 + Icges(6,6) * t57;
t16 = Icges(6,3) * t68 + t72 * t66;
t17 = -Icges(6,3) * t66 + t72 * t68;
t87 = -t66 * (-t16 * t107 + t62 * t17) + t68 * (-t17 * t107 + t63 * t16);
t38 = t57 * rSges(6,1) - t56 * rSges(6,2);
t86 = pkin(4) * t67 + t38;
t85 = -rSges(5,1) * t65 - t101;
t14 = t86 * t66;
t15 = t86 * t68;
t83 = t14 * t66 + t15 * t68;
t77 = Icges(5,1) * t65 + t94;
t75 = Icges(5,2) * t67 + t95;
t73 = Icges(5,5) * t65 + Icges(5,6) * t67;
t10 = t89 - t90;
t60 = t68 * qJ(2);
t9 = t60 + (-rSges(6,3) + t69) * t68 + (-t84 + t98 - t104) * t66;
t71 = m(6) * (t66 * t10 + t68 * t9);
t12 = t60 + t103 * t68 + (t85 + t98) * t66;
t13 = t103 * t66 + t89 + t97;
t70 = m(5) * (t68 * t12 + t66 * t13);
t48 = t68 * rSges(2,1) - t66 * rSges(2,2);
t47 = t67 * rSges(5,1) - t65 * rSges(5,2);
t46 = -t66 * rSges(2,1) - t68 * rSges(2,2);
t33 = -t68 * rSges(3,2) + t66 * rSges(3,3) + t96;
t32 = t68 * rSges(3,3) + t60 + (rSges(3,2) - pkin(1)) * t66;
t27 = -Icges(5,3) * t66 + t73 * t68;
t26 = Icges(5,3) * t68 + t73 * t66;
t25 = t66 * rSges(4,2) + t68 * rSges(4,3) + t89;
t24 = t68 * rSges(4,2) + t60 + (-rSges(4,3) + t98) * t66;
t22 = t68 * rSges(6,3) + t84 * t66;
t11 = t85 * t62 - t68 * t97;
t8 = -t66 * t22 - t68 * t23;
t3 = t90 * t68 + (-t66 * t104 + t68 * t69 - t22) * t66;
t1 = [-t56 * t36 + t57 * t37 - t65 * (-Icges(5,2) * t65 + t94) + t67 * (Icges(5,1) * t67 - t95) + Icges(3,1) + Icges(4,1) + Icges(2,3) + m(6) * (t10 ^ 2 + t9 ^ 2) + m(2) * (t46 ^ 2 + t48 ^ 2) + m(3) * (t32 ^ 2 + t33 ^ 2) + m(4) * (t24 ^ 2 + t25 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2); m(6) * (-t68 * t10 + t66 * t9) + m(3) * (t66 * t32 - t68 * t33) + m(4) * (t66 * t24 - t68 * t25) + m(5) * (t66 * t12 - t68 * t13); m(3) * t51 + t91; t71 + m(4) * (t68 * t24 + t66 * t25) + t70; 0; t91; (-t65 * (-Icges(5,6) * t66 + t75 * t68) + t67 * (-Icges(5,5) * t66 + t77 * t68)) * t106 + (-t65 * (Icges(5,6) * t68 + t75 * t66) + t67 * (Icges(5,5) * t68 + t77 * t66)) * t105 + m(6) * (t14 * t10 + t15 * t9) + t47 * t70 + (t62 / 0.2e1 + t63 / 0.2e1) * (Icges(5,5) * t67 - Icges(5,6) * t65) + t88; m(6) * (-t14 * t68 + t15 * t66); m(6) * t83 + t47 * t40; m(5) * (t51 * t47 ^ 2 + t11 ^ 2) - t66 * (-t26 * t107 + t62 * t27) + t68 * (-t27 * t107 + t63 * t26) + m(6) * (t14 ^ 2 + t15 ^ 2 + t3 ^ 2) + t87; t38 * t71 + t88; 0; t38 * t39; m(6) * (t8 * t3 + t83 * t38) + t87; m(6) * (t51 * t38 ^ 2 + t8 ^ 2) + t87;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
