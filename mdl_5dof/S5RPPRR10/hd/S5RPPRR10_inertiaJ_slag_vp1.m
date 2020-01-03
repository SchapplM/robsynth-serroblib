% Calculate joint inertia matrix for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR10_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR10_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR10_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:03:59
% EndTime: 2019-12-31 18:04:01
% DurationCPUTime: 0.62s
% Computational Cost: add. (1480->196), mult. (2354->288), div. (0->0), fcn. (2510->8), ass. (0->93)
t84 = cos(qJ(4));
t69 = t84 * pkin(4) + pkin(3);
t117 = pkin(3) - t69;
t80 = sin(pkin(8));
t82 = sin(qJ(4));
t110 = t80 * t82;
t102 = pkin(4) * t110;
t81 = cos(pkin(8));
t85 = cos(qJ(1));
t108 = t81 * t85;
t79 = qJ(4) + qJ(5);
t71 = sin(t79);
t72 = cos(t79);
t58 = -t81 * t71 + t80 * t72;
t48 = t58 * t85;
t90 = t80 * t71 + t81 * t72;
t49 = t90 * t85;
t83 = sin(qJ(1));
t23 = t49 * rSges(6,1) + t48 * rSges(6,2) - t83 * rSges(6,3);
t86 = -pkin(7) - pkin(6);
t116 = -t85 * t102 - t69 * t108 - t83 * t86 - t23;
t115 = t81 ^ 2;
t78 = t85 ^ 2;
t114 = m(4) / 0.2e1;
t113 = -t83 / 0.2e1;
t112 = t85 / 0.2e1;
t111 = -rSges(5,3) - pkin(6);
t109 = t81 * t82;
t36 = t58 * rSges(6,1) - rSges(6,2) * t90;
t107 = -pkin(4) * t109 - t117 * t80 + t36;
t61 = t80 * t84 - t109;
t55 = t61 * t85;
t89 = t81 * t84 + t110;
t56 = t89 * t85;
t106 = t56 * rSges(5,1) + t55 * rSges(5,2);
t105 = t85 * pkin(1) + t83 * qJ(2);
t104 = t83 ^ 2 + t78;
t103 = qJ(3) * t80;
t101 = pkin(3) * t108;
t46 = t58 * t83;
t47 = t90 * t83;
t18 = Icges(6,4) * t47 + Icges(6,2) * t46 + Icges(6,6) * t85;
t19 = Icges(6,4) * t49 + Icges(6,2) * t48 - Icges(6,6) * t83;
t20 = Icges(6,1) * t47 + Icges(6,4) * t46 + Icges(6,5) * t85;
t21 = Icges(6,1) * t49 + Icges(6,4) * t48 - Icges(6,5) * t83;
t33 = Icges(6,5) * t58 - Icges(6,6) * t90;
t34 = Icges(6,4) * t58 - Icges(6,2) * t90;
t35 = Icges(6,1) * t58 - Icges(6,4) * t90;
t99 = (-t19 * t90 + t58 * t21 - t83 * t33 + t48 * t34 + t49 * t35) * t113 + (-t18 * t90 + t58 * t20 + t85 * t33 + t46 * t34 + t47 * t35) * t112;
t16 = Icges(6,5) * t47 + Icges(6,6) * t46 + Icges(6,3) * t85;
t17 = Icges(6,5) * t49 + Icges(6,6) * t48 - Icges(6,3) * t83;
t1 = t85 * (-(t85 * t17 + t46 * t19 + t47 * t21) * t83 + (t85 * t16 + t46 * t18 + t47 * t20) * t85);
t2 = -(-t83 * t17 + t48 * t19 + t49 * t21) * t83 + (-t83 * t16 + t48 * t18 + t49 * t20) * t85;
t98 = -t83 * t2 + t1;
t97 = t104 * t80;
t96 = pkin(2) * t108 + t85 * t103 + t105;
t95 = t114 + m(5) / 0.2e1 + m(6) / 0.2e1;
t94 = rSges(3,1) * t81 - rSges(3,2) * t80;
t53 = t61 * t83;
t54 = t89 * t83;
t93 = t54 * rSges(5,1) + t53 * rSges(5,2);
t92 = -t47 * rSges(6,1) - t46 * rSges(6,2);
t14 = t107 * t83;
t15 = t107 * t85;
t91 = t14 * t83 + t15 * t85;
t10 = t96 - t116;
t74 = t85 * qJ(2);
t9 = t74 + (-rSges(6,3) + t86) * t85 + (-pkin(1) + (-pkin(2) - t69) * t81 + (-pkin(4) * t82 - qJ(3)) * t80) * t83 + t92;
t88 = m(6) * (t10 * t83 + t85 * t9);
t12 = t74 + t111 * t85 + (-t103 - pkin(1) + (-pkin(2) - pkin(3)) * t81) * t83 - t93;
t13 = t111 * t83 + t101 + t106 + t96;
t87 = m(5) * (t12 * t85 + t13 * t83);
t65 = t85 * rSges(2,1) - t83 * rSges(2,2);
t64 = -t83 * rSges(2,1) - t85 * rSges(2,2);
t43 = t83 * rSges(3,3) + t94 * t85 + t105;
t42 = t85 * rSges(3,3) + t74 + (-pkin(1) - t94) * t83;
t41 = t61 * rSges(5,1) - rSges(5,2) * t89;
t40 = Icges(5,1) * t61 - Icges(5,4) * t89;
t39 = Icges(5,4) * t61 - Icges(5,2) * t89;
t38 = Icges(5,5) * t61 - Icges(5,6) * t89;
t31 = t83 * rSges(4,2) + (rSges(4,1) * t81 + rSges(4,3) * t80) * t85 + t96;
t30 = t85 * rSges(4,2) + t74 + (-pkin(1) + (-rSges(4,1) - pkin(2)) * t81 + (-rSges(4,3) - qJ(3)) * t80) * t83;
t29 = Icges(5,1) * t56 + Icges(5,4) * t55 - Icges(5,5) * t83;
t28 = Icges(5,1) * t54 + Icges(5,4) * t53 + Icges(5,5) * t85;
t27 = Icges(5,4) * t56 + Icges(5,2) * t55 - Icges(5,6) * t83;
t26 = Icges(5,4) * t54 + Icges(5,2) * t53 + Icges(5,6) * t85;
t25 = Icges(5,5) * t56 + Icges(5,6) * t55 - Icges(5,3) * t83;
t24 = Icges(5,5) * t54 + Icges(5,6) * t53 + Icges(5,3) * t85;
t22 = t85 * rSges(6,3) - t92;
t11 = -t85 * t106 - t83 * t93;
t8 = -t83 * t22 - t85 * t23;
t3 = (t101 + t116) * t85 + (t85 * t86 - t22 + (t117 * t81 - t102) * t83) * t83;
t4 = [-t90 * t34 + t58 * t35 - t89 * t39 + t61 * t40 + Icges(2,3) + (Icges(3,2) + Icges(4,3)) * t115 + ((Icges(3,1) + Icges(4,1)) * t80 + 0.2e1 * (Icges(3,4) - Icges(4,5)) * t81) * t80 + m(6) * (t10 ^ 2 + t9 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(4) * (t30 ^ 2 + t31 ^ 2) + m(3) * (t42 ^ 2 + t43 ^ 2) + m(2) * (t64 ^ 2 + t65 ^ 2); m(6) * (-t85 * t10 + t83 * t9) + m(5) * (t83 * t12 - t85 * t13) + m(4) * (t83 * t30 - t85 * t31) + m(3) * (t83 * t42 - t85 * t43); 0.2e1 * (m(3) / 0.2e1 + t95) * t104; 0.2e1 * (t88 / 0.2e1 + t87 / 0.2e1 + (t30 * t85 + t31 * t83) * t114) * t80; 0; 0.2e1 * t95 * (t104 * t80 ^ 2 + t115); m(6) * (t14 * t10 + t15 * t9) + t41 * t87 + t99 + (-t27 * t89 + t61 * t29 - t83 * t38 + t55 * t39 + t56 * t40) * t113 + (-t26 * t89 + t61 * t28 + t85 * t38 + t53 * t39 + t54 * t40) * t112; m(6) * (-t14 * t85 + t15 * t83); m(5) * (-t11 * t81 + t41 * t97) + m(6) * (-t3 * t81 + t91 * t80); (t85 * t24 + t53 * t26 + t54 * t28) * t78 + t1 + (-t2 + (-t83 * t25 + t55 * t27 + t56 * t29) * t83 + (t83 * t24 - t85 * t25 - t55 * t26 - t53 * t27 - t56 * t28 - t54 * t29) * t85) * t83 + m(5) * (t104 * t41 ^ 2 + t11 ^ 2) + m(6) * (t14 ^ 2 + t15 ^ 2 + t3 ^ 2); t36 * t88 + t99; 0; m(6) * (t36 * t97 - t8 * t81); m(6) * (t8 * t3 + t91 * t36) + t98; m(6) * (t104 * t36 ^ 2 + t8 ^ 2) + t98;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;
