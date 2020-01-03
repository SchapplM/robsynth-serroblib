% Calculate joint inertia matrix for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR9_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR9_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR9_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:09
% EndTime: 2019-12-31 16:56:11
% DurationCPUTime: 0.79s
% Computational Cost: add. (866->158), mult. (2093->258), div. (0->0), fcn. (2202->6), ass. (0->88)
t76 = cos(qJ(3));
t117 = Icges(4,5) * t76;
t73 = sin(qJ(3));
t116 = Icges(4,6) * t73;
t115 = t117 / 0.2e1 - t116 / 0.2e1;
t74 = sin(qJ(1));
t70 = t74 ^ 2;
t77 = cos(qJ(1));
t71 = t77 ^ 2;
t96 = t70 + t71;
t114 = (rSges(4,1) * t73 + rSges(4,2) * t76) * t77;
t111 = t77 / 0.2e1;
t110 = -pkin(1) - pkin(5);
t109 = pkin(3) * t73;
t72 = sin(qJ(4));
t75 = cos(qJ(4));
t38 = Icges(5,6) * t73 + (Icges(5,4) * t75 - Icges(5,2) * t72) * t76;
t108 = t72 * t38;
t107 = t73 * t74;
t106 = t74 * t72;
t105 = t74 * t75;
t104 = t74 * t76;
t103 = t76 * t77;
t102 = t77 * t72;
t101 = t77 * t75;
t35 = Icges(5,3) * t73 + (Icges(5,5) * t75 - Icges(5,6) * t72) * t76;
t41 = Icges(5,5) * t73 + (Icges(5,1) * t75 - Icges(5,4) * t72) * t76;
t100 = t76 * t75 * t41 + t73 * t35;
t44 = t73 * rSges(5,3) + (rSges(5,1) * t75 - rSges(5,2) * t72) * t76;
t99 = t76 * pkin(3) + t73 * pkin(6) + t44;
t49 = -t73 * t106 + t101;
t50 = t73 * t105 + t102;
t98 = t50 * rSges(5,1) + t49 * rSges(5,2);
t97 = t77 * pkin(1) + t74 * qJ(2);
t93 = Icges(5,5) * t76;
t92 = Icges(5,6) * t76;
t91 = Icges(5,3) * t76;
t90 = rSges(4,1) * t107 + rSges(4,2) * t104 + t77 * rSges(4,3);
t89 = t77 * pkin(5) + t97;
t12 = -t35 * t104 + t49 * t38 + t50 * t41;
t21 = Icges(5,5) * t50 + Icges(5,6) * t49 - t74 * t91;
t23 = Icges(5,4) * t50 + Icges(5,2) * t49 - t74 * t92;
t25 = Icges(5,1) * t50 + Icges(5,4) * t49 - t74 * t93;
t9 = t73 * t21 + (-t23 * t72 + t25 * t75) * t76;
t88 = -t9 / 0.2e1 - t12 / 0.2e1;
t51 = t73 * t102 + t105;
t52 = -t73 * t101 + t106;
t22 = Icges(5,5) * t52 + Icges(5,6) * t51 + t77 * t91;
t24 = Icges(5,4) * t52 + Icges(5,2) * t51 + t77 * t92;
t26 = Icges(5,1) * t52 + Icges(5,4) * t51 + t77 * t93;
t10 = t73 * t22 + (-t24 * t72 + t26 * t75) * t76;
t13 = t35 * t103 + t51 * t38 + t52 * t41;
t87 = t10 / 0.2e1 + t13 / 0.2e1;
t86 = (-rSges(5,3) - pkin(6)) * t76;
t84 = -t52 * rSges(5,1) - t51 * rSges(5,2);
t79 = Icges(4,5) * t73 + Icges(4,6) * t76;
t66 = t77 * qJ(2);
t31 = t66 + t114 + (-rSges(4,3) + t110) * t74;
t32 = t89 + t90;
t78 = m(4) * (t74 * t31 - t77 * t32);
t64 = pkin(3) * t107;
t59 = t77 * rSges(2,1) - t74 * rSges(2,2);
t58 = t76 * rSges(4,1) - t73 * rSges(4,2);
t57 = -t74 * rSges(2,1) - t77 * rSges(2,2);
t46 = -t77 * rSges(3,2) + t74 * rSges(3,3) + t97;
t45 = t77 * rSges(3,3) + t66 + (rSges(3,2) - pkin(1)) * t74;
t36 = Icges(4,3) * t77 + t79 * t74;
t30 = t99 * t77;
t29 = t99 * t74;
t28 = rSges(5,3) * t103 - t84;
t27 = -rSges(5,3) * t104 + t98;
t20 = -t74 * t90 + (t74 * rSges(4,3) - t114) * t77;
t19 = t74 * t86 + t64 + t89 + t98;
t18 = t66 + t110 * t74 + (t86 + t109) * t77 + t84;
t17 = t44 * t103 - t73 * t28;
t16 = t44 * t104 + t73 * t27;
t15 = (-t76 * t108 + t100) * t73;
t14 = (-t27 * t77 - t28 * t74) * t76;
t11 = (pkin(6) * t104 - t27 - t64) * t74 + (t28 + (pkin(6) * t76 - t109) * t77) * t77;
t8 = t22 * t103 + t51 * t24 + t52 * t26;
t7 = t21 * t103 + t51 * t23 + t52 * t25;
t6 = -t22 * t104 + t49 * t24 + t50 * t26;
t5 = -t21 * t104 + t49 * t23 + t50 * t25;
t4 = t7 * t77 + t8 * t74;
t3 = t5 * t77 + t6 * t74;
t2 = t13 * t73 + (-t7 * t74 + t77 * t8) * t76;
t1 = t12 * t73 + (-t5 * t74 + t6 * t77) * t76;
t33 = [Icges(3,1) + Icges(2,3) + (Icges(4,1) * t76 - t108) * t76 + m(2) * (t57 ^ 2 + t59 ^ 2) + m(3) * (t45 ^ 2 + t46 ^ 2) + m(4) * (t31 ^ 2 + t32 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2) + t100 + (-0.2e1 * Icges(4,4) * t76 + Icges(4,2) * t73) * t73; m(3) * (t74 * t45 - t77 * t46) + t78 + m(5) * (t74 * t18 - t77 * t19); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1) * t96; m(5) * (t29 * t18 - t30 * t19) + t58 * t78 + (t71 / 0.2e1 + t70 / 0.2e1) * (-t116 + t117) + (t115 * t77 - t88) * t77 + (t115 * t74 + t87) * t74; m(5) * (t29 * t74 + t30 * t77) + m(4) * t96 * t58; m(4) * (t96 * t58 ^ 2 + t20 ^ 2) + m(5) * (t11 ^ 2 + t29 ^ 2 + t30 ^ 2) + (t71 * t36 + t3) * t77 + (t74 * t36 * t77 + t4 + t96 * (Icges(4,3) * t74 - t79 * t77)) * t74; t15 + m(5) * (t16 * t19 + t17 * t18) + (t88 * t74 + t87 * t77) * t76; m(5) * (-t16 * t77 + t17 * t74); m(5) * (t14 * t11 - t16 * t30 + t17 * t29) + t1 * t111 + t74 * t2 / 0.2e1 + t73 * (t10 * t74 + t9 * t77) / 0.2e1 + (-t74 * t3 / 0.2e1 + t4 * t111) * t76; m(5) * (t14 ^ 2 + t16 ^ 2 + t17 ^ 2) + t73 * t15 + (-t74 * t1 + t77 * t2 + t73 * (t10 * t77 - t74 * t9)) * t76;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t33(1), t33(2), t33(4), t33(7); t33(2), t33(3), t33(5), t33(8); t33(4), t33(5), t33(6), t33(9); t33(7), t33(8), t33(9), t33(10);];
Mq = res;
