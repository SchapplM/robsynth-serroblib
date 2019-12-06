% Calculate joint inertia matrix for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:40
% EndTime: 2019-12-05 17:07:43
% DurationCPUTime: 0.54s
% Computational Cost: add. (2225->144), mult. (1314->206), div. (0->0), fcn. (1144->8), ass. (0->83)
t82 = pkin(9) + qJ(2);
t79 = qJ(3) + t82;
t74 = sin(t79);
t75 = cos(t79);
t123 = t74 * t75;
t71 = t74 ^ 2;
t72 = t75 ^ 2;
t122 = t74 / 0.2e1;
t121 = -t75 / 0.2e1;
t84 = sin(qJ(4));
t85 = cos(qJ(4));
t62 = t84 * rSges(5,1) + t85 * rSges(5,2);
t120 = m(5) * t62;
t83 = qJ(4) + qJ(5);
t80 = sin(t83);
t81 = cos(t83);
t50 = t80 * rSges(6,1) + t81 * rSges(6,2);
t119 = m(6) * t50;
t77 = sin(t82);
t118 = pkin(2) * t77;
t117 = rSges(5,1) * t85;
t116 = rSges(6,1) * t81;
t115 = rSges(5,2) * t84;
t114 = rSges(6,2) * t80;
t113 = t75 * rSges(6,3) + t74 * t114;
t88 = t74 * rSges(6,3) + (-t114 + t116) * t75;
t6 = t74 * (t74 * t116 - t113) + t75 * t88;
t112 = t75 * rSges(5,3) + t74 * t115;
t111 = -t75 * pkin(3) - t74 * pkin(7);
t110 = t71 + t72;
t109 = Icges(5,4) * t84;
t108 = Icges(5,4) * t85;
t107 = Icges(6,4) * t80;
t106 = Icges(6,4) * t81;
t47 = Icges(6,5) * t80 + Icges(6,6) * t81;
t92 = -Icges(6,2) * t80 + t106;
t94 = Icges(6,1) * t81 - t107;
t48 = Icges(6,2) * t81 + t107;
t49 = Icges(6,1) * t80 + t106;
t97 = -t48 * t80 + t49 * t81;
t105 = (t81 * (Icges(6,6) * t74 + t92 * t75) + t80 * (Icges(6,5) * t74 + t94 * t75) + t74 * t47 + t97 * t75) * t122 + (t81 * (-Icges(6,6) * t75 + t92 * t74) + t80 * (-Icges(6,5) * t75 + t94 * t74) - t75 * t47 + t97 * t74) * t121;
t90 = Icges(6,5) * t81 - Icges(6,6) * t80;
t24 = -Icges(6,3) * t75 + t90 * t74;
t25 = Icges(6,3) * t74 + t90 * t75;
t104 = -t75 * (-t25 * t123 + t72 * t24) + t74 * (-t24 * t123 + t71 * t25);
t103 = -pkin(4) * t84 - t50;
t43 = t75 * rSges(4,1) - t74 * rSges(4,2);
t60 = Icges(5,2) * t85 + t109;
t61 = Icges(5,1) * t84 + t108;
t102 = t81 * t48 + t80 * t49 + t85 * t60 + t84 * t61 + Icges(4,3);
t42 = -t74 * rSges(4,1) - t75 * rSges(4,2);
t96 = -t60 * t84 + t61 * t85;
t95 = Icges(5,1) * t85 - t109;
t93 = -Icges(5,2) * t84 + t108;
t91 = Icges(5,5) * t85 - Icges(5,6) * t84;
t89 = t74 * rSges(5,3) + (-t115 + t117) * t75;
t59 = Icges(5,5) * t84 + Icges(5,6) * t85;
t87 = t105 + (t85 * (Icges(5,6) * t74 + t93 * t75) + t84 * (Icges(5,5) * t74 + t95 * t75) + t74 * t59 + t96 * t75) * t122 + (t85 * (-Icges(5,6) * t75 + t93 * t74) + t84 * (-Icges(5,5) * t75 + t95 * t74) - t75 * t59 + t96 * t74) * t121;
t21 = t89 - t111;
t69 = t75 * pkin(7);
t20 = t69 + (-pkin(3) - t117) * t74 + t112;
t76 = t85 * pkin(4) + pkin(3);
t53 = t75 * t76;
t86 = -pkin(8) - pkin(7);
t17 = -t74 * t86 + t53 + t88;
t16 = -t75 * t86 + (-t76 - t116) * t74 + t113;
t78 = cos(t82);
t73 = pkin(2) * t78;
t45 = t78 * rSges(3,1) - t77 * rSges(3,2);
t44 = -t77 * rSges(3,1) - t78 * rSges(3,2);
t39 = t43 + t73;
t38 = t42 - t118;
t33 = Icges(5,3) * t74 + t91 * t75;
t32 = -Icges(5,3) * t75 + t91 * t74;
t31 = t103 * t75;
t30 = t103 * t74;
t19 = t21 + t73;
t18 = t20 - t118;
t15 = t17 + t73;
t14 = t16 - t118;
t11 = t74 * (t74 * t117 - t112) + t75 * t89;
t3 = t75 * (t53 + t111) + (t69 + (-pkin(3) + t76) * t74) * t74 + t6;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; Icges(3,3) + m(6) * (t14 ^ 2 + t15 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2) + m(4) * (t38 ^ 2 + t39 ^ 2) + m(3) * (t44 ^ 2 + t45 ^ 2) + t102; 0; m(6) * (t16 * t14 + t17 * t15) + m(5) * (t20 * t18 + t21 * t19) + m(4) * (t42 * t38 + t43 * t39) + t102; m(6) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t42 ^ 2 + t43 ^ 2) + t102; m(5) * t11 + m(6) * t3; m(6) * (t31 * t14 + t30 * t15) + (-t18 * t75 - t19 * t74) * t120 + t87; m(6) * (t31 * t16 + t30 * t17) + (-t20 * t75 - t21 * t74) * t120 + t87; m(5) * (t110 * t62 ^ 2 + t11 ^ 2) + t74 * (-t32 * t123 + t71 * t33) - t75 * (-t33 * t123 + t72 * t32) + m(6) * (t3 ^ 2 + t30 ^ 2 + t31 ^ 2) + t104; m(6) * t6; (-t14 * t75 - t15 * t74) * t119 + t105; (-t16 * t75 - t17 * t74) * t119 + t105; m(6) * (t6 * t3 + (-t30 * t74 - t31 * t75) * t50) + t104; m(6) * (t110 * t50 ^ 2 + t6 ^ 2) + t104;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
