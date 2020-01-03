% Calculate joint inertia matrix for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR2_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR2_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR2_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:08
% EndTime: 2019-12-31 17:23:09
% DurationCPUTime: 0.52s
% Computational Cost: add. (1501->138), mult. (1252->203), div. (0->0), fcn. (1096->8), ass. (0->82)
t80 = qJ(1) + qJ(2);
t75 = sin(t80);
t77 = cos(t80);
t122 = t75 * t77;
t72 = t75 ^ 2;
t73 = t77 ^ 2;
t121 = t75 / 0.2e1;
t120 = -t77 / 0.2e1;
t81 = sin(qJ(3));
t83 = cos(qJ(3));
t58 = t81 * rSges(4,1) + t83 * rSges(4,2);
t119 = m(4) * t58;
t79 = qJ(3) + qJ(4);
t74 = sin(t79);
t76 = cos(t79);
t46 = t74 * rSges(5,1) + t76 * rSges(5,2);
t118 = m(5) * t46;
t82 = sin(qJ(1));
t117 = t82 * pkin(1);
t116 = rSges(4,1) * t83;
t115 = rSges(5,1) * t76;
t114 = rSges(4,2) * t81;
t113 = rSges(5,2) * t74;
t112 = t77 * rSges(5,3) + t75 * t113;
t87 = t75 * rSges(5,3) + (-t113 + t115) * t77;
t8 = t75 * (t75 * t115 - t112) + t77 * t87;
t111 = t77 * rSges(4,3) + t75 * t114;
t110 = -t77 * pkin(2) - t75 * pkin(6);
t109 = t72 + t73;
t108 = Icges(4,4) * t81;
t107 = Icges(4,4) * t83;
t106 = Icges(5,4) * t74;
t105 = Icges(5,4) * t76;
t43 = Icges(5,5) * t74 + Icges(5,6) * t76;
t91 = -Icges(5,2) * t74 + t105;
t93 = Icges(5,1) * t76 - t106;
t44 = Icges(5,2) * t76 + t106;
t45 = Icges(5,1) * t74 + t105;
t96 = -t44 * t74 + t45 * t76;
t104 = (t76 * (Icges(5,6) * t75 + t91 * t77) + t74 * (Icges(5,5) * t75 + t93 * t77) + t75 * t43 + t96 * t77) * t121 + (t76 * (-Icges(5,6) * t77 + t91 * t75) + t74 * (-Icges(5,5) * t77 + t93 * t75) - t77 * t43 + t96 * t75) * t120;
t89 = Icges(5,5) * t76 - Icges(5,6) * t74;
t24 = -Icges(5,3) * t77 + t89 * t75;
t25 = Icges(5,3) * t75 + t89 * t77;
t103 = -t77 * (-t25 * t122 + t73 * t24) + t75 * (-t24 * t122 + t72 * t25);
t102 = -pkin(3) * t81 - t46;
t48 = t77 * rSges(3,1) - t75 * rSges(3,2);
t56 = Icges(4,2) * t83 + t108;
t57 = Icges(4,1) * t81 + t107;
t101 = t76 * t44 + t74 * t45 + t83 * t56 + t81 * t57 + Icges(3,3);
t47 = -t75 * rSges(3,1) - t77 * rSges(3,2);
t95 = -t56 * t81 + t57 * t83;
t94 = Icges(4,1) * t83 - t108;
t92 = -Icges(4,2) * t81 + t107;
t90 = Icges(4,5) * t83 - Icges(4,6) * t81;
t88 = t75 * rSges(4,3) + (-t114 + t116) * t77;
t55 = Icges(4,5) * t81 + Icges(4,6) * t83;
t86 = t104 + (t83 * (Icges(4,6) * t75 + t92 * t77) + t81 * (Icges(4,5) * t75 + t94 * t77) + t75 * t55 + t95 * t77) * t121 + (t83 * (-Icges(4,6) * t77 + t92 * t75) + t81 * (-Icges(4,5) * t77 + t94 * t75) - t77 * t55 + t95 * t75) * t120;
t21 = t88 - t110;
t69 = t77 * pkin(6);
t20 = t69 + (-pkin(2) - t116) * t75 + t111;
t71 = t83 * pkin(3) + pkin(2);
t51 = t77 * t71;
t85 = -pkin(7) - pkin(6);
t17 = -t75 * t85 + t51 + t87;
t16 = -t77 * t85 + (-t71 - t115) * t75 + t112;
t84 = cos(qJ(1));
t78 = t84 * pkin(1);
t60 = t84 * rSges(2,1) - t82 * rSges(2,2);
t59 = -t82 * rSges(2,1) - t84 * rSges(2,2);
t41 = t48 + t78;
t40 = t47 - t117;
t33 = Icges(4,3) * t75 + t90 * t77;
t32 = -Icges(4,3) * t77 + t90 * t75;
t31 = t102 * t77;
t30 = t102 * t75;
t19 = t21 + t78;
t18 = t20 - t117;
t15 = t17 + t78;
t14 = t16 - t117;
t11 = t75 * (t75 * t116 - t111) + t77 * t88;
t3 = t77 * (t51 + t110) + (t69 + (-pkin(2) + t71) * t75) * t75 + t8;
t1 = [Icges(2,3) + m(2) * (t59 ^ 2 + t60 ^ 2) + m(3) * (t40 ^ 2 + t41 ^ 2) + m(4) * (t18 ^ 2 + t19 ^ 2) + m(5) * (t14 ^ 2 + t15 ^ 2) + t101; m(3) * (t47 * t40 + t48 * t41) + m(4) * (t20 * t18 + t21 * t19) + m(5) * (t16 * t14 + t17 * t15) + t101; m(5) * (t16 ^ 2 + t17 ^ 2) + m(4) * (t20 ^ 2 + t21 ^ 2) + m(3) * (t47 ^ 2 + t48 ^ 2) + t101; m(5) * (t31 * t14 + t30 * t15) + (-t18 * t77 - t19 * t75) * t119 + t86; m(5) * (t31 * t16 + t30 * t17) + (-t20 * t77 - t21 * t75) * t119 + t86; m(4) * (t109 * t58 ^ 2 + t11 ^ 2) + t75 * (-t122 * t32 + t72 * t33) - t77 * (-t122 * t33 + t73 * t32) + m(5) * (t3 ^ 2 + t30 ^ 2 + t31 ^ 2) + t103; (-t14 * t77 - t15 * t75) * t118 + t104; (-t16 * t77 - t17 * t75) * t118 + t104; m(5) * (t8 * t3 + (-t30 * t75 - t31 * t77) * t46) + t103; m(5) * (t109 * t46 ^ 2 + t8 ^ 2) + t103;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
