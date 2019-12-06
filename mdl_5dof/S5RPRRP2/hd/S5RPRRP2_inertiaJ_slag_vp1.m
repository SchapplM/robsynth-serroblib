% Calculate joint inertia matrix for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:01:23
% EndTime: 2019-12-05 18:01:26
% DurationCPUTime: 0.62s
% Computational Cost: add. (1385->112), mult. (970->146), div. (0->0), fcn. (820->8), ass. (0->65)
t123 = Icges(5,5) + Icges(6,5);
t122 = Icges(5,6) + Icges(6,6);
t75 = cos(qJ(4));
t118 = t122 * t75;
t73 = sin(qJ(4));
t119 = t123 * t73;
t116 = t118 + t119;
t71 = qJ(1) + pkin(8);
t70 = qJ(3) + t71;
t65 = sin(t70);
t63 = t65 ^ 2;
t66 = cos(t70);
t64 = t66 ^ 2;
t121 = -t122 * t73 + t123 * t75;
t120 = Icges(5,3) + Icges(6,3);
t117 = rSges(6,1) + pkin(4);
t114 = t120 * t66 - t121 * t65;
t113 = t120 * t65 + t121 * t66;
t56 = t73 * rSges(5,1) + t75 * rSges(5,2);
t110 = m(5) * t56;
t74 = sin(qJ(1));
t109 = t74 * pkin(1);
t76 = cos(qJ(1));
t108 = t76 * pkin(1);
t107 = rSges(5,1) * t75;
t106 = t65 * t73;
t105 = t66 * t73;
t104 = rSges(6,2) * t106 + t66 * rSges(6,3);
t103 = rSges(5,2) * t106 + t66 * rSges(5,3);
t72 = -qJ(5) - pkin(7);
t102 = -rSges(6,2) * t105 - t65 * t72;
t101 = t64 + t63;
t96 = -pkin(3) - t107;
t95 = t75 * rSges(6,2) + t117 * t73;
t38 = -t66 * rSges(4,1) + t65 * rSges(4,2);
t94 = -t117 * t75 - pkin(3);
t93 = -pkin(3) - t94;
t68 = sin(t71);
t92 = -pkin(2) * t68 - t109;
t69 = cos(t71);
t91 = -pkin(2) * t69 - t108;
t90 = Icges(4,3) + (Icges(5,2) + Icges(6,2)) * t75 ^ 2 + ((Icges(5,1) + Icges(6,1)) * t73 + (2 * Icges(5,4) + 2 * Icges(6,4)) * t75) * t73;
t37 = -t65 * rSges(4,1) - t66 * rSges(4,2);
t77 = t116 * t64 + (t118 / 0.2e1 + t119 / 0.2e1 + t116 / 0.2e1) * t63;
t62 = t66 * pkin(7);
t17 = t96 * t65 + t103 + t62;
t16 = -t65 * rSges(6,3) + t94 * t66 - t102;
t15 = t94 * t65 - t66 * t72 + t104;
t46 = rSges(5,2) * t105;
t18 = t46 + t96 * t66 + (-rSges(5,3) - pkin(7)) * t65;
t58 = -t76 * rSges(2,1) + t74 * rSges(2,2);
t57 = -t74 * rSges(2,1) - t76 * rSges(2,2);
t36 = -t69 * rSges(3,1) + t68 * rSges(3,2) - t108;
t35 = -t68 * rSges(3,1) - t69 * rSges(3,2) - t109;
t34 = t38 + t91;
t33 = t37 + t92;
t32 = t95 * t66;
t31 = t95 * t65;
t14 = t18 + t91;
t13 = t17 + t92;
t12 = t16 + t91;
t11 = t15 + t92;
t6 = t66 * (t65 * rSges(5,3) + t66 * t107 - t46) - t65 * (-t65 * t107 + t103);
t1 = (t93 * t66 + t102) * t66 + (t62 + t93 * t65 + (rSges(6,3) - pkin(7) + t72) * t66 - t104) * t65;
t2 = [Icges(2,3) + Icges(3,3) + m(2) * (t57 ^ 2 + t58 ^ 2) + m(3) * (t35 ^ 2 + t36 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2) + m(5) * (t13 ^ 2 + t14 ^ 2) + m(6) * (t11 ^ 2 + t12 ^ 2) + t90; 0; m(3) + m(4) + m(5) + m(6); m(4) * (t37 * t33 + t38 * t34) + m(5) * (t17 * t13 + t18 * t14) + m(6) * (t15 * t11 + t16 * t12) + t90; 0; m(6) * (t15 ^ 2 + t16 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2) + m(4) * (t37 ^ 2 + t38 ^ 2) + t90; m(6) * (-t32 * t11 + t31 * t12) + (-t13 * t66 + t14 * t65) * t110 + t77; m(5) * t6 + m(6) * t1; m(6) * (-t32 * t15 + t31 * t16) + (-t17 * t66 + t18 * t65) * t110 + t77; m(5) * (t101 * t56 ^ 2 + t6 ^ 2) + m(6) * (t1 ^ 2 + t31 ^ 2 + t32 ^ 2) + t113 * t65 * t63 + (t114 * t64 + (t113 * t66 + t114 * t65) * t65) * t66; m(6) * (t65 * t11 + t66 * t12); 0; m(6) * (t65 * t15 + t66 * t16); m(6) * (t66 * t31 - t65 * t32); m(6) * t101;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
