% Calculate joint inertia matrix for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:39:48
% EndTime: 2019-12-05 16:39:51
% DurationCPUTime: 0.57s
% Computational Cost: add. (1361->102), mult. (940->138), div. (0->0), fcn. (796->6), ass. (0->59)
t123 = Icges(5,5) + Icges(6,5);
t122 = Icges(5,6) + Icges(6,6);
t76 = sin(qJ(4));
t120 = t123 * t76;
t77 = cos(qJ(4));
t121 = t122 * t77;
t115 = t120 + t121;
t74 = pkin(8) + qJ(2);
t73 = qJ(3) + t74;
t68 = sin(t73);
t65 = t68 ^ 2;
t69 = cos(t73);
t66 = t69 ^ 2;
t119 = -t122 * t76 + t123 * t77;
t116 = Icges(5,3) + Icges(6,3);
t113 = t116 * t69 - t119 * t68;
t112 = t116 * t68 + t119 * t69;
t56 = t76 * rSges(5,1) + t77 * rSges(5,2);
t109 = m(5) * t56;
t71 = sin(t74);
t108 = pkin(2) * t71;
t107 = rSges(5,1) * t77;
t106 = t68 * t76;
t105 = t69 * t76;
t104 = t69 * t77;
t103 = -rSges(6,2) * t106 - t69 * rSges(6,3);
t102 = rSges(5,2) * t106 + t69 * rSges(5,3);
t101 = t69 * pkin(3) + t68 * pkin(7);
t100 = t65 + t66;
t95 = -t77 * rSges(6,2) + (-rSges(6,1) - pkin(4)) * t76;
t70 = t77 * pkin(4) + pkin(3);
t94 = -rSges(6,1) * t77 - t70;
t36 = t69 * rSges(4,1) - t68 * rSges(4,2);
t93 = Icges(4,3) + (Icges(5,2) + Icges(6,2)) * t77 ^ 2 + ((Icges(5,1) + Icges(6,1)) * t76 + (2 * Icges(5,4) + 2 * Icges(6,4)) * t77) * t76;
t35 = -t68 * rSges(4,1) - t69 * rSges(4,2);
t80 = rSges(5,1) * t104 - rSges(5,2) * t105 + t68 * rSges(5,3);
t79 = rSges(6,1) * t104 - rSges(6,2) * t105 + t68 * rSges(6,3) + t69 * t70;
t78 = t115 * t66 + (t121 / 0.2e1 + t120 / 0.2e1 + t115 / 0.2e1) * t65;
t18 = t80 + t101;
t63 = t69 * pkin(7);
t17 = t63 + (-pkin(3) - t107) * t68 + t102;
t75 = -qJ(5) - pkin(7);
t14 = -t68 * t75 + t79;
t13 = t94 * t68 - t69 * t75 - t103;
t72 = cos(t74);
t67 = pkin(2) * t72;
t38 = t72 * rSges(3,1) - t71 * rSges(3,2);
t37 = -t71 * rSges(3,1) - t72 * rSges(3,2);
t34 = t36 + t67;
t33 = t35 - t108;
t32 = t95 * t69;
t31 = t95 * t68;
t16 = t18 + t67;
t15 = t17 - t108;
t12 = t14 + t67;
t11 = t13 - t108;
t6 = t68 * (t68 * t107 - t102) + t69 * t80;
t1 = (t79 - t101) * t69 + (t63 + (-pkin(3) - t94) * t68 + t103) * t68;
t2 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; Icges(3,3) + m(6) * (t11 ^ 2 + t12 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2) + m(3) * (t37 ^ 2 + t38 ^ 2) + t93; 0; m(6) * (t13 * t11 + t14 * t12) + m(5) * (t17 * t15 + t18 * t16) + m(4) * (t35 * t33 + t36 * t34) + t93; m(6) * (t13 ^ 2 + t14 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2) + m(4) * (t35 ^ 2 + t36 ^ 2) + t93; m(5) * t6 + m(6) * t1; m(6) * (t32 * t11 + t31 * t12) + (-t15 * t69 - t16 * t68) * t109 + t78; m(6) * (t32 * t13 + t31 * t14) + (-t17 * t69 - t18 * t68) * t109 + t78; m(5) * (t100 * t56 ^ 2 + t6 ^ 2) + m(6) * (t1 ^ 2 + t31 ^ 2 + t32 ^ 2) + t112 * t68 * t65 + (t113 * t66 + (t112 * t69 + t113 * t68) * t68) * t69; 0; m(6) * (t68 * t11 - t69 * t12); m(6) * (t68 * t13 - t69 * t14); m(6) * (-t69 * t31 + t68 * t32); m(6) * t100;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
