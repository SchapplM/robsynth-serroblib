% Calculate joint inertia matrix for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:22:00
% EndTime: 2019-12-05 18:22:02
% DurationCPUTime: 0.65s
% Computational Cost: add. (1446->119), mult. (1014->152), div. (0->0), fcn. (852->8), ass. (0->67)
t126 = Icges(5,5) + Icges(6,5);
t125 = Icges(5,6) + Icges(6,6);
t78 = cos(qJ(4));
t121 = t125 * t78;
t76 = sin(qJ(4));
t122 = t126 * t76;
t119 = t121 + t122;
t74 = qJ(1) + qJ(2);
t71 = pkin(8) + t74;
t68 = sin(t71);
t65 = t68 ^ 2;
t69 = cos(t71);
t66 = t69 ^ 2;
t124 = -t125 * t76 + t126 * t78;
t123 = Icges(5,3) + Icges(6,3);
t120 = rSges(6,1) + pkin(4);
t117 = t123 * t69 - t124 * t68;
t116 = t123 * t68 + t124 * t69;
t58 = t76 * rSges(5,1) + t78 * rSges(5,2);
t113 = m(5) * t58;
t72 = sin(t74);
t112 = pkin(2) * t72;
t73 = cos(t74);
t111 = pkin(2) * t73;
t77 = sin(qJ(1));
t110 = t77 * pkin(1);
t79 = cos(qJ(1));
t109 = t79 * pkin(1);
t108 = rSges(5,1) * t78;
t107 = t68 * t76;
t106 = t69 * t76;
t105 = -rSges(6,2) * t107 - t69 * rSges(6,3);
t104 = rSges(5,2) * t107 + t69 * rSges(5,3);
t75 = -qJ(5) - pkin(7);
t103 = -rSges(6,2) * t106 - t68 * t75;
t102 = t66 + t65;
t97 = -pkin(3) - t108;
t96 = t78 * rSges(6,2) + t120 * t76;
t40 = -t73 * rSges(3,1) + t72 * rSges(3,2);
t95 = -t120 * t78 - pkin(3);
t94 = -pkin(3) - t95;
t39 = -t72 * rSges(3,1) - t73 * rSges(3,2);
t36 = -t69 * rSges(4,1) + t68 * rSges(4,2) - t111;
t81 = Icges(3,3) + Icges(4,3) + (Icges(5,2) + Icges(6,2)) * t78 ^ 2 + ((Icges(5,1) + Icges(6,1)) * t76 + (2 * Icges(5,4) + 2 * Icges(6,4)) * t78) * t76;
t80 = t119 * t66 + (t121 / 0.2e1 + t122 / 0.2e1 + t119 / 0.2e1) * t65;
t35 = -t68 * rSges(4,1) - t69 * rSges(4,2) - t112;
t64 = t69 * pkin(7);
t17 = t97 * t68 + t104 - t112 + t64;
t13 = t95 * t68 - t69 * t75 - t105 - t112;
t14 = -t68 * rSges(6,3) + t95 * t69 - t103 - t111;
t48 = rSges(5,2) * t106;
t18 = -t111 + t48 + t97 * t69 + (-rSges(5,3) - pkin(7)) * t68;
t60 = -t79 * rSges(2,1) + t77 * rSges(2,2);
t59 = -t77 * rSges(2,1) - t79 * rSges(2,2);
t38 = t40 - t109;
t37 = t39 - t110;
t34 = t36 - t109;
t33 = t35 - t110;
t32 = t96 * t69;
t31 = t96 * t68;
t16 = t18 - t109;
t15 = t17 - t110;
t12 = t14 - t109;
t11 = t13 - t110;
t6 = t69 * (t68 * rSges(5,3) + t69 * t108 - t48) - t68 * (-t68 * t108 + t104);
t1 = (t94 * t69 + t103) * t69 + (t64 + t94 * t68 + (rSges(6,3) - pkin(7) + t75) * t69 + t105) * t68;
t2 = [Icges(2,3) + m(2) * (t59 ^ 2 + t60 ^ 2) + m(3) * (t37 ^ 2 + t38 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2) + m(6) * (t11 ^ 2 + t12 ^ 2) + t81; m(3) * (t39 * t37 + t40 * t38) + m(4) * (t35 * t33 + t36 * t34) + m(5) * (t17 * t15 + t18 * t16) + m(6) * (t13 * t11 + t14 * t12) + t81; m(6) * (t13 ^ 2 + t14 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2) + m(4) * (t35 ^ 2 + t36 ^ 2) + m(3) * (t39 ^ 2 + t40 ^ 2) + t81; 0; 0; m(4) + m(5) + m(6); m(6) * (-t32 * t11 + t31 * t12) + (-t15 * t69 + t16 * t68) * t113 + t80; m(6) * (-t32 * t13 + t31 * t14) + (-t17 * t69 + t18 * t68) * t113 + t80; m(5) * t6 + m(6) * t1; m(5) * (t102 * t58 ^ 2 + t6 ^ 2) + m(6) * (t1 ^ 2 + t31 ^ 2 + t32 ^ 2) + t116 * t68 * t65 + (t117 * t66 + (t116 * t69 + t117 * t68) * t68) * t69; m(6) * (t68 * t11 + t69 * t12); m(6) * (t68 * t13 + t69 * t14); 0; m(6) * (t69 * t31 - t68 * t32); m(6) * t102;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
