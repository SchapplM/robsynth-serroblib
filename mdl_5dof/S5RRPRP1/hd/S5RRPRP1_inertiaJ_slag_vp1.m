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
% m [6x1]
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
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:19:32
% EndTime: 2022-01-20 10:19:33
% DurationCPUTime: 0.62s
% Computational Cost: add. (1446->116), mult. (1014->153), div. (0->0), fcn. (852->8), ass. (0->67)
t132 = Icges(5,5) + Icges(6,5);
t131 = Icges(5,6) + Icges(6,6);
t82 = sin(qJ(4));
t129 = t132 * t82;
t84 = cos(qJ(4));
t130 = t131 * t84;
t124 = t129 + t130;
t80 = qJ(1) + qJ(2);
t76 = pkin(8) + t80;
t72 = sin(t76);
t69 = t72 ^ 2;
t73 = cos(t76);
t70 = t73 ^ 2;
t128 = -t131 * t82 + t132 * t84;
t125 = Icges(5,3) + Icges(6,3);
t122 = t125 * t73 - t128 * t72;
t121 = t125 * t72 + t128 * t73;
t58 = t82 * rSges(5,1) + t84 * rSges(5,2);
t118 = m(5) * t58;
t77 = sin(t80);
t117 = pkin(2) * t77;
t83 = sin(qJ(1));
t116 = t83 * pkin(1);
t115 = rSges(5,1) * t84;
t114 = t72 * t82;
t113 = t73 * t82;
t112 = t73 * t84;
t111 = -rSges(6,2) * t114 - t73 * rSges(6,3);
t110 = rSges(5,2) * t114 + t73 * rSges(5,3);
t109 = t73 * pkin(3) + t72 * pkin(7);
t108 = t69 + t70;
t103 = -t84 * rSges(6,2) + (-rSges(6,1) - pkin(4)) * t82;
t75 = t84 * pkin(4) + pkin(3);
t102 = -rSges(6,1) * t84 - t75;
t78 = cos(t80);
t40 = t78 * rSges(3,1) - t77 * rSges(3,2);
t74 = pkin(2) * t78;
t36 = t73 * rSges(4,1) - t72 * rSges(4,2) + t74;
t39 = -t77 * rSges(3,1) - t78 * rSges(3,2);
t89 = rSges(5,1) * t112 - rSges(5,2) * t113 + t72 * rSges(5,3);
t88 = Icges(3,3) + Icges(4,3) + (Icges(5,2) + Icges(6,2)) * t84 ^ 2 + ((Icges(5,1) + Icges(6,1)) * t82 + (2 * Icges(5,4) + 2 * Icges(6,4)) * t84) * t82;
t87 = rSges(6,1) * t112 - rSges(6,2) * t113 + t72 * rSges(6,3) + t73 * t75;
t86 = t124 * t70 + (t130 / 0.2e1 + t129 / 0.2e1 + t124 / 0.2e1) * t69;
t35 = -t72 * rSges(4,1) - t73 * rSges(4,2) - t117;
t18 = t74 + t89 + t109;
t81 = -qJ(5) - pkin(7);
t14 = -t72 * t81 + t74 + t87;
t67 = t73 * pkin(7);
t17 = -t117 + t67 + (-pkin(3) - t115) * t72 + t110;
t13 = t102 * t72 - t73 * t81 - t111 - t117;
t85 = cos(qJ(1));
t79 = t85 * pkin(1);
t60 = t85 * rSges(2,1) - t83 * rSges(2,2);
t59 = -t83 * rSges(2,1) - t85 * rSges(2,2);
t38 = t40 + t79;
t37 = t39 - t116;
t34 = t36 + t79;
t33 = t35 - t116;
t32 = t103 * t73;
t31 = t103 * t72;
t16 = t18 + t79;
t15 = t17 - t116;
t12 = t14 + t79;
t11 = t13 - t116;
t6 = t72 * (t72 * t115 - t110) + t73 * t89;
t1 = (t87 - t109) * t73 + (t67 + (-pkin(3) - t102) * t72 + t111) * t72;
t2 = [Icges(2,3) + m(6) * (t11 ^ 2 + t12 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2) + m(3) * (t37 ^ 2 + t38 ^ 2) + m(2) * (t59 ^ 2 + t60 ^ 2) + t88; m(6) * (t13 * t11 + t14 * t12) + m(5) * (t17 * t15 + t18 * t16) + m(4) * (t35 * t33 + t36 * t34) + m(3) * (t39 * t37 + t40 * t38) + t88; m(6) * (t13 ^ 2 + t14 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2) + m(4) * (t35 ^ 2 + t36 ^ 2) + m(3) * (t39 ^ 2 + t40 ^ 2) + t88; 0; 0; m(4) + m(5) + m(6); m(6) * (t32 * t11 + t31 * t12) + (-t15 * t73 - t16 * t72) * t118 + t86; m(6) * (t32 * t13 + t31 * t14) + (-t17 * t73 - t18 * t72) * t118 + t86; m(5) * t6 + m(6) * t1; m(5) * (t108 * t58 ^ 2 + t6 ^ 2) + m(6) * (t1 ^ 2 + t31 ^ 2 + t32 ^ 2) + t121 * t72 * t69 + (t122 * t70 + (t121 * t73 + t122 * t72) * t72) * t73; m(6) * (t72 * t11 - t73 * t12); m(6) * (t72 * t13 - t73 * t14); 0; m(6) * (-t73 * t31 + t72 * t32); m(6) * t108;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
