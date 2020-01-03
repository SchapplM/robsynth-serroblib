% Calculate joint inertia matrix for
% S5RRPRP3
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
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:50:55
% EndTime: 2019-12-31 19:50:58
% DurationCPUTime: 0.75s
% Computational Cost: add. (1653->128), mult. (1294->188), div. (0->0), fcn. (1098->8), ass. (0->71)
t87 = pkin(8) + qJ(4);
t83 = cos(t87);
t147 = t83 ^ 2;
t146 = Icges(6,4) + Icges(5,5);
t145 = Icges(5,6) - Icges(6,6);
t82 = sin(t87);
t143 = t146 * t82;
t144 = t145 * t83;
t138 = t143 + t144;
t88 = qJ(1) + qJ(2);
t84 = sin(t88);
t80 = t84 ^ 2;
t85 = cos(t88);
t81 = t85 ^ 2;
t142 = -t145 * t82 + t146 * t83;
t140 = Icges(6,2) + Icges(5,3);
t132 = rSges(6,1) + pkin(4);
t137 = t132 * t83;
t136 = rSges(6,3) + qJ(5);
t135 = qJ(3) + rSges(4,3);
t89 = sin(pkin(8));
t90 = cos(pkin(8));
t133 = -rSges(4,1) * t90 + rSges(4,2) * t89 - pkin(2);
t131 = t140 * t85 - t142 * t84;
t130 = t140 * t84 + t142 * t85;
t53 = t82 * rSges(5,1) + t83 * rSges(5,2);
t127 = m(5) * t53;
t126 = m(6) * t82;
t92 = sin(qJ(1));
t125 = t92 * pkin(1);
t123 = rSges(5,1) * t83;
t121 = t82 * t85;
t120 = t83 * t85;
t91 = -pkin(7) - qJ(3);
t119 = t85 * t91;
t118 = -t132 * t82 + t136 * t83;
t117 = t84 * t82 * rSges(5,2) + t85 * rSges(5,3);
t116 = t80 + t81;
t111 = qJ(5) * t82;
t55 = t85 * rSges(3,1) - t84 * rSges(3,2);
t77 = t90 * pkin(3) + pkin(2);
t110 = t85 * t77 - t84 * t91;
t109 = t84 * rSges(6,2) + rSges(6,3) * t121 + t85 * t111 + t132 * t120;
t54 = -t84 * rSges(3,1) - t85 * rSges(3,2);
t96 = rSges(5,1) * t120 - rSges(5,2) * t121 + t84 * rSges(5,3);
t95 = t138 * t81 + (t138 / 0.2e1 + t144 / 0.2e1 + t143 / 0.2e1) * t80;
t24 = -t133 * t85 + t135 * t84;
t13 = t109 + t110;
t23 = t133 * t84 + t135 * t85;
t94 = Icges(4,2) * t90 ^ 2 + Icges(3,3) + (Icges(4,1) * t89 + 0.2e1 * Icges(4,4) * t90) * t89 + (Icges(5,2) + Icges(6,3)) * t147 + ((Icges(5,1) + Icges(6,1)) * t82 + (2 * Icges(5,4) - 2 * Icges(6,5)) * t83) * t82;
t18 = t110 + t96;
t17 = -t119 + (-t77 - t123) * t84 + t117;
t75 = t85 * rSges(6,2);
t12 = -t119 + t75 + (-t136 * t82 - t137 - t77) * t84;
t93 = cos(qJ(1));
t86 = t93 * pkin(1);
t65 = t93 * rSges(2,1) - t92 * rSges(2,2);
t64 = -t92 * rSges(2,1) - t93 * rSges(2,2);
t41 = t55 + t86;
t40 = t54 - t125;
t22 = t118 * t85;
t21 = t118 * t84;
t20 = t24 + t86;
t19 = t23 - t125;
t16 = t18 + t86;
t15 = t17 - t125;
t14 = t84 * (t123 * t84 - t117) + t85 * t96;
t11 = t13 + t86;
t10 = t12 - t125;
t1 = t109 * t85 + (-t75 + (rSges(6,3) * t82 + t111 + t137) * t84) * t84;
t2 = [Icges(2,3) + m(6) * (t10 ^ 2 + t11 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2) + m(4) * (t19 ^ 2 + t20 ^ 2) + m(3) * (t40 ^ 2 + t41 ^ 2) + m(2) * (t64 ^ 2 + t65 ^ 2) + t94; m(6) * (t12 * t10 + t13 * t11) + m(5) * (t17 * t15 + t18 * t16) + m(4) * (t23 * t19 + t24 * t20) + m(3) * (t54 * t40 + t55 * t41) + t94; m(6) * (t12 ^ 2 + t13 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2) + m(4) * (t23 ^ 2 + t24 ^ 2) + m(3) * (t54 ^ 2 + t55 ^ 2) + t94; m(6) * (t84 * t10 - t85 * t11) + m(5) * (t84 * t15 - t85 * t16) + m(4) * (t84 * t19 - t85 * t20); m(6) * (t84 * t12 - t85 * t13) + m(5) * (t84 * t17 - t85 * t18) + m(4) * (t84 * t23 - t85 * t24); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t116; m(6) * (t22 * t10 + t21 * t11) + (-t15 * t85 - t16 * t84) * t127 + t95; m(6) * (t22 * t12 + t21 * t13) + (-t17 * t85 - t18 * t84) * t127 + t95; m(6) * (-t21 * t85 + t22 * t84); m(5) * (t116 * t53 ^ 2 + t14 ^ 2) + m(6) * (t1 ^ 2 + t21 ^ 2 + t22 ^ 2) + t131 * t85 * t81 + (t130 * t80 + (t130 * t85 + t131 * t84) * t85) * t84; (t10 * t85 + t11 * t84) * t126; (t12 * t85 + t13 * t84) * t126; 0; m(6) * (-t83 * t1 + (t21 * t84 + t22 * t85) * t82); m(6) * (t116 * t82 ^ 2 + t147);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
