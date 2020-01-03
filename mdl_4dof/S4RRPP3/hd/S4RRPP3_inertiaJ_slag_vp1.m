% Calculate joint inertia matrix for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP3_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP3_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:23
% EndTime: 2019-12-31 16:57:25
% DurationCPUTime: 0.87s
% Computational Cost: add. (705->109), mult. (958->158), div. (0->0), fcn. (838->6), ass. (0->58)
t66 = qJ(2) + pkin(6);
t61 = cos(t66);
t132 = t61 ^ 2;
t71 = sin(qJ(1));
t67 = t71 ^ 2;
t73 = cos(qJ(1));
t68 = t73 ^ 2;
t104 = t67 + t68;
t130 = Icges(5,4) + Icges(4,5);
t129 = Icges(4,6) - Icges(5,6);
t60 = sin(t66);
t70 = sin(qJ(2));
t72 = cos(qJ(2));
t125 = Icges(3,5) * t72 - Icges(3,6) * t70 - t129 * t60 + t130 * t61;
t124 = Icges(5,2) + Icges(3,3) + Icges(4,3);
t120 = t71 * pkin(5);
t117 = rSges(5,1) + pkin(3);
t119 = t117 * t61;
t118 = rSges(5,3) + qJ(4);
t116 = t124 * t73 - t125 * t71;
t115 = t124 * t71 + t125 * t73;
t112 = pkin(2) * t70;
t111 = rSges(3,1) * t72;
t110 = rSges(3,2) * t70;
t109 = t60 * t73;
t108 = t61 * t73;
t107 = t73 * rSges(3,3);
t58 = pkin(2) * t72 + pkin(1);
t55 = t73 * t58;
t65 = t73 * pkin(5);
t106 = t71 * (t65 + (-pkin(1) + t58) * t71) + t73 * (-t73 * pkin(1) - t120 + t55);
t105 = rSges(3,3) * t71 + t111 * t73;
t97 = qJ(4) * t60;
t96 = -rSges(4,1) * t60 - rSges(4,2) * t61 - t112;
t69 = -qJ(3) - pkin(5);
t95 = -t71 * t69 + t55;
t94 = -t117 * t60 + t118 * t61 - t112;
t92 = t71 * rSges(5,2) + rSges(5,3) * t109 + t108 * t117 + t73 * t97;
t91 = -t110 + t111;
t90 = rSges(4,1) * t61 - rSges(4,2) * t60;
t74 = rSges(4,1) * t108 - rSges(4,2) * t109 + rSges(4,3) * t71;
t49 = rSges(2,1) * t73 - rSges(2,2) * t71;
t48 = -rSges(2,1) * t71 - rSges(2,2) * t73;
t47 = rSges(3,1) * t70 + rSges(3,2) * t72;
t15 = t96 * t73;
t14 = t96 * t71;
t13 = t120 + (pkin(1) - t110) * t73 + t105;
t12 = t107 + t65 + (-pkin(1) - t91) * t71;
t9 = t74 + t95;
t8 = (rSges(4,3) - t69) * t73 + (-t58 - t90) * t71;
t7 = t94 * t73;
t6 = t94 * t71;
t5 = t73 * (-t110 * t73 + t105) + (t71 * t91 - t107) * t71;
t4 = t92 + t95;
t3 = (rSges(5,2) - t69) * t73 + (-t118 * t60 - t119 - t58) * t71;
t2 = t73 * t74 + (-t73 * rSges(4,3) + t71 * t90) * t71 + t106;
t1 = t92 * t73 + (-t73 * rSges(5,2) + (rSges(5,3) * t60 + t119 + t97) * t71) * t71 + t106;
t10 = [Icges(3,2) * t72 ^ 2 + Icges(2,3) + m(2) * (t48 ^ 2 + t49 ^ 2) + m(3) * (t12 ^ 2 + t13 ^ 2) + m(4) * (t8 ^ 2 + t9 ^ 2) + m(5) * (t3 ^ 2 + t4 ^ 2) + (Icges(4,2) + Icges(5,3)) * t132 + (Icges(3,1) * t70 + 0.2e1 * Icges(3,4) * t72) * t70 + (0.2e1 * (Icges(4,4) - Icges(5,5)) * t61 + (Icges(4,1) + Icges(5,1)) * t60) * t60; m(4) * (t14 * t9 + t15 * t8) + m(5) * (t3 * t7 + t4 * t6) + m(3) * (-t12 * t73 - t13 * t71) * t47 + t104 * (Icges(3,5) * t70 + Icges(3,6) * t72 + t129 * t61 + t130 * t60); m(5) * (t1 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(4) * (t14 ^ 2 + t15 ^ 2 + t2 ^ 2) + m(3) * (t104 * t47 ^ 2 + t5 ^ 2) + t115 * t71 * t67 + (t116 * t68 + (t115 * t73 + t116 * t71) * t71) * t73; m(4) * (t71 * t8 - t73 * t9) + m(5) * (t3 * t71 - t4 * t73); m(5) * (-t6 * t73 + t7 * t71) + m(4) * (-t14 * t73 + t15 * t71); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * t104; m(5) * (t3 * t73 + t4 * t71) * t60; m(5) * (-t61 * t1 + (t6 * t71 + t7 * t73) * t60); 0; m(5) * (t104 * t60 ^ 2 + t132);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t10(1), t10(2), t10(4), t10(7); t10(2), t10(3), t10(5), t10(8); t10(4), t10(5), t10(6), t10(9); t10(7), t10(8), t10(9), t10(10);];
Mq = res;
