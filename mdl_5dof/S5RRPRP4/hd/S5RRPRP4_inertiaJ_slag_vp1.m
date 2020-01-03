% Calculate joint inertia matrix for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:37
% EndTime: 2019-12-31 19:52:39
% DurationCPUTime: 0.68s
% Computational Cost: add. (1223->121), mult. (1260->167), div. (0->0), fcn. (1054->6), ass. (0->70)
t84 = cos(qJ(4));
t143 = t84 ^ 2;
t142 = Icges(6,4) + Icges(5,5);
t141 = Icges(5,6) - Icges(6,6);
t82 = sin(qJ(4));
t137 = t141 * t82;
t138 = t142 * t84;
t135 = t138 - t137;
t81 = qJ(1) + qJ(2);
t77 = sin(t81);
t75 = t77 ^ 2;
t78 = cos(t81);
t76 = t78 ^ 2;
t140 = t141 * t84 + t142 * t82;
t139 = Icges(6,2) + Icges(5,3);
t136 = rSges(6,1) + pkin(4);
t133 = rSges(6,3) + qJ(5);
t121 = t77 * t82;
t128 = t84 * t133;
t71 = t78 * rSges(6,2);
t132 = t136 * t121 - t128 * t77 + t71;
t120 = t78 * t82;
t131 = -t136 * t120 + t78 * t128;
t130 = t139 * t78 + t140 * t77;
t129 = t139 * t77 - t140 * t78;
t125 = -pkin(2) - pkin(7);
t124 = m(6) * t84;
t83 = sin(qJ(1));
t123 = t83 * pkin(1);
t122 = rSges(5,2) * t84;
t119 = t133 * t82 + t136 * t84;
t118 = rSges(5,1) * t120 + t78 * t122;
t116 = t78 * pkin(2) + t77 * qJ(3);
t115 = t75 + t76;
t108 = rSges(5,1) * t121 + t78 * rSges(5,3) + t77 * t122;
t107 = t78 * pkin(7) + t116;
t40 = t78 * rSges(3,1) - t77 * rSges(3,2);
t66 = t78 * qJ(3);
t8 = t66 + (-rSges(6,2) + t125) * t77 - t131;
t2 = t8 - t123;
t85 = cos(qJ(1));
t79 = t85 * pkin(1);
t9 = t107 + t132;
t3 = t79 + t9;
t104 = t77 * t2 - t78 * t3;
t103 = t77 * t8 - t78 * t9;
t39 = -t77 * rSges(3,1) - t78 * rSges(3,2);
t21 = t119 * t77;
t22 = t119 * t78;
t102 = t21 * t77 + t22 * t78;
t23 = t78 * rSges(4,3) + t66 + (rSges(4,2) - pkin(2)) * t77;
t24 = -t78 * rSges(4,2) + t77 * rSges(4,3) + t116;
t18 = t107 + t108;
t17 = t66 + (-rSges(5,3) + t125) * t77 + t118;
t15 = t17 - t123;
t16 = t79 + t18;
t89 = m(5) * (t77 * t15 - t78 * t16);
t88 = m(5) * (t77 * t17 - t78 * t18);
t87 = t135 * t76 + (t138 / 0.2e1 - t137 / 0.2e1 + t135 / 0.2e1) * t75;
t86 = Icges(4,1) + Icges(3,3) + (Icges(5,1) + Icges(6,1)) * t143 + ((Icges(5,2) + Icges(6,3)) * t82 + (-2 * Icges(5,4) + 2 * Icges(6,5)) * t84) * t82;
t56 = t85 * rSges(2,1) - t83 * rSges(2,2);
t55 = t84 * rSges(5,1) - t82 * rSges(5,2);
t52 = -t83 * rSges(2,1) - t85 * rSges(2,2);
t38 = t40 + t79;
t37 = t39 - t123;
t20 = t24 + t79;
t19 = t23 - t123;
t10 = t78 * (t77 * rSges(5,3) - t118) - t77 * t108;
t1 = t131 * t78 + (-t132 + t71) * t77;
t4 = [Icges(2,3) + m(6) * (t2 ^ 2 + t3 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2) + m(4) * (t19 ^ 2 + t20 ^ 2) + m(3) * (t37 ^ 2 + t38 ^ 2) + m(2) * (t52 ^ 2 + t56 ^ 2) + t86; m(6) * (t8 * t2 + t9 * t3) + m(5) * (t17 * t15 + t18 * t16) + m(4) * (t23 * t19 + t24 * t20) + m(3) * (t39 * t37 + t40 * t38) + t86; m(6) * (t8 ^ 2 + t9 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2) + m(4) * (t23 ^ 2 + t24 ^ 2) + m(3) * (t39 ^ 2 + t40 ^ 2) + t86; m(6) * t104 + t89 + m(4) * (t77 * t19 - t78 * t20); m(6) * t103 + t88 + m(4) * (t77 * t23 - t78 * t24); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t115; m(6) * (t21 * t2 - t22 * t3) + t55 * t89 + t87; m(6) * (t21 * t8 - t22 * t9) + t55 * t88 + t87; m(5) * t115 * t55 + m(6) * t102; m(5) * (t115 * t55 ^ 2 + t10 ^ 2) + m(6) * (t1 ^ 2 + t21 ^ 2 + t22 ^ 2) + t129 * t77 * t75 + (t130 * t76 + (t129 * t78 + t130 * t77) * t77) * t78; -t104 * t124; -t103 * t124; -t115 * t124; m(6) * (t82 * t1 - t102 * t84); m(6) * (t115 * t143 + t82 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;
