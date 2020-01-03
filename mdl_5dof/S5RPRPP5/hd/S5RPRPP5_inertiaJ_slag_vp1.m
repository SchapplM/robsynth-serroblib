% Calculate joint inertia matrix for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:08
% EndTime: 2019-12-31 18:16:10
% DurationCPUTime: 0.94s
% Computational Cost: add. (563->145), mult. (1248->201), div. (0->0), fcn. (1064->4), ass. (0->81)
t153 = Icges(6,4) + Icges(5,5);
t152 = Icges(4,5) + Icges(5,4);
t151 = -Icges(6,2) - Icges(5,3);
t150 = Icges(4,6) + Icges(6,6);
t149 = Icges(4,1) + Icges(5,1) + Icges(6,1);
t73 = sin(qJ(3));
t148 = t153 * t73;
t75 = cos(qJ(3));
t147 = (-Icges(4,4) + t153) * t75;
t146 = t151 * t75 + t148;
t145 = (Icges(5,6) - t150) * t75 + (Icges(6,5) - t152) * t73;
t144 = t149 * t73 - t147;
t143 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t134 = Icges(5,6) / 0.2e1;
t142 = t134 - Icges(4,6) / 0.2e1 - Icges(6,6) / 0.2e1;
t136 = -Icges(6,5) / 0.2e1;
t141 = t136 + Icges(4,5) / 0.2e1 + Icges(5,4) / 0.2e1;
t74 = sin(qJ(1));
t140 = t74 / 0.2e1;
t76 = cos(qJ(1));
t139 = -t76 / 0.2e1;
t132 = rSges(6,1) + pkin(4);
t100 = (-rSges(5,3) - qJ(4)) * t75;
t63 = t76 * t73 * pkin(3);
t65 = t76 * qJ(2);
t115 = t63 + t65;
t120 = rSges(5,1) * t73;
t122 = -pkin(1) - pkin(6);
t5 = (t100 + t120) * t76 + (-rSges(5,2) + t122) * t74 + t115;
t119 = t73 * t74;
t67 = t76 * rSges(5,2);
t116 = rSges(5,1) * t119 + t67;
t114 = t76 * pkin(1) + t74 * qJ(2);
t103 = t76 * pkin(6) + t114;
t62 = pkin(3) * t119;
t99 = t62 + t103;
t6 = t100 * t74 + t116 + t99;
t129 = m(5) * (t74 * t5 - t76 * t6);
t101 = (-rSges(6,2) - qJ(4)) * t75;
t106 = -rSges(6,3) - qJ(5);
t126 = t132 * t73;
t3 = (t101 + t126) * t76 + (-t106 + t122) * t74 + t115;
t117 = t132 * t119;
t4 = t74 * t101 + t106 * t76 + t117 + t99;
t128 = m(6) * (t74 * t3 - t76 * t4);
t127 = (rSges(4,1) * t73 + rSges(4,2) * t75) * t76;
t70 = t74 ^ 2;
t72 = t76 ^ 2;
t56 = t70 + t72;
t125 = t143 * t76 - t145 * t74;
t124 = t143 * t74 + t145 * t76;
t123 = m(5) / 0.2e1;
t118 = t74 * t75;
t107 = qJ(4) * t75;
t105 = t123 + m(6) / 0.2e1;
t104 = rSges(4,1) * t119 + rSges(4,2) * t118 + t76 * rSges(4,3);
t102 = t73 * rSges(6,2) + t132 * t75;
t51 = t75 * pkin(3) + t73 * qJ(4);
t37 = t74 * t51;
t8 = t102 * t74 + t37;
t9 = (-t102 - t51) * t76;
t95 = t8 * t74 - t9 * t76;
t53 = t75 * rSges(5,1) + t73 * rSges(5,3);
t12 = t74 * t53 + t37;
t13 = (-t51 - t53) * t76;
t93 = t12 * t74 - t13 * t76;
t10 = t65 + t127 + (-rSges(4,3) + t122) * t74;
t11 = t103 + t104;
t77 = m(4) * (t74 * t10 - t76 * t11);
t55 = t76 * rSges(2,1) - t74 * rSges(2,2);
t54 = t75 * rSges(4,1) - t73 * rSges(4,2);
t50 = -t74 * rSges(2,1) - t76 * rSges(2,2);
t39 = m(6) * t56;
t36 = -t74 * t107 + t62;
t34 = t76 * (t76 * t107 - t63);
t33 = -t76 * rSges(3,2) + t74 * rSges(3,3) + t114;
t32 = t76 * rSges(3,3) + t65 + (rSges(3,2) - pkin(1)) * t74;
t7 = -t74 * t104 + (t74 * rSges(4,3) - t127) * t76;
t2 = t34 + (rSges(5,3) * t75 - t120) * t72 + (rSges(5,3) * t118 - t116 - t36 + t67) * t74;
t1 = t34 + (rSges(6,2) * t118 - t117 - t36) * t74 + (rSges(6,2) * t75 - t126) * t72;
t14 = [Icges(3,1) + Icges(2,3) + m(2) * (t50 ^ 2 + t55 ^ 2) + m(3) * (t32 ^ 2 + t33 ^ 2) + m(4) * (t10 ^ 2 + t11 ^ 2) + m(5) * (t5 ^ 2 + t6 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2) + (-Icges(4,4) * t73 + t149 * t75 + t148) * t75 + ((Icges(4,2) - t151) * t73 + t147) * t73; m(3) * (t74 * t32 - t76 * t33) + t77 + t129 + t128; t39 + 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t123) * t56; m(5) * (t12 * t5 + t13 * t6) + m(6) * (t8 * t3 + t9 * t4) + t54 * t77 + ((t144 * t140 + t141 * t76) * t76 + (t144 * t139 + t141 * t74) * t74) * t75 + ((t146 * t140 + t142 * t76) * t76 + (t146 * t139 + t142 * t74) * t74) * t73 + t56 * ((t136 + t152 / 0.2e1) * t75 + (t134 - t150 / 0.2e1) * t73); m(4) * t56 * t54 + m(5) * t93 + m(6) * t95; m(6) * (t1 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(5) * (t12 ^ 2 + t13 ^ 2 + t2 ^ 2) + m(4) * (t56 * t54 ^ 2 + t7 ^ 2) + t124 * t74 * t70 + (t125 * t72 + (t124 * t76 + t125 * t74) * t74) * t76; 0.2e1 * (-t129 / 0.2e1 - t128 / 0.2e1) * t75; -0.2e1 * t105 * t56 * t75; m(6) * (t73 * t1 - t95 * t75) + m(5) * (t73 * t2 - t93 * t75); 0.2e1 * t105 * (t56 * t75 ^ 2 + t73 ^ 2); m(6) * (-t76 * t3 - t74 * t4); 0; m(6) * (-t74 * t9 - t76 * t8); 0; t39;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t14(1), t14(2), t14(4), t14(7), t14(11); t14(2), t14(3), t14(5), t14(8), t14(12); t14(4), t14(5), t14(6), t14(9), t14(13); t14(7), t14(8), t14(9), t14(10), t14(14); t14(11), t14(12), t14(13), t14(14), t14(15);];
Mq = res;
