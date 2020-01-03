% Calculate joint inertia matrix for
% S5RPRRP9
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
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP9_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP9_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP9_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:48:37
% EndTime: 2019-12-31 18:48:41
% DurationCPUTime: 1.69s
% Computational Cost: add. (2188->180), mult. (1886->271), div. (0->0), fcn. (1681->8), ass. (0->99)
t182 = Icges(5,4) - Icges(6,5);
t181 = Icges(5,1) + Icges(6,1);
t180 = Icges(5,2) + Icges(6,3);
t96 = pkin(8) + qJ(3);
t90 = qJ(4) + t96;
t86 = cos(t90);
t179 = t182 * t86;
t85 = sin(t90);
t178 = t182 * t85;
t177 = Icges(6,4) + Icges(5,5);
t176 = Icges(5,6) - Icges(6,6);
t175 = -t180 * t85 + t179;
t174 = t181 * t86 - t178;
t102 = sin(qJ(1));
t103 = cos(qJ(1));
t173 = t175 * t102 - t176 * t103;
t172 = t176 * t102 + t175 * t103;
t171 = t174 * t102 - t177 * t103;
t170 = t177 * t102 + t174 * t103;
t169 = Icges(6,2) + Icges(5,3);
t168 = -t176 * t85 + t177 * t86;
t155 = rSges(6,3) + qJ(5);
t164 = rSges(6,1) + pkin(4);
t167 = t155 * t85 + t164 * t86;
t166 = t180 * t86 + t178;
t165 = t181 * t85 + t179;
t163 = -t168 * t102 + t103 * t169;
t162 = t102 * t169 + t168 * t103;
t160 = -t170 * t86 + t172 * t85;
t159 = -t171 * t86 + t173 * t85;
t97 = t102 ^ 2;
t157 = t176 * t86 + t177 * t85;
t156 = t102 * rSges(6,2);
t154 = t165 * t86 - t166 * t85;
t133 = t103 * t86;
t134 = t103 * t85;
t153 = t164 * t133 + t155 * t134 + t156;
t98 = t103 ^ 2;
t152 = t163 * t98 + (t160 * t102 + (-t159 + t162) * t103) * t102;
t151 = (t162 * t97 + ((-t160 + t163) * t102 + t159 * t103) * t103) * t102;
t88 = sin(t96);
t150 = pkin(3) * t88;
t149 = t102 / 0.2e1;
t148 = -t103 / 0.2e1;
t100 = cos(pkin(8));
t87 = t100 * pkin(2) + pkin(1);
t89 = cos(t96);
t147 = rSges(4,1) * t89;
t146 = rSges(4,2) * t88;
t101 = -pkin(6) - qJ(2);
t74 = pkin(3) * t89 + t87;
t68 = t103 * t74;
t145 = t103 * (-t103 * t87 + t68) + (t74 - t87) * t97;
t106 = rSges(5,1) * t133 - rSges(5,2) * t134 + t102 * rSges(5,3);
t125 = rSges(5,1) * t86 - rSges(5,2) * t85;
t18 = t102 * (-t103 * rSges(5,3) + t125 * t102) + t103 * t106;
t144 = t155 * t86 - t164 * t85;
t142 = t102 * rSges(4,3) + t103 * t147;
t141 = t97 + t98;
t140 = Icges(4,4) * t88;
t139 = Icges(4,4) * t89;
t132 = rSges(3,3) + qJ(2);
t67 = t85 * rSges(5,1) + t86 * rSges(5,2);
t129 = -t67 - t150;
t95 = -pkin(7) + t101;
t128 = -t102 * t95 + t68;
t7 = (t153 - t156) * t103 + t167 * t97;
t127 = t144 - t150;
t126 = -t146 + t147;
t116 = Icges(4,1) * t89 - t140;
t113 = -Icges(4,2) * t88 + t139;
t110 = Icges(4,5) * t89 - Icges(4,6) * t88;
t107 = t152 * t103 + t151;
t99 = sin(pkin(8));
t105 = rSges(3,1) * t100 - rSges(3,2) * t99 + pkin(1);
t104 = (t157 * t102 + t154 * t103 + t170 * t85 + t172 * t86) * t149 + (t154 * t102 - t157 * t103 + t171 * t85 + t173 * t86) * t148;
t81 = t103 * rSges(2,1) - t102 * rSges(2,2);
t80 = -t102 * rSges(2,1) - t103 * rSges(2,2);
t73 = t88 * rSges(4,1) + t89 * rSges(4,2);
t53 = Icges(4,3) * t102 + t110 * t103;
t52 = -Icges(4,3) * t103 + t110 * t102;
t37 = t132 * t102 + t105 * t103;
t36 = -t105 * t102 + t132 * t103;
t31 = t129 * t103;
t30 = t129 * t102;
t29 = -t102 * t101 + (t87 - t146) * t103 + t142;
t28 = (rSges(4,3) - t101) * t103 + (-t126 - t87) * t102;
t27 = t144 * t103;
t26 = t144 * t102;
t23 = t106 + t128;
t22 = (rSges(5,3) - t95) * t103 + (-t125 - t74) * t102;
t21 = t127 * t103;
t20 = t127 * t102;
t19 = t103 * (-t103 * t146 + t142) + (-t103 * rSges(4,3) + t126 * t102) * t102;
t17 = t128 + t153;
t16 = (rSges(6,2) - t95) * t103 + (-t167 - t74) * t102;
t6 = t18 + t145;
t5 = t7 + t145;
t1 = [Icges(3,1) * t99 ^ 2 + t89 * (Icges(4,2) * t89 + t140) + t88 * (Icges(4,1) * t88 + t139) + Icges(2,3) + (0.2e1 * Icges(3,4) * t99 + Icges(3,2) * t100) * t100 + t166 * t86 + t165 * t85 + m(5) * (t22 ^ 2 + t23 ^ 2) + m(6) * (t16 ^ 2 + t17 ^ 2) + m(4) * (t28 ^ 2 + t29 ^ 2) + m(3) * (t36 ^ 2 + t37 ^ 2) + m(2) * (t80 ^ 2 + t81 ^ 2); m(5) * (t102 * t22 - t103 * t23) + m(6) * (t102 * t16 - t103 * t17) + m(4) * (t102 * t28 - t103 * t29) + m(3) * (t102 * t36 - t103 * t37); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t141; (t89 * (Icges(4,6) * t102 + t113 * t103) + t88 * (Icges(4,5) * t102 + t116 * t103)) * t149 + (t89 * (-Icges(4,6) * t103 + t113 * t102) + t88 * (-Icges(4,5) * t103 + t116 * t102)) * t148 + m(5) * (t31 * t22 + t30 * t23) + m(6) * (t21 * t16 + t20 * t17) + m(4) * (-t102 * t29 - t103 * t28) * t73 + (t97 / 0.2e1 + t98 / 0.2e1) * (Icges(4,5) * t88 + Icges(4,6) * t89) + t104; m(5) * (t31 * t102 - t30 * t103) + m(6) * (t21 * t102 - t20 * t103); m(6) * (t20 ^ 2 + t21 ^ 2 + t5 ^ 2) + m(5) * (t30 ^ 2 + t31 ^ 2 + t6 ^ 2) + m(4) * (t141 * t73 ^ 2 + t19 ^ 2) + t102 * t97 * t53 + t151 + (-t98 * t52 + (-t102 * t52 + t103 * t53) * t102 + t152) * t103; m(6) * (t27 * t16 + t26 * t17) + m(5) * (-t102 * t23 - t103 * t22) * t67 + t104; m(6) * (t27 * t102 - t26 * t103); m(6) * (t26 * t20 + t27 * t21 + t7 * t5) + m(5) * (t18 * t6 + (-t102 * t30 - t103 * t31) * t67) + t107; m(5) * (t141 * t67 ^ 2 + t18 ^ 2) + m(6) * (t26 ^ 2 + t27 ^ 2 + t7 ^ 2) + t107; m(6) * (t102 * t17 + t103 * t16) * t85; 0; m(6) * (-t86 * t5 + (t102 * t20 + t103 * t21) * t85); m(6) * (-t86 * t7 + (t102 * t26 + t103 * t27) * t85); m(6) * (t141 * t85 ^ 2 + t86 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
