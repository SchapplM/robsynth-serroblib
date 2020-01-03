% Calculate joint inertia matrix for
% S5RPRRP4
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
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:49:15
% EndTime: 2020-01-03 11:49:24
% DurationCPUTime: 3.11s
% Computational Cost: add. (4343->322), mult. (6023->450), div. (0->0), fcn. (6406->8), ass. (0->151)
t149 = qJ(3) + qJ(4);
t139 = sin(t149);
t140 = cos(t149);
t150 = sin(pkin(8));
t151 = cos(pkin(8));
t102 = -Icges(6,5) * t151 + (Icges(6,1) * t140 - Icges(6,4) * t139) * t150;
t103 = -Icges(5,5) * t151 + (Icges(5,1) * t140 - Icges(5,4) * t139) * t150;
t215 = t102 + t103;
t155 = cos(qJ(1));
t153 = sin(qJ(1));
t182 = t151 * t153;
t114 = -t139 * t182 - t140 * t155;
t115 = -t139 * t155 + t140 * t182;
t176 = t153 * t150;
t62 = Icges(6,5) * t115 + Icges(6,6) * t114 + Icges(6,3) * t176;
t64 = Icges(5,5) * t115 + Icges(5,6) * t114 + Icges(5,3) * t176;
t214 = t62 + t64;
t181 = t151 * t155;
t116 = t139 * t181 - t140 * t153;
t117 = -t139 * t153 - t140 * t181;
t172 = t155 * t150;
t63 = Icges(6,5) * t117 + Icges(6,6) * t116 - Icges(6,3) * t172;
t65 = Icges(5,5) * t117 + Icges(5,6) * t116 - Icges(5,3) * t172;
t213 = -t63 - t65;
t66 = Icges(6,4) * t115 + Icges(6,2) * t114 + Icges(6,6) * t176;
t68 = Icges(5,4) * t115 + Icges(5,2) * t114 + Icges(5,6) * t176;
t212 = t66 + t68;
t67 = Icges(6,4) * t117 + Icges(6,2) * t116 - Icges(6,6) * t172;
t69 = Icges(5,4) * t117 + Icges(5,2) * t116 - Icges(5,6) * t172;
t211 = -t67 - t69;
t70 = Icges(6,1) * t115 + Icges(6,4) * t114 + Icges(6,5) * t176;
t72 = Icges(5,1) * t115 + Icges(5,4) * t114 + Icges(5,5) * t176;
t210 = t70 + t72;
t71 = Icges(6,1) * t117 + Icges(6,4) * t116 - Icges(6,5) * t172;
t73 = Icges(5,1) * t117 + Icges(5,4) * t116 - Icges(5,5) * t172;
t209 = t71 + t73;
t98 = -Icges(6,3) * t151 + (Icges(6,5) * t140 - Icges(6,6) * t139) * t150;
t99 = -Icges(5,3) * t151 + (Icges(5,5) * t140 - Icges(5,6) * t139) * t150;
t208 = -t98 - t99;
t204 = t215 * t140 * t150;
t100 = -Icges(6,6) * t151 + (Icges(6,4) * t140 - Icges(6,2) * t139) * t150;
t101 = -Icges(5,6) * t151 + (Icges(5,4) * t140 - Icges(5,2) * t139) * t150;
t206 = -t100 - t101;
t205 = t206 * t139;
t198 = t150 * t205 + t208 * t151 + t204;
t207 = t198 * t151;
t154 = cos(qJ(3));
t138 = t154 * pkin(3) + pkin(2);
t127 = pkin(4) * t140 + t138;
t203 = t115 * rSges(6,1) + t114 * rSges(6,2) + rSges(6,3) * t176 + t127 * t182;
t152 = sin(qJ(3));
t192 = pkin(3) * t152;
t128 = pkin(4) * t139 + t192;
t202 = -rSges(6,1) * t117 - rSges(6,2) * t116 + t153 * t128;
t201 = t214 * t151 + (t212 * t139 - t210 * t140) * t150;
t200 = t213 * t151 + (t211 * t139 + t209 * t140) * t150;
t199 = t206 * t114 - t215 * t115 + t208 * t176;
t125 = t138 * t182;
t156 = -pkin(7) - pkin(6);
t145 = -qJ(5) + t156;
t168 = t145 - t156;
t197 = -t125 + (-t128 + t192) * t155 - t168 * t176 + t203;
t179 = t152 * t153;
t171 = pkin(3) * t179 + t138 * t181;
t185 = t127 * t151;
t196 = (t150 * t168 - t185) * t155 + t171 - rSges(6,3) * t172 - t202;
t147 = t153 ^ 2;
t148 = t155 ^ 2;
t195 = ((t212 * t114 + t210 * t115 + t214 * t176) * t176 + t199 * t151 + (t211 * t114 - t209 * t115 + t213 * t176) * t172) * t176;
t194 = t207 + (t201 * t153 + t200 * t155) * t150;
t193 = pkin(2) * t151;
t191 = pkin(6) + t156;
t189 = t196 * t151;
t75 = t115 * rSges(5,1) + t114 * rSges(5,2) + rSges(5,3) * t176;
t160 = -t117 * rSges(5,1) - t116 * rSges(5,2);
t77 = -rSges(5,3) * t172 - t160;
t34 = t75 * t172 + t77 * t176;
t178 = t152 * t155;
t87 = -pkin(3) * t178 + t125 + (-t150 * t191 - t193) * t153;
t170 = pkin(2) * t181 + pkin(6) * t172;
t88 = t156 * t172 + t170 - t171;
t187 = t87 * t172 + t88 * t176;
t186 = (rSges(6,3) - t168) * t151 + (-rSges(6,1) * t140 + rSges(6,2) * t139 - t127 + t138) * t150;
t111 = -Icges(4,6) * t151 + (Icges(4,4) * t154 - Icges(4,2) * t152) * t150;
t180 = t152 * t111;
t175 = t153 * t154;
t174 = t153 * t155;
t173 = t154 * t155;
t169 = t155 * pkin(1) + t153 * qJ(2);
t167 = t148 + t147;
t121 = -t151 * t179 - t173;
t122 = t151 * t175 - t178;
t89 = t122 * rSges(4,1) + t121 * rSges(4,2) + rSges(4,3) * t176;
t164 = t150 * t186;
t105 = -rSges(5,3) * t151 + (rSges(5,1) * t140 - rSges(5,2) * t139) * t150;
t97 = t191 * t151 + (-pkin(2) + t138) * t150;
t163 = t150 * (-t105 - t97);
t10 = t197 * t172 + t196 * t176;
t162 = t150 * (-t97 + t186);
t161 = rSges(3,1) * t151 - rSges(3,2) * t150;
t32 = t100 * t116 + t102 * t117 - t172 * t98;
t33 = t101 * t116 + t103 * t117 - t172 * t99;
t158 = (-t199 - t201) * t176 / 0.2e1 - (t32 + t33 + t200) * t172 / 0.2e1;
t123 = t151 * t178 - t175;
t124 = -t151 * t173 - t179;
t90 = rSges(4,1) * t124 + rSges(4,2) * t123 - rSges(4,3) * t172;
t3 = -t32 * t151 + (t116 * t66 + t117 * t70 - t172 * t62) * t176 - (t116 * t67 + t117 * t71 - t172 * t63) * t172;
t4 = -t33 * t151 + (t116 * t68 + t117 * t72 - t172 * t64) * t176 - (t116 * t69 + t117 * t73 - t172 * t65) * t172;
t157 = t194 * t151 + (-t3 - t4) * t172 + t195;
t142 = t153 * pkin(1);
t130 = rSges(2,1) * t155 - rSges(2,2) * t153;
t129 = rSges(2,1) * t153 + rSges(2,2) * t155;
t113 = -rSges(4,3) * t151 + (rSges(4,1) * t154 - rSges(4,2) * t152) * t150;
t112 = -Icges(4,5) * t151 + (Icges(4,1) * t154 - Icges(4,4) * t152) * t150;
t110 = -Icges(4,3) * t151 + (Icges(4,5) * t154 - Icges(4,6) * t152) * t150;
t96 = t150 * t154 * t112;
t95 = rSges(3,3) * t153 + t155 * t161 + t169;
t94 = t142 + (-rSges(3,3) - qJ(2)) * t155 + t161 * t153;
t86 = Icges(4,1) * t124 + Icges(4,4) * t123 - Icges(4,5) * t172;
t85 = Icges(4,1) * t122 + Icges(4,4) * t121 + Icges(4,5) * t176;
t84 = Icges(4,4) * t124 + Icges(4,2) * t123 - Icges(4,6) * t172;
t83 = Icges(4,4) * t122 + Icges(4,2) * t121 + Icges(4,6) * t176;
t82 = Icges(4,5) * t124 + Icges(4,6) * t123 - Icges(4,3) * t172;
t81 = Icges(4,5) * t122 + Icges(4,6) * t121 + Icges(4,3) * t176;
t80 = t151 * t88;
t61 = t151 * t77;
t55 = -t90 + t169 + t170;
t54 = -t155 * qJ(2) + t142 + (pkin(6) * t150 + t193) * t153 + t89;
t51 = -t113 * t172 + t151 * t90;
t50 = -t113 * t176 - t151 * t89;
t46 = -t151 * t110 - t150 * t180 + t96;
t45 = -t105 * t172 + t61;
t44 = -t105 * t176 - t151 * t75;
t43 = (rSges(5,3) - t156) * t172 + t160 + t169 + t171;
t42 = -t156 * t176 + t125 + t142 + (-qJ(2) - t192) * t155 + t75;
t39 = (t153 * t90 + t155 * t89) * t150;
t38 = -t110 * t172 + t111 * t123 + t112 * t124;
t37 = t110 * t176 + t111 * t121 + t112 * t122;
t36 = (t185 + (rSges(6,3) - t145) * t150) * t155 + t169 + t202;
t35 = -t145 * t176 + t142 + (-qJ(2) - t128) * t155 + t203;
t25 = -t151 * t82 + (-t152 * t84 + t154 * t86) * t150;
t24 = -t151 * t81 + (-t152 * t83 + t154 * t85) * t150;
t23 = t155 * t163 + t61 + t80;
t22 = (-t75 - t87) * t151 + t153 * t163;
t13 = t155 * t164 + t189;
t12 = -t151 * t197 + t153 * t164;
t11 = t187 + t34;
t9 = t155 * t162 + t189 + t80;
t8 = (-t87 - t197) * t151 + t153 * t162;
t7 = t10 + t187;
t1 = [Icges(2,3) + t96 + (Icges(3,2) * t151 - t110 + t208) * t151 + (Icges(3,1) * t150 + 0.2e1 * Icges(3,4) * t151 - t180 + t205) * t150 + m(2) * (t129 ^ 2 + t130 ^ 2) + m(3) * (t94 ^ 2 + t95 ^ 2) + m(4) * (t54 ^ 2 + t55 ^ 2) + m(5) * (t42 ^ 2 + t43 ^ 2) + m(6) * (t35 ^ 2 + t36 ^ 2) + t204; m(3) * (-t153 * t94 - t155 * t95) + m(4) * (-t153 * t54 - t155 * t55) + m(5) * (-t153 * t42 - t155 * t43) + m(6) * (-t153 * t35 - t155 * t36); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t167; (-t46 - t198) * t151 + m(4) * (t50 * t54 + t51 * t55) + m(5) * (t22 * t42 + t23 * t43) + m(6) * (t35 * t8 + t36 * t9) + ((-t25 / 0.2e1 - t38 / 0.2e1) * t155 + (t24 / 0.2e1 + t37 / 0.2e1) * t153) * t150 + t158; m(4) * (-t153 * t50 - t155 * t51) + m(5) * (-t153 * t22 - t155 * t23) + m(6) * (-t153 * t8 - t155 * t9); (t46 * t151 + t194) * t151 + (-t155 * t4 - t155 * t3 + (t153 * ((t121 * t83 + t122 * t85) * t153 - (t121 * t84 + t122 * t86) * t155) - t155 * ((t123 * t83 + t124 * t85) * t153 - (t123 * t84 + t124 * t86) * t155) + (t153 * (t147 * t81 - t174 * t82) - t155 * (t148 * t82 - t174 * t81)) * t150) * t150 + ((t25 + t38) * t155 + (-t24 - t37) * t153) * t151) * t150 + m(6) * (t7 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(5) * (t11 ^ 2 + t22 ^ 2 + t23 ^ 2) + m(4) * (t39 ^ 2 + t50 ^ 2 + t51 ^ 2) + t195; -t207 + m(5) * (t42 * t44 + t43 * t45) + m(6) * (t12 * t35 + t13 * t36) + t158; m(5) * (-t153 * t44 - t155 * t45) + m(6) * (-t12 * t153 - t13 * t155); m(6) * (t10 * t7 + t12 * t8 + t13 * t9) + m(5) * (t11 * t34 + t22 * t44 + t23 * t45) + t157; m(5) * (t34 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t10 ^ 2 + t12 ^ 2 + t13 ^ 2) + t157; m(6) * (t153 * t36 - t155 * t35) * t150; 0; m(6) * (-t151 * t7 + (t153 * t9 - t155 * t8) * t150); m(6) * (-t10 * t151 + (-t12 * t155 + t13 * t153) * t150); m(6) * (t150 ^ 2 * t167 + t151 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
