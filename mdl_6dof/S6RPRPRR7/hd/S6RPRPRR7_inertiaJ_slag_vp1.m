% Calculate joint inertia matrix for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:54:58
% EndTime: 2019-03-09 03:55:02
% DurationCPUTime: 2.04s
% Computational Cost: add. (4963->327), mult. (4791->471), div. (0->0), fcn. (4903->10), ass. (0->161)
t242 = Icges(4,3) + Icges(5,3);
t156 = qJ(3) + pkin(10);
t146 = sin(t156);
t147 = cos(t156);
t161 = sin(qJ(3));
t164 = cos(qJ(3));
t241 = Icges(4,5) * t161 + Icges(5,5) * t146 + Icges(4,6) * t164 + Icges(5,6) * t147;
t148 = qJ(5) + t156;
t143 = sin(t148);
t144 = cos(t148);
t160 = sin(qJ(6));
t163 = cos(qJ(6));
t71 = rSges(7,3) * t143 + (rSges(7,1) * t163 - rSges(7,2) * t160) * t144;
t240 = pkin(5) * t144 + pkin(9) * t143 + t71;
t162 = sin(qJ(1));
t165 = cos(qJ(1));
t239 = t242 * t162 - t241 * t165;
t238 = t241 * t162 + t242 * t165;
t237 = (rSges(6,1) * t143 + rSges(6,2) * t144) * t165;
t236 = (rSges(4,1) * t161 + rSges(4,2) * t164) * t165;
t157 = t162 ^ 2;
t158 = t165 ^ 2;
t235 = t162 / 0.2e1;
t234 = t165 / 0.2e1;
t136 = t157 + t158;
t120 = m(6) * t136;
t233 = pkin(3) * t161;
t232 = pkin(3) * t164;
t231 = pkin(5) * t143;
t159 = -qJ(4) - pkin(7);
t206 = t162 * t163;
t209 = t160 * t165;
t102 = t143 * t209 + t206;
t204 = t163 * t165;
t207 = t162 * t160;
t103 = -t143 * t204 + t207;
t187 = -t103 * rSges(7,1) - t102 * rSges(7,2);
t211 = t144 * t165;
t55 = rSges(7,3) * t211 - t187;
t230 = t165 * t55 + t158 * (pkin(9) * t144 - t231);
t213 = t143 * t162;
t125 = pkin(5) * t213;
t212 = t144 * t162;
t100 = -t143 * t207 + t204;
t101 = t143 * t206 + t209;
t226 = t101 * rSges(7,1) + t100 * rSges(7,2);
t54 = -rSges(7,3) * t212 + t226;
t229 = pkin(9) * t212 - t125 - t54;
t68 = Icges(7,3) * t143 + (Icges(7,5) * t163 - Icges(7,6) * t160) * t144;
t70 = Icges(7,5) * t143 + (Icges(7,1) * t163 - Icges(7,4) * t160) * t144;
t228 = t144 * t163 * t70 + t143 * t68;
t201 = t162 * t159 + t165 * t233;
t122 = pkin(4) * t146 + t233;
t155 = -pkin(8) + t159;
t203 = -t165 * t122 - t162 * t155;
t90 = t165 * (-t162 * pkin(7) - t201);
t227 = t165 * (t201 + t203) + t90;
t46 = t240 * t162;
t225 = rSges(5,1) * t146;
t69 = Icges(7,6) * t143 + (Icges(7,4) * t163 - Icges(7,2) * t160) * t144;
t224 = t160 * t69;
t152 = t165 * rSges(5,3);
t48 = Icges(7,5) * t101 + Icges(7,6) * t100 - Icges(7,3) * t212;
t50 = Icges(7,4) * t101 + Icges(7,2) * t100 - Icges(7,6) * t212;
t52 = Icges(7,1) * t101 + Icges(7,4) * t100 - Icges(7,5) * t212;
t21 = t143 * t48 + (-t160 * t50 + t163 * t52) * t144;
t223 = t21 * t165;
t49 = Icges(7,5) * t103 + Icges(7,6) * t102 + Icges(7,3) * t211;
t51 = Icges(7,4) * t103 + Icges(7,2) * t102 + Icges(7,6) * t211;
t53 = Icges(7,1) * t103 + Icges(7,4) * t102 + Icges(7,5) * t211;
t22 = t143 * t49 + (-t160 * t51 + t163 * t53) * t144;
t222 = t22 * t162;
t208 = t161 * t162;
t140 = pkin(3) * t208;
t106 = t140 + (-pkin(7) - t159) * t165;
t113 = t162 * t122;
t221 = -t106 - t113 + t140 - (-t155 + t159) * t165;
t219 = Icges(4,4) * t161;
t218 = Icges(4,4) * t164;
t217 = Icges(5,4) * t146;
t216 = Icges(5,4) * t147;
t215 = Icges(6,4) * t143;
t214 = Icges(6,4) * t144;
t210 = t147 * t162;
t205 = t162 * t164;
t141 = pkin(3) * t205;
t202 = pkin(4) * t210 + t141;
t200 = t165 * pkin(1) + t162 * qJ(2);
t150 = t165 * qJ(2);
t199 = t150 - t203;
t198 = t120 + (m(5) + m(7)) * t136;
t79 = rSges(6,1) * t213 + rSges(6,2) * t212 + t165 * rSges(6,3);
t197 = -rSges(5,2) * t210 - t162 * t225 - t152;
t196 = rSges(4,1) * t208 + rSges(4,2) * t205 + t165 * rSges(4,3);
t195 = (-rSges(7,3) - pkin(9)) * t144;
t17 = t102 * t50 + t103 * t52 + t211 * t48;
t18 = t102 * t51 + t103 * t53 + t211 * t49;
t10 = t18 * t162 + t165 * t17;
t171 = Icges(6,5) * t143 + Icges(6,6) * t144;
t73 = Icges(6,3) * t165 + t162 * t171;
t74 = Icges(6,3) * t162 - t165 * t171;
t15 = t100 * t50 + t101 * t52 - t212 * t48;
t16 = t100 * t51 + t101 * t53 - t212 * t49;
t9 = t15 * t165 + t16 * t162;
t194 = (t158 * t73 + t9) * t165 + (t157 * t74 + t10 + (t162 * t73 + t165 * t74) * t165) * t162;
t27 = t100 * t69 + t101 * t70 - t212 * t68;
t3 = t27 * t143 + (-t15 * t162 + t16 * t165) * t144;
t28 = t102 * t69 + t103 * t70 + t211 * t68;
t4 = t28 * t143 + (-t162 * t17 + t165 * t18) * t144;
t193 = t3 * t234 + t4 * t235 + t143 * (t222 + t223) / 0.2e1 - t9 * t212 / 0.2e1 + t10 * t211 / 0.2e1;
t191 = -pkin(4) * t147 - t232;
t189 = rSges(5,2) * t147 + t225;
t111 = rSges(6,1) * t144 - rSges(6,2) * t143;
t62 = t111 * t162 + t202;
t63 = (-t111 + t191) * t165;
t180 = t62 * t162 - t165 * t63;
t179 = Icges(4,1) * t161 + t218;
t178 = Icges(5,1) * t146 + t216;
t177 = Icges(6,1) * t143 + t214;
t176 = Icges(4,2) * t164 + t219;
t175 = Icges(5,2) * t147 + t217;
t174 = Icges(6,2) * t144 + t215;
t109 = -Icges(6,2) * t143 + t214;
t110 = Icges(6,1) * t144 - t215;
t170 = t109 * t144 + t110 * t143;
t169 = -t155 * t165 + t113 + t200;
t66 = t150 + t236 + (-rSges(4,3) - pkin(1) - pkin(7)) * t162;
t67 = pkin(7) * t165 + t196 + t200;
t168 = m(4) * (t162 * t66 - t165 * t67);
t43 = t237 + (-rSges(6,3) - pkin(1)) * t162 + t199;
t44 = t169 + t79;
t167 = m(6) * (t162 * t43 - t165 * t44);
t108 = Icges(6,5) * t144 - Icges(6,6) * t143;
t166 = t223 / 0.2e1 + t222 / 0.2e1 + (t162 * t108 - t143 * (Icges(6,6) * t162 - t165 * t174) + t144 * (Icges(6,5) * t162 - t165 * t177) - t165 * t170 + t28) * t235 + (t108 * t165 - t143 * (Icges(6,6) * t165 + t162 * t174) + t144 * (Icges(6,5) * t165 + t162 * t177) + t162 * t170 + t27) * t234;
t132 = rSges(2,1) * t165 - t162 * rSges(2,2);
t131 = rSges(4,1) * t164 - rSges(4,2) * t161;
t130 = -t162 * rSges(2,1) - rSges(2,2) * t165;
t118 = rSges(5,1) * t147 - rSges(5,2) * t146;
t105 = -rSges(3,2) * t165 + t162 * rSges(3,3) + t200;
t104 = rSges(3,3) * t165 + t150 + (rSges(3,2) - pkin(1)) * t162;
t81 = (-t118 - t232) * t165;
t80 = t118 * t162 + t141;
t72 = t165 * (t162 * rSges(6,3) - t237);
t58 = -t159 * t165 + t140 - t197 + t200;
t57 = t150 + t189 * t165 + (-rSges(5,3) - pkin(1)) * t162 + t201;
t56 = -t162 * t196 + (t162 * rSges(4,3) - t236) * t165;
t47 = t240 * t165;
t42 = -t162 * t79 + t72;
t41 = (t191 - t240) * t165;
t40 = t202 + t46;
t35 = t90 - t189 * t158 + (-t106 + t197 + t152) * t162;
t34 = -t143 * t55 + t211 * t71;
t33 = t143 * t54 + t212 * t71;
t32 = t162 * t195 + t125 + t169 + t226;
t31 = -t162 * pkin(1) + (t195 + t231) * t165 + t187 + t199;
t30 = (-t144 * t224 + t228) * t143;
t29 = (-t162 * t55 - t165 * t54) * t144;
t24 = t162 * t229 + t230;
t23 = t72 + (-t79 + t221) * t162 + t227;
t12 = (t221 + t229) * t162 + t227 + t230;
t1 = [-t143 * t109 - t146 * (-Icges(5,2) * t146 + t216) + t147 * (Icges(5,1) * t147 - t217) - t161 * (-Icges(4,2) * t161 + t218) + t164 * (Icges(4,1) * t164 - t219) + Icges(3,1) + Icges(2,3) + (t110 - t224) * t144 + m(7) * (t31 ^ 2 + t32 ^ 2) + m(6) * (t43 ^ 2 + t44 ^ 2) + m(5) * (t57 ^ 2 + t58 ^ 2) + m(4) * (t66 ^ 2 + t67 ^ 2) + m(3) * (t104 ^ 2 + t105 ^ 2) + m(2) * (t130 ^ 2 + t132 ^ 2) + t228; m(7) * (t162 * t31 - t165 * t32) + t167 + m(5) * (t162 * t57 - t165 * t58) + t168 + m(3) * (t162 * t104 - t105 * t165); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1) * t136 + t198; m(7) * (t31 * t40 + t32 * t41) + m(6) * (t43 * t62 + t44 * t63) + m(5) * (t57 * t80 + t58 * t81) + t166 + t131 * t168 + (-t146 * (Icges(5,6) * t162 - t165 * t175) + t147 * (Icges(5,5) * t162 - t165 * t178) - t161 * (Icges(4,6) * t162 - t165 * t176) + t164 * (Icges(4,5) * t162 - t165 * t179)) * t235 + (-t146 * (Icges(5,6) * t165 + t162 * t175) + t147 * (Icges(5,5) * t165 + t162 * t178) - t161 * (Icges(4,6) * t165 + t162 * t176) + t164 * (Icges(4,5) * t165 + t162 * t179)) * t234 + (Icges(4,5) * t164 + Icges(5,5) * t147 - Icges(4,6) * t161 - Icges(5,6) * t146) * (t158 / 0.2e1 + t157 / 0.2e1); m(5) * (t80 * t162 - t165 * t81) + m(6) * t180 + m(7) * (t40 * t162 - t165 * t41) + m(4) * t136 * t131; m(7) * (t12 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(6) * (t23 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(5) * (t35 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(4) * (t131 ^ 2 * t136 + t56 ^ 2) + t194 + t239 * t162 * t157 + (t238 * t158 + (t238 * t162 + t165 * t239) * t162) * t165; m(7) * (t162 * t32 + t165 * t31) + m(6) * (t162 * t44 + t165 * t43) + m(5) * (t162 * t58 + t165 * t57); 0; m(7) * (t162 * t41 + t165 * t40) + m(6) * (t162 * t63 + t165 * t62) + m(5) * (t162 * t81 + t165 * t80); t198; m(7) * (t31 * t46 - t32 * t47) + t111 * t167 + t166; m(7) * (t46 * t162 + t165 * t47) + t111 * t120; m(7) * (t12 * t24 + t40 * t46 - t41 * t47) + m(6) * (t111 * t180 + t42 * t23) + t194; m(7) * (-t47 * t162 + t165 * t46); m(6) * (t111 ^ 2 * t136 + t42 ^ 2) + m(7) * (t24 ^ 2 + t46 ^ 2 + t47 ^ 2) + t194; m(7) * (t31 * t34 + t32 * t33) + t30 + ((t22 / 0.2e1 + t28 / 0.2e1) * t165 + (-t21 / 0.2e1 - t27 / 0.2e1) * t162) * t144; m(7) * (t34 * t162 - t165 * t33); m(7) * (t12 * t29 + t33 * t41 + t34 * t40) + t193; m(7) * (t33 * t162 + t165 * t34); m(7) * (t24 * t29 - t33 * t47 + t34 * t46) + t193; t143 * t30 + m(7) * (t29 ^ 2 + t33 ^ 2 + t34 ^ 2) + (-t162 * t3 + t165 * t4 + t143 * (-t162 * t21 + t165 * t22)) * t144;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
