% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:30
% EndTime: 2019-12-31 21:49:34
% DurationCPUTime: 1.73s
% Computational Cost: add. (2789->303), mult. (4190->325), div. (0->0), fcn. (2234->12), ass. (0->184)
t114 = qJDD(1) + qJDD(2);
t107 = qJDD(3) + t114;
t120 = sin(qJ(3));
t124 = cos(qJ(3));
t195 = qJD(3) * t124;
t125 = cos(qJ(2));
t220 = pkin(1) * qJD(1);
t178 = qJD(2) * t220;
t121 = sin(qJ(2));
t188 = qJDD(1) * t121;
t245 = -pkin(1) * t188 - t125 * t178;
t235 = pkin(1) * t125;
t104 = qJDD(1) * t235;
t182 = t121 * t220;
t246 = -pkin(2) * t114 + qJD(3) * t182 + t121 * t178 - t104;
t115 = qJD(1) + qJD(2);
t181 = t125 * t220;
t64 = pkin(2) * t115 + t181;
t157 = t246 * t120 + t245 * t124 - t64 * t195;
t11 = pkin(8) * t107 - t157;
t119 = sin(qJ(4));
t116 = t119 ^ 2;
t123 = cos(qJ(4));
t117 = t123 ^ 2;
t197 = t116 + t117;
t247 = t197 * t11;
t204 = t107 * t117;
t205 = t107 * t116;
t244 = t204 + t205;
t187 = qJDD(4) * qJ(5);
t108 = qJD(3) + t115;
t40 = t120 * t64 + t124 * t182;
t33 = pkin(8) * t108 + t40;
t213 = t119 * t33;
t9 = t123 * t11;
t4 = t187 + t9 + (qJD(5) - t213) * qJD(4);
t191 = qJD(4) * t123;
t208 = qJDD(4) * pkin(4);
t240 = t33 * t191 - t208;
t7 = t119 * t11;
t5 = qJDD(5) + t7 + t240;
t243 = t5 * t119 + t4 * t123;
t118 = qJ(1) + qJ(2);
t112 = qJ(3) + t118;
t100 = cos(t112);
t99 = sin(t112);
t222 = g(1) * t100 + g(2) * t99;
t241 = t197 * t33;
t111 = cos(t118);
t110 = sin(t118);
t96 = g(1) * t110;
t239 = -g(2) * t111 + t96;
t192 = qJD(4) * t119;
t166 = qJD(5) + t213;
t219 = qJD(4) * pkin(4);
t21 = t166 - t219;
t194 = qJD(4) * qJ(5);
t211 = t123 * t33;
t22 = t194 + t211;
t238 = t21 * t191 - t22 * t192 + t243;
t193 = qJD(4) * t108;
t198 = -t116 + t117;
t199 = t123 * t107;
t237 = 0.2e1 * t119 * t199 + 0.2e1 * t198 * t193;
t236 = pkin(3) * t99;
t91 = g(1) * t99;
t234 = pkin(2) * t110;
t233 = pkin(2) * t124;
t232 = pkin(3) * t107;
t231 = pkin(3) * t108;
t230 = pkin(4) * t123;
t127 = qJD(4) ^ 2;
t229 = pkin(8) * t127;
t228 = g(2) * t100;
t196 = qJD(3) * t120;
t180 = pkin(2) * t196;
t190 = qJD(5) * t119;
t54 = pkin(4) * t192 - qJ(5) * t191 - t190;
t41 = t54 + t180;
t200 = t121 * t124;
t145 = t120 * t125 + t200;
t51 = t145 * t220;
t226 = t41 - t51;
t225 = t54 - t40;
t207 = t100 * t119;
t212 = t119 * t99;
t224 = g(1) * t212 - g(2) * t207;
t223 = t100 * pkin(3) + t99 * pkin(8);
t103 = pkin(2) + t235;
t56 = pkin(1) * t200 + t120 * t103;
t221 = g(1) * t111 + g(2) * t110;
t150 = qJ(5) * t119 + t230;
t66 = -pkin(3) - t150;
t218 = t107 * t66;
t201 = t120 * t121;
t23 = t103 * t195 + (-t121 * t196 + (t124 * t125 - t201) * qJD(2)) * pkin(1);
t217 = t108 * t23;
t24 = t103 * t196 + (t145 * qJD(2) + t121 * t195) * pkin(1);
t216 = t108 * t24;
t215 = t108 * t40;
t214 = t108 * t66;
t50 = pkin(8) + t56;
t210 = t127 * t50;
t209 = pkin(8) * qJDD(4);
t101 = pkin(2) * t120 + pkin(8);
t206 = t101 * t127;
t203 = t108 * t119;
t202 = t108 * t123;
t189 = qJDD(4) * t50;
t186 = qJDD(4) * t101;
t86 = t120 * t182;
t39 = t124 * t64 - t86;
t74 = t123 * t91;
t185 = t39 * t192 + t40 * t202 + t74;
t52 = t124 * t181 - t86;
t184 = t52 * t192 + t51 * t202 + t74;
t98 = pkin(2) * t111;
t183 = t98 + t223;
t179 = pkin(2) * t195;
t88 = t100 * pkin(8);
t175 = t88 - t236;
t174 = t108 * t196;
t173 = t108 * t195;
t14 = t245 * t120 - t246 * t124 - t64 * t196;
t12 = -t14 - t232;
t172 = -t12 - t228;
t55 = -pkin(1) * t201 + t103 * t124;
t35 = -t55 + t66;
t171 = t108 * t35 - t23;
t32 = -t39 - t231;
t169 = t12 * t119 + t32 * t191 - t224;
t168 = qJ(5) * t207 + t100 * t230 + t223;
t167 = t108 * t197;
t165 = qJD(1) * (-qJD(2) + t115);
t164 = qJD(2) * (-qJD(1) - t115);
t160 = t191 * t203;
t159 = t104 + t239;
t158 = t98 + t168;
t156 = -t51 + t180;
t122 = sin(qJ(1));
t154 = -pkin(1) * t122 - t234;
t153 = t229 - t232;
t126 = cos(qJ(1));
t152 = g(1) * t122 - g(2) * t126;
t151 = g(1) * t207 + g(2) * t212 - g(3) * t123 - t7;
t149 = pkin(4) * t119 - qJ(5) * t123;
t148 = t21 * t119 + t22 * t123;
t147 = t197 * t217 + t244 * t50 - t222;
t102 = -pkin(3) - t233;
t146 = t102 * t107 + t206;
t144 = t154 + t88;
t1 = (t149 * qJD(4) - t190) * t108 + t218 - t14;
t143 = -t1 - t218 - t229;
t59 = t66 - t233;
t142 = -t107 * t59 - t1 - t206;
t141 = t108 * t59 - t179;
t140 = t157 + t222;
t139 = t102 * t108 - t179;
t138 = -qJDD(5) + t151;
t49 = -pkin(3) - t55;
t137 = t107 * t49 + t210 + t216;
t136 = t14 + t91 - t228;
t135 = t66 * t91;
t16 = t24 + t54;
t134 = -t107 * t35 - t108 * t16 - t1 - t210;
t133 = -t189 + (t108 * t49 - t23) * qJD(4);
t132 = t244 * pkin(8) - t39 * t167 - t222;
t131 = t197 * pkin(2) * t173 + t244 * t101 - t52 * t167 - t222;
t130 = -g(1) * t88 - t135;
t129 = (-t119 * t22 + t123 * t21) * qJD(4) + t243;
t113 = t126 * pkin(1);
t106 = t108 ^ 2;
t83 = t119 * t107;
t68 = qJDD(4) * t123 - t127 * t119;
t67 = qJDD(4) * t119 + t123 * t127;
t65 = t119 * t106 * t123;
t53 = t198 * t106;
t47 = t149 * t108;
t46 = -0.2e1 * t160 + t204;
t45 = 0.2e1 * t160 + t205;
t26 = t32 * t192;
t17 = -t39 + t214;
t15 = t17 * t192;
t2 = [0, 0, 0, 0, 0, qJDD(1), t152, g(1) * t126 + g(2) * t122, 0, 0, 0, 0, 0, 0, 0, t114, (t114 * t125 + t121 * t164) * pkin(1) + t159, ((-qJDD(1) - t114) * t121 + t125 * t164) * pkin(1) + t221, 0, (t152 + (t121 ^ 2 + t125 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t107, t107 * t55 + t136 - t216, -t107 * t56 + t140 - t217, 0, -t157 * t56 + t40 * t23 + t14 * t55 - t39 * t24 - g(1) * t154 - g(2) * (t98 + t113), t45, t237, t67, t46, t68, 0, t26 + t74 + t133 * t119 + (-t137 + t172) * t123, t119 * t137 + t123 * t133 + t169, t147 + t247, t12 * t49 + t32 * t24 - g(1) * (t144 - t236) - g(2) * (t113 + t183) + t23 * t241 + t50 * t247, t45, t67, -t237, 0, -t68, t46, t15 + t74 + (qJD(4) * t171 - t189) * t119 + (t134 - t228) * t123, t147 + t238, (t189 + (-t17 - t171) * qJD(4)) * t123 + t134 * t119 + t224, t1 * t35 + t17 * t16 - g(1) * t144 - g(2) * (t113 + t158) - t135 + t148 * t23 + t129 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, t121 * pkin(1) * t165 + t159, (t125 * t165 - t188) * pkin(1) + t221, 0, 0, 0, 0, 0, 0, 0, t107, t108 * t51 + (t107 * t124 - t174) * pkin(2) + t136, t108 * t52 + (-t107 * t120 - t173) * pkin(2) + t140, 0, t39 * t51 - t40 * t52 + (-t120 * t157 + t124 * t14 + (-t120 * t39 + t124 * t40) * qJD(3) + t239) * pkin(2), t45, t237, t67, t46, t68, 0, t26 + (qJD(4) * t139 - t186) * t119 + (-pkin(2) * t174 - t146 + t172) * t123 + t184, (-t186 + (t139 + t52) * qJD(4)) * t123 + (t108 * t156 + t146) * t119 + t169, t131 + t247, t12 * t102 - g(1) * (t175 - t234) - g(2) * t183 + t156 * t32 + t101 * t247 + (t179 - t52) * t241, t45, t67, -t237, 0, -t68, t46, t15 + (qJD(4) * t141 - t186) * t119 + (-t108 * t41 + t142 - t228) * t123 + t184, t131 + t238, (t186 + (-t141 - t17 - t52) * qJD(4)) * t123 + (-t108 * t226 + t142) * t119 + t224, t1 * t59 - g(2) * t158 - t148 * t52 + t226 * t17 + (t148 * t195 + t96) * pkin(2) + t129 * t101 + t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, t136 + t215, t108 * t39 + t140, 0, 0, t45, t237, t67, t46, t68, 0, t26 + (-pkin(3) * t193 - t209) * t119 + (-t153 + t172) * t123 + t185, (-t209 + (t39 - t231) * qJD(4)) * t123 + (t153 - t215) * t119 + t169, t132 + t247, -t12 * pkin(3) + pkin(8) * t247 - g(1) * t175 - g(2) * t223 - t39 * t241 - t32 * t40, t45, t67, -t237, 0, -t68, t46, t15 + (t193 * t66 - t209) * t119 + (-t108 * t54 + t143 - t228) * t123 + t185, t132 + t238, (t209 + (-t17 - t39 - t214) * qJD(4)) * t123 + (-t108 * t225 + t143) * t119 + t224, pkin(8) * t129 - g(2) * t168 + t1 * t66 - t148 * t39 + t17 * t225 + t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t53, t83, t65, t199, qJDD(4), -t203 * t32 + t151, g(3) * t119 - t9 + (-t108 * t32 + t222) * t123, 0, 0, -t65, t83, t53, qJDD(4), -t199, t65, 0.2e1 * t208 + (-t119 * t17 + t123 * t47) * t108 + t138, -t149 * t107 + ((t22 - t194) * t119 + (qJD(5) - t21 - t219) * t123) * t108, 0.2e1 * t187 + 0.2e1 * qJD(4) * qJD(5) + t9 + (t108 * t47 - g(3)) * t119 + (t108 * t17 - t222) * t123, -t5 * pkin(4) - g(3) * t150 + t4 * qJ(5) + t222 * t149 + t166 * t22 - t17 * t47 - t21 * t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t65, t83, -t106 * t116 - t127, -qJD(4) * t22 + t17 * t203 - t138 + t240;];
tau_reg = t2;
