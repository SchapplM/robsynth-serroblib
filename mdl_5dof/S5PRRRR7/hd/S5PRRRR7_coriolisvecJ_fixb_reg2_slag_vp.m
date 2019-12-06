% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:12:59
% EndTime: 2019-12-05 17:13:08
% DurationCPUTime: 2.27s
% Computational Cost: add. (3677->281), mult. (9236->398), div. (0->0), fcn. (6653->8), ass. (0->163)
t138 = sin(qJ(5));
t142 = cos(qJ(5));
t143 = cos(qJ(4));
t144 = cos(qJ(3));
t197 = t143 * t144;
t173 = qJD(2) * t197;
t139 = sin(qJ(4));
t140 = sin(qJ(3));
t189 = qJD(2) * t140;
t174 = t139 * t189;
t104 = -t173 + t174;
t217 = t104 * pkin(8);
t141 = sin(qJ(2));
t181 = t141 * qJD(1);
t205 = qJD(2) * pkin(6);
t124 = t181 + t205;
t165 = pkin(7) * qJD(2) + t124;
t96 = t165 * t144;
t87 = t143 * t96;
t95 = t165 * t140;
t88 = qJD(3) * pkin(3) - t95;
t46 = t139 * t88 + t87;
t32 = t46 - t217;
t203 = t142 * t32;
t135 = qJD(3) + qJD(4);
t114 = t139 * t144 + t143 * t140;
t106 = t114 * qJD(2);
t100 = t106 * pkin(8);
t85 = t139 * t96;
t45 = t143 * t88 - t85;
t31 = -t100 + t45;
t29 = t135 * pkin(4) + t31;
t14 = t138 * t29 + t203;
t184 = qJD(4) * t143;
t185 = qJD(4) * t139;
t145 = cos(qJ(2));
t180 = t145 * qJD(1);
t187 = qJD(3) * t140;
t77 = -t124 * t187 + (-pkin(7) * t187 + t144 * t180) * qJD(2);
t186 = qJD(3) * t144;
t78 = -t124 * t186 + (-pkin(7) * t186 - t140 * t180) * qJD(2);
t166 = -t139 * t78 - t143 * t77 - t88 * t184 + t96 * t185;
t223 = t135 * t114;
t65 = t223 * qJD(2);
t11 = -t65 * pkin(8) - t166;
t167 = -t139 * t77 + t143 * t78;
t22 = -qJD(4) * t46 + t167;
t199 = t139 * t140;
t158 = t135 * t199;
t179 = qJD(2) * qJD(3);
t170 = t144 * t179;
t193 = -qJD(4) * t173 - t143 * t170;
t64 = qJD(2) * t158 + t193;
t12 = t64 * pkin(8) + t22;
t171 = -t138 * t11 + t142 * t12;
t2 = -qJD(5) * t14 + t171;
t157 = t138 * t104 - t142 * t106;
t133 = -t144 * pkin(3) - pkin(2);
t110 = t133 * qJD(2) - t180;
t74 = t104 * pkin(4) + t110;
t212 = t74 * t157;
t228 = t2 + t212;
t183 = qJD(5) * t138;
t169 = t138 * t12 - t32 * t183;
t1 = (qJD(5) * t29 + t11) * t142 + t169;
t57 = t142 * t104 + t138 * t106;
t211 = t74 * t57;
t227 = t211 - t1;
t220 = -pkin(7) - pkin(6);
t175 = qJD(3) * t220;
t118 = t140 * t175;
t119 = t144 * t175;
t120 = t220 * t140;
t121 = t220 * t144;
t113 = -t197 + t199;
t152 = t113 * t145;
t208 = qJD(1) * t152 + t143 * t118 + t139 * t119 + t120 * t184 + t121 * t185;
t153 = t114 * t145;
t81 = t139 * t120 - t143 * t121;
t207 = qJD(1) * t153 - t81 * qJD(4) - t139 * t118 + t143 * t119;
t226 = -pkin(8) * t223 + t208;
t75 = -t143 * t186 - t144 * t184 + t158;
t225 = t75 * pkin(8) + t207;
t224 = t57 * t157;
t10 = t157 ^ 2 - t57 ^ 2;
t134 = qJD(5) + t135;
t182 = qJD(5) * t142;
t19 = t104 * t182 + t106 * t183 + t138 * t65 + t142 * t64;
t7 = t57 * t134 - t19;
t149 = qJD(5) * t157 + t138 * t64 - t142 * t65;
t8 = -t134 * t157 + t149;
t222 = -0.2e1 * t179;
t98 = t113 * t141;
t146 = qJD(3) ^ 2;
t147 = qJD(2) ^ 2;
t221 = t141 * (t146 + t147);
t177 = pkin(3) * t187;
t154 = t177 - t181;
t80 = t143 * t120 + t139 * t121;
t51 = -t114 * pkin(8) + t80;
t52 = -t113 * pkin(8) + t81;
t27 = -t138 * t52 + t142 * t51;
t219 = qJD(5) * t27 + t138 * t225 + t142 * t226;
t28 = t138 * t51 + t142 * t52;
t218 = -qJD(5) * t28 - t138 * t226 + t142 * t225;
t216 = t106 * pkin(4);
t132 = t143 * pkin(3) + pkin(4);
t198 = t139 * t142;
t47 = t139 * t95 - t87;
t35 = t47 + t217;
t48 = -t143 * t95 - t85;
t36 = -t100 + t48;
t210 = -t138 * t36 + t142 * t35 + t132 * t183 - (-t139 * t182 + (-t138 * t143 - t198) * qJD(4)) * pkin(3);
t200 = t138 * t139;
t209 = t138 * t35 + t142 * t36 - t132 * t182 - (-t139 * t183 + (t142 * t143 - t200) * qJD(4)) * pkin(3);
t206 = qJD(2) * pkin(2);
t204 = t138 * t32;
t202 = t106 * t104;
t201 = t110 * t106;
t196 = t146 * t140;
t195 = t146 * t144;
t194 = t147 * t145;
t111 = (t177 + t181) * qJD(2);
t136 = t140 ^ 2;
t137 = t144 ^ 2;
t192 = t136 - t137;
t191 = t136 + t137;
t188 = qJD(2) * t141;
t178 = pkin(3) * t189;
t176 = t140 * t147 * t144;
t172 = -pkin(4) * t134 - t29;
t163 = pkin(4) * t223 + t154;
t125 = -t180 - t206;
t162 = -t125 - t180;
t161 = t145 * t222;
t160 = t140 * t170;
t13 = t142 * t29 - t204;
t159 = -t13 * t57 - t14 * t157;
t97 = t114 * t141;
t49 = t138 * t98 - t142 * t97;
t50 = -t138 * t97 - t142 * t98;
t73 = -t138 * t113 + t142 * t114;
t156 = qJD(2) * t162;
t155 = t110 * t104 + t166;
t151 = qJD(3) * (-t162 - t206);
t102 = pkin(3) * t198 + t138 * t132;
t101 = -pkin(3) * t200 + t142 * t132;
t92 = t113 * pkin(4) + t133;
t83 = t178 + t216;
t72 = t142 * t113 + t138 * t114;
t44 = -t104 ^ 2 + t106 ^ 2;
t41 = t65 * pkin(4) + t111;
t40 = t106 * t135 - t65;
t39 = -t193 + (t104 - t174) * t135;
t34 = -qJD(2) * t153 + t135 * t98;
t33 = -qJD(2) * t152 - t141 * t223;
t24 = qJD(5) * t73 - t138 * t75 + t142 * t223;
t23 = t113 * t182 + t114 * t183 + t138 * t223 + t142 * t75;
t16 = t142 * t31 - t204;
t15 = -t138 * t31 - t203;
t6 = -qJD(5) * t50 - t138 * t33 + t142 * t34;
t5 = qJD(5) * t49 + t138 * t34 + t142 * t33;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147 * t141, -t194, 0, 0, 0, 0, 0, 0, 0, 0, t140 * t161 - t144 * t221, t140 * t221 + t144 * t161, t191 * t194, (t125 * t141 + (-t181 + (t124 + t181) * t191) * t145) * qJD(2), 0, 0, 0, 0, 0, 0, t104 * t188 + t34 * t135 - t145 * t65, t106 * t188 - t33 * t135 + t145 * t64, -t33 * t104 - t34 * t106 - t97 * t64 + t98 * t65, t110 * t188 - t111 * t145 + t166 * t98 - t22 * t97 + t46 * t33 + t45 * t34, 0, 0, 0, 0, 0, 0, t6 * t134 + t145 * t149 + t57 * t188, -t5 * t134 + t145 * t19 - t157 * t188, t149 * t50 + t157 * t6 + t49 * t19 - t5 * t57, t1 * t50 + t13 * t6 + t14 * t5 - t41 * t145 + t188 * t74 + t2 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t160, t192 * t222, t195, -0.2e1 * t160, -t196, 0, -pkin(6) * t195 + t140 * t151, pkin(6) * t196 + t144 * t151, 0, ((-t125 - t206) * t141 + (-t124 + t205) * t145 * t191) * qJD(1), -t106 * t75 - t64 * t114, t75 * t104 - t106 * t223 + t64 * t113 - t114 * t65, -t75 * t135, t104 * t223 + t65 * t113, -t223 * t135, 0, t154 * t104 + t110 * t223 + t111 * t113 + t133 * t65 + t207 * t135, t154 * t106 - t110 * t75 + t111 * t114 - t133 * t64 - t208 * t135, -t208 * t104 - t207 * t106 + t113 * t166 - t22 * t114 - t223 * t46 + t45 * t75 + t80 * t64 - t81 * t65, t154 * t110 + t111 * t133 - t166 * t81 + t207 * t45 + t208 * t46 + t22 * t80, t157 * t23 - t19 * t73, t149 * t73 + t157 * t24 + t19 * t72 + t23 * t57, -t23 * t134, -t149 * t72 + t57 * t24, -t24 * t134, 0, t218 * t134 - t149 * t92 + t163 * t57 + t74 * t24 + t41 * t72, -t219 * t134 - t157 * t163 - t92 * t19 - t74 * t23 + t41 * t73, -t1 * t72 + t13 * t23 - t14 * t24 + t149 * t28 + t157 * t218 + t27 * t19 - t2 * t73 - t219 * t57, t1 * t28 + t218 * t13 + t219 * t14 + t163 * t74 + t2 * t27 + t41 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176, t192 * t147, 0, t176, 0, 0, t140 * t156, t144 * t156, 0, 0, t202, t44, t39, -t202, t40, 0, -t104 * t178 - t201 - t47 * t135 + (-t87 + (-pkin(3) * t135 - t88) * t139) * qJD(4) + t167, t48 * t135 + (-t106 * t189 - t135 * t184) * pkin(3) + t155, (t46 + t47) * t106 + (-t45 + t48) * t104 + (-t139 * t65 + t143 * t64 + (-t104 * t143 + t106 * t139) * qJD(4)) * pkin(3), -t45 * t47 - t46 * t48 + (-t110 * t189 - t139 * t166 + t143 * t22 + (-t139 * t45 + t143 * t46) * qJD(4)) * pkin(3), -t224, t10, t7, t224, t8, 0, -t210 * t134 - t83 * t57 + t228, t209 * t134 + t157 * t83 + t227, t101 * t19 + t102 * t149 - t157 * t210 + t209 * t57 + t159, t1 * t102 + t2 * t101 - t210 * t13 - t209 * t14 - t74 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, t44, t39, -t202, t40, 0, t46 * t135 - t201 + t22, t45 * t135 + t155, 0, 0, -t224, t10, t7, t224, t8, 0, -t57 * t216 - t15 * t134 + t212 + (t138 * t172 - t203) * qJD(5) + t171, t157 * t216 + t16 * t134 + t211 + (qJD(5) * t172 - t11) * t142 - t169, -t15 * t157 + t16 * t57 + (t138 * t149 + t142 * t19 + (-t138 * t157 - t142 * t57) * qJD(5)) * pkin(4) + t159, -t13 * t15 - t14 * t16 + (t1 * t138 - t106 * t74 + t142 * t2 + (-t13 * t138 + t14 * t142) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t224, t10, t7, t224, t8, 0, t14 * t134 + t228, t13 * t134 + t227, 0, 0;];
tauc_reg = t3;
