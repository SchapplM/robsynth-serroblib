% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR9_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:08:00
% EndTime: 2019-12-31 19:08:10
% DurationCPUTime: 3.97s
% Computational Cost: add. (7051->342), mult. (18848->459), div. (0->0), fcn. (14609->8), ass. (0->182)
t131 = cos(qJ(5));
t128 = sin(qJ(5));
t176 = qJD(5) * t131;
t127 = cos(pkin(9));
t126 = sin(pkin(9));
t130 = sin(qJ(3));
t185 = t130 * t126;
t213 = cos(qJ(3));
t147 = -t213 * t127 + t185;
t103 = t147 * qJD(1);
t184 = t130 * t127;
t110 = t213 * t126 + t184;
t104 = t110 * qJD(1);
t129 = sin(qJ(4));
t212 = cos(qJ(4));
t164 = qJD(4) * t212;
t178 = qJD(4) * t129;
t165 = qJD(3) * t213;
t117 = t127 * qJD(1) * t165;
t167 = qJD(1) * t185;
t96 = qJD(3) * t167 - t117;
t221 = t110 * qJD(3);
t97 = qJD(1) * t221;
t145 = -t103 * t164 - t104 * t178 - t129 * t97 - t212 * t96;
t146 = t129 * t103 - t212 * t104;
t173 = qJD(3) + qJD(4);
t66 = t128 * t173 - t131 * t146;
t188 = qJD(5) * t66;
t27 = t128 * t145 + t188;
t157 = t131 * t173;
t64 = -t128 * t146 - t157;
t199 = -t128 * t27 - t64 * t176;
t77 = -t212 * t103 - t129 * t104;
t233 = qJD(5) - t77;
t244 = t128 * t233;
t225 = t66 * t244;
t236 = t131 * t77;
t177 = qJD(5) * t128;
t26 = -qJD(5) * t157 - t131 * t145 - t146 * t177;
t245 = -t131 * t26 + t236 * t64 + t199 - t225;
t200 = pkin(6) + qJ(2);
t114 = t200 * t126;
t111 = qJD(1) * t114;
t115 = t200 * t127;
t112 = qJD(1) * t115;
t82 = -t130 * t111 + t213 * t112;
t63 = -t103 * pkin(7) + t82;
t171 = t212 * t63;
t186 = t130 * t112;
t81 = -t213 * t111 - t186;
t62 = -t104 * pkin(7) + t81;
t61 = qJD(3) * pkin(3) + t62;
t33 = t129 * t61 + t171;
t31 = t173 * pkin(8) + t33;
t121 = -t127 * pkin(2) - pkin(1);
t113 = t121 * qJD(1) + qJD(2);
t83 = t103 * pkin(3) + t113;
t34 = -pkin(4) * t77 + pkin(8) * t146 + t83;
t150 = t128 * t31 - t131 * t34;
t243 = t233 * t150;
t224 = t77 * t173;
t242 = t145 - t224;
t203 = t146 ^ 2;
t204 = t77 ^ 2;
t241 = t203 - t204;
t10 = t128 * t34 + t131 * t31;
t174 = qJD(1) * qJD(2);
t163 = t126 * t174;
t166 = qJD(2) * t213;
t119 = t127 * t166;
t189 = qJD(1) * t119 - t111 * t165;
t58 = (-qJD(3) * t112 - t163) * t130 + t189;
t51 = -t97 * pkin(7) + t58;
t142 = t110 * qJD(2);
t137 = qJD(1) * t142;
t59 = -t82 * qJD(3) - t137;
t52 = t96 * pkin(7) + t59;
t134 = -t129 * t52 - t61 * t164 + t63 * t178 - t212 * t51;
t215 = t97 * pkin(3);
t161 = t129 * t96 - t212 * t97;
t231 = t146 * qJD(4);
t45 = -t161 - t231;
t18 = t45 * pkin(4) - pkin(8) * t145 + t215;
t3 = -qJD(5) * t10 + t128 * t134 + t131 * t18;
t226 = t233 * t10 + t3;
t23 = t26 * t128;
t239 = -t23 + (t176 - t236) * t66;
t41 = t128 * t45;
t198 = t176 * t233 + t41;
t206 = t66 * t146;
t238 = -t233 * t236 + t198 + t206;
t192 = t129 * t63;
t32 = t212 * t61 - t192;
t30 = -t173 * pkin(4) - t32;
t237 = t30 * t77;
t205 = t146 * t64;
t235 = t233 * t146;
t202 = t77 * t146;
t234 = t83 * t146;
t183 = t146 * qJD(3);
t232 = -t183 + t161;
t230 = -t83 * t77 + t134;
t29 = t30 * t176;
t48 = t212 * t52;
t162 = t129 * t51 - t48;
t8 = t33 * qJD(4) + t162;
t229 = -t10 * t146 + t8 * t128 + t29;
t28 = t30 * t177;
t228 = -t150 * t146 + t28;
t47 = -pkin(4) * t146 - t77 * pkin(8);
t2 = -t150 * qJD(5) + t128 * t18 - t131 * t134;
t227 = t2 + t243;
t43 = t131 * t45;
t223 = t177 * t233 - t43;
t85 = -t130 * t114 + t213 * t115;
t222 = t10 * t131 + t128 * t150;
t220 = t104 ^ 2;
t219 = qJD(3) ^ 2;
t84 = -t213 * t114 - t130 * t115;
t71 = -t110 * pkin(7) + t84;
t72 = -t147 * pkin(7) + t85;
t38 = t129 * t72 - t212 * t71;
t218 = t8 * t38;
t144 = t129 * t147;
t80 = t212 * t110 - t144;
t217 = t8 * t80;
t210 = t104 * pkin(3);
t1 = t2 * t131;
t138 = t212 * t147;
t79 = t129 * t110 + t138;
t208 = t45 * t79;
t207 = t66 * t64;
t201 = t80 * t45;
t196 = pkin(3) * qJD(4);
t193 = t128 * t64;
t25 = t27 * t131;
t187 = t104 * t103;
t181 = -t114 * t165 + t119;
t180 = t126 ^ 2 + t127 ^ 2;
t179 = qJD(3) * t130;
t175 = t126 * qJD(2);
t169 = t80 * t177;
t168 = t80 * t176;
t158 = t180 * qJD(1) ^ 2;
t156 = pkin(3) * t164;
t35 = t129 * t62 + t171;
t155 = pkin(3) * t178 - t35;
t49 = t110 * t178 + t129 * t221 + t173 * t138;
t154 = -t30 * t49 + t217;
t153 = -t233 * t49 + t201;
t152 = t10 * t128 - t131 * t150;
t39 = t129 * t71 + t212 * t72;
t89 = t147 * pkin(3) + t121;
t40 = t79 * pkin(4) - t80 * pkin(8) + t89;
t20 = t128 * t40 + t131 * t39;
t19 = -t128 * t39 + t131 * t40;
t149 = t244 * t77 - t223;
t148 = 0.2e1 * t180 * t174;
t143 = t147 * qJD(3);
t139 = pkin(3) * t221;
t122 = t129 * pkin(3) + pkin(8);
t136 = -t122 * t45 - t156 * t233 - t237;
t135 = -t152 * qJD(5) - t3 * t128 + t1;
t133 = pkin(7) * t143 - qJD(2) * t184 + t114 * t179 - t115 * t165 - t126 * t166;
t123 = -t212 * pkin(3) - pkin(4);
t100 = t103 ^ 2;
t68 = -t85 * qJD(3) - t142;
t67 = (-qJD(3) * t115 - t175) * t130 + t181;
t54 = -pkin(7) * t221 - t115 * t179 - t130 * t175 + t181;
t50 = -qJD(4) * t144 + t110 * t164 - t129 * t143 + t212 * t221;
t37 = t210 + t47;
t36 = t212 * t62 - t192;
t21 = t50 * pkin(4) + t49 * pkin(8) + t139;
t16 = t128 * t47 + t131 * t32;
t15 = -t128 * t32 + t131 * t47;
t14 = t128 * t37 + t131 * t36;
t13 = -t128 * t36 + t131 * t37;
t12 = t39 * qJD(4) + t129 * t54 - t212 * t133;
t11 = t129 * t133 + t71 * t164 - t72 * t178 + t212 * t54;
t5 = -t20 * qJD(5) - t128 * t11 + t131 * t21;
t4 = t19 * qJD(5) + t131 * t11 + t128 * t21;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, qJ(2) * t148, -t104 * t143 - t96 * t110, t103 * t143 - t104 * t221 - t110 * t97 + t96 * t147, -t147 * t219, t103 * t221 + t97 * t147, -t110 * t219, 0, t68 * qJD(3) + t113 * t221 + t121 * t97, -t67 * qJD(3) - t113 * t143 - t121 * t96, -t67 * t103 - t68 * t104 - t59 * t110 - t82 * t221 + t84 * t96 - t85 * t97 + (t81 * qJD(3) - t58) * t147, t58 * t85 + t59 * t84 + t82 * t67 + t81 * t68, t145 * t80 + t146 * t49, -t145 * t79 + t146 * t50 - t49 * t77 - t201, -t49 * t173, -t50 * t77 + t208, -t50 * t173, 0, -t12 * t173 - t139 * t77 + t79 * t215 + t89 * t45 + t83 * t50, -t11 * t173 - t139 * t146 + t145 * t89 + t80 * t215 - t83 * t49, t11 * t77 - t12 * t146 + t134 * t79 + t145 * t38 + t32 * t49 - t33 * t50 - t39 * t45 + t217, t33 * t11 - t32 * t12 - t134 * t39 + t139 * t83 + t89 * t215 + t218, -t66 * t169 + (-t26 * t80 - t49 * t66) * t131, (t128 * t66 + t131 * t64) * t49 + (t23 - t25 + (-t131 * t66 + t193) * qJD(5)) * t80, t131 * t153 - t169 * t233 - t26 * t79 + t66 * t50, t64 * t168 + (t27 * t80 - t49 * t64) * t128, -t128 * t153 - t168 * t233 - t27 * t79 - t64 * t50, t233 * t50 + t208, t12 * t64 + t128 * t154 - t150 * t50 + t19 * t45 + t233 * t5 + t38 * t27 + t29 * t80 + t3 * t79, -t10 * t50 + t12 * t66 + t131 * t154 - t2 * t79 - t20 * t45 - t233 * t4 - t38 * t26 - t28 * t80, t19 * t26 - t20 * t27 - t4 * t64 - t5 * t66 + t152 * t49 + (-t222 * qJD(5) - t2 * t128 - t3 * t131) * t80, t10 * t4 + t30 * t12 - t150 * t5 + t3 * t19 + t2 * t20 + t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t158, -qJ(2) * t158, 0, 0, 0, 0, 0, 0, 0.2e1 * t104 * qJD(3), t117 + (-t103 - t167) * qJD(3), -t100 - t220, t82 * t103 + t81 * t104, 0, 0, 0, 0, 0, 0, -t161 - t183 - 0.2e1 * t231, t145 + t224, -t203 - t204, -t146 * t32 - t33 * t77 + t215, 0, 0, 0, 0, 0, 0, t149 + t205, -t131 * t233 ^ 2 + t206 - t41, (t64 * t77 + t26) * t131 + t225 + t199, t227 * t128 + t226 * t131 + t146 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t187, -t100 + t220, t117 + (t103 - t167) * qJD(3), -t187, 0, 0, -t113 * t104 - t137, t130 * t163 + t113 * t103 + (t81 + t186) * qJD(3) - t189, 0, 0, t202, t241, t242, -t202, t232, 0, -t63 * t164 + t48 + t77 * t210 + t234 + t35 * t173 + (-qJD(4) * t61 - t173 * t196 - t51) * t129, t36 * t173 + (t104 * t146 - t173 * t164) * pkin(3) + t230, t32 * t77 - t33 * t146 + t35 * t146 - t36 * t77 + (-t212 * t145 - t129 * t45 + (-t129 * t146 + t212 * t77) * qJD(4)) * pkin(3), t32 * t35 - t33 * t36 + (-t212 * t8 - t104 * t83 - t129 * t134 + (-t129 * t32 + t212 * t33) * qJD(4)) * pkin(3), t239, t245, t238, t244 * t64 - t25, t149 - t205, t235, t123 * t27 - t13 * t233 + t155 * t64 + (-qJD(5) * t122 * t233 - t8) * t131 + t136 * t128 + t228, -t123 * t26 + (t122 * t177 + t14) * t233 + t155 * t66 + t136 * t131 + t229, t13 * t66 + t14 * t64 + t1 + (-t64 * t156 - t122 * t27 - t77 * t150 + (t122 * t66 + t150) * qJD(5)) * t131 + (t66 * t156 + t10 * t77 - t122 * t26 - t3 + (t122 * t64 - t10) * qJD(5)) * t128, -t10 * t14 + t8 * t123 + t150 * t13 - t30 * t35 + (t129 * t30 + t222 * t212) * t196 + t135 * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, t241, t242, -t202, t232, 0, t33 * qJD(3) - t162 + t234, t32 * t173 + t230, 0, 0, t239, t245, t238, t193 * t233 - t25, -t233 * t244 - t205 + t43, t235, -pkin(4) * t27 - pkin(8) * t198 - t128 * t237 - t8 * t131 - t15 * t233 - t33 * t64 + t228, pkin(4) * t26 + t223 * pkin(8) + t16 * t233 - t236 * t30 - t33 * t66 + t229, t15 * t66 + t16 * t64 + t1 + (t243 + (-t27 + t188) * pkin(8)) * t131 + ((qJD(5) * t64 - t26) * pkin(8) - t226) * t128, -t8 * pkin(4) + pkin(8) * t135 - t10 * t16 + t15 * t150 - t30 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t207, -t64 ^ 2 + t66 ^ 2, t233 * t64 - t26, -t207, t233 * t66 - t27, t45, -t30 * t66 + t226, t30 * t64 - t227, 0, 0;];
tauc_reg = t6;
