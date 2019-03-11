% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPPR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:42:39
% EndTime: 2019-03-09 02:42:46
% DurationCPUTime: 2.72s
% Computational Cost: add. (4202->319), mult. (9920->411), div. (0->0), fcn. (6771->8), ass. (0->187)
t127 = sin(pkin(10));
t129 = cos(pkin(10));
t132 = sin(qJ(3));
t134 = cos(qJ(3));
t107 = t127 * t134 + t129 * t132;
t143 = qJD(1) * t107;
t240 = qJD(6) + t143;
t133 = cos(qJ(6));
t131 = sin(qJ(6));
t184 = t131 * qJD(3);
t190 = qJD(1) * t134;
t174 = t129 * t190;
t191 = qJD(1) * t132;
t98 = t127 * t191 - t174;
t72 = -t133 * t98 + t184;
t169 = t240 * t72;
t186 = qJD(6) * t133;
t100 = t107 * qJD(3);
t86 = qJD(1) * t100;
t44 = qJD(6) * t184 - t131 * t86 - t98 * t186;
t243 = t44 - t169;
t181 = t134 * qJD(4);
t188 = qJD(3) * t132;
t118 = sin(pkin(9)) * pkin(1) + pkin(7);
t109 = t118 * qJD(1);
t124 = t134 * qJD(2);
t121 = qJD(3) * t124;
t82 = -t109 * t188 + t121;
t59 = (-qJ(4) * t188 + t181) * qJD(1) + t82;
t182 = t132 * qJD(4);
t164 = qJ(4) * qJD(1) + t109;
t183 = t132 * qJD(2);
t76 = t134 * t164 + t183;
t60 = -qJD(1) * t182 - t76 * qJD(3);
t18 = t127 * t59 - t129 * t60;
t179 = qJD(1) * qJD(3);
t173 = t132 * t179;
t111 = t127 * t173;
t172 = t134 * t179;
t87 = t129 * t172 - t111;
t13 = t87 * pkin(5) + t18;
t115 = pkin(3) * t173;
t171 = -t87 * qJ(5) + t115;
t146 = -qJD(5) * t143 + t171;
t226 = pkin(4) + pkin(8);
t14 = t226 * t86 + t146;
t67 = t127 * t76;
t75 = -t132 * t164 + t124;
t70 = qJD(3) * pkin(3) + t75;
t32 = t129 * t70 - t67;
t161 = qJD(5) - t32;
t222 = t143 * pkin(5);
t15 = -t226 * qJD(3) + t161 + t222;
t120 = -cos(pkin(9)) * pkin(1) - pkin(2);
t108 = -t134 * pkin(3) + t120;
t192 = qJD(1) * t108;
t97 = qJD(4) + t192;
t140 = -qJ(5) * t143 + t97;
t25 = t226 * t98 + t140;
t6 = t131 * t15 + t133 * t25;
t2 = -qJD(6) * t6 + t133 * t13 - t131 * t14;
t232 = t240 * t6 + t2;
t153 = t131 * t25 - t133 * t15;
t1 = -qJD(6) * t153 + t131 * t13 + t133 * t14;
t157 = t153 * t240 + t1;
t74 = t133 * qJD(3) + t131 * t98;
t45 = qJD(6) * t74 - t133 * t86;
t245 = t240 * t74 - t45;
t235 = t133 * t240;
t244 = t74 * t235;
t167 = t131 * t240;
t81 = t133 * t87;
t149 = -t167 * t240 + t81;
t200 = t129 * t134;
t106 = t127 * t132 - t200;
t152 = -t100 * t143 - t106 * t87;
t103 = qJD(3) * t200 - t127 * t188;
t233 = -t103 * t98 - t107 * t86;
t242 = -t152 + t233;
t241 = t152 + t233;
t227 = t98 ^ 2;
t96 = t143 ^ 2;
t238 = -t227 - t96;
t237 = -t227 + t96;
t61 = t87 * t107;
t236 = t103 * t143 + t61;
t36 = t129 * t75 - t67;
t195 = -qJD(5) + t36;
t234 = 0.2e1 * t143;
t119 = -t129 * pkin(3) - pkin(4);
t114 = -pkin(8) + t119;
t224 = t98 * pkin(5);
t211 = t129 * t76;
t33 = t127 * t70 + t211;
t30 = -qJD(3) * qJ(5) - t33;
t17 = -t30 - t224;
t231 = t114 * t87 + t17 * t240;
t203 = t100 * t131;
t210 = t131 * t87;
t230 = -t106 * (t186 * t240 + t210) - t240 * t203;
t187 = qJD(6) * t131;
t202 = t100 * t133;
t209 = t133 * t44;
t229 = -t106 * (t74 * t187 + t209) + t74 * t202;
t228 = qJD(3) * (-t98 + t174) - t111;
t225 = t86 * pkin(4);
t223 = pkin(3) * t132;
t196 = qJ(4) + t118;
t104 = t196 * t132;
t105 = t196 * t134;
t54 = t129 * t104 + t127 * t105;
t221 = t18 * t54;
t220 = t74 * t72;
t219 = t74 * t98;
t218 = t98 * t72;
t217 = t74 * t103 - t44 * t107;
t19 = t127 * t60 + t129 * t59;
t38 = t133 * t45;
t208 = t18 * t106;
t41 = t98 * pkin(4) + t140;
t207 = t41 * t143;
t206 = t45 * t131;
t205 = t98 * t143;
t199 = t132 * t109;
t135 = qJD(3) ^ 2;
t198 = t135 * t132;
t197 = t135 * t134;
t194 = t222 - t195;
t193 = t132 ^ 2 - t134 ^ 2;
t110 = qJD(1) * t120;
t189 = qJD(3) * t100;
t185 = t103 * qJD(3);
t178 = -t72 * t203 + (-t186 * t72 - t206) * t106;
t123 = pkin(3) * t188;
t122 = pkin(3) * t191;
t136 = qJD(1) ^ 2;
t175 = t132 * t136 * t134;
t35 = t127 * t75 + t211;
t168 = qJD(3) * t196;
t77 = -t132 * t168 + t181;
t78 = -t134 * t168 - t182;
t39 = t127 * t77 - t129 * t78;
t170 = t98 * qJ(5) + t122;
t163 = t132 * t172;
t16 = -qJD(3) * qJD(5) - t19;
t11 = -t86 * pkin(5) - t16;
t162 = -qJD(6) * t114 * t240 + t11;
t159 = t131 * t6 - t133 * t153;
t158 = t131 * t153 + t133 * t6;
t156 = t98 * t100 + t86 * t106;
t155 = -t103 * t72 - t107 * t45;
t40 = t127 * t78 + t129 * t77;
t141 = -t107 * qJ(5) + t108;
t37 = t226 * t106 + t141;
t42 = t107 * pkin(5) + t54;
t10 = t131 * t42 + t133 * t37;
t9 = -t131 * t37 + t133 * t42;
t55 = -t127 * t104 + t129 * t105;
t89 = t134 * t109 + t183;
t151 = 0.2e1 * qJD(3) * t110;
t150 = -t35 * qJD(3) + t18;
t145 = -t103 * qJ(5) - t107 * qJD(5) + t123;
t144 = t240 * t202 + (-t187 * t240 + t81) * t106;
t142 = -t235 * t240 - t210;
t56 = -t111 + (t98 + t174) * qJD(3);
t139 = qJD(6) * t158 + t1 * t131 + t2 * t133;
t138 = t18 * t107 + t143 * t39 - t40 * t98 + t54 * t87 - t55 * t86;
t83 = t89 * qJD(3);
t88 = t124 - t199;
t137 = t83 * t132 + t82 * t134 + (-t132 * t89 - t134 * t88) * qJD(3);
t116 = t127 * pkin(3) + qJ(5);
t50 = t106 * pkin(4) + t141;
t48 = pkin(4) * t143 + t170;
t43 = -t106 * pkin(5) + t55;
t34 = t100 * pkin(4) + t145;
t31 = t143 * t226 + t170;
t28 = -qJD(3) * pkin(4) + t161;
t26 = t146 + t225;
t24 = -t100 * pkin(5) + t40;
t23 = t103 * pkin(5) + t39;
t21 = t35 - t224;
t20 = t226 * t100 + t145;
t8 = t131 * t21 + t133 * t31;
t7 = -t131 * t31 + t133 * t21;
t4 = -qJD(6) * t10 - t131 * t20 + t133 * t23;
t3 = qJD(6) * t9 + t131 * t23 + t133 * t20;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t163, -0.2e1 * t193 * t179, t197, -0.2e1 * t163, -t198, 0, -t118 * t197 + t132 * t151, t118 * t198 + t134 * t151, t137, t137 * t118, t236, t241, t185, t156, -t189, 0, t97 * t100 + t108 * t86 + (-t39 + (qJD(1) * t106 + t98) * t223) * qJD(3), t97 * t103 + t108 * t87 + (t234 * t223 - t40) * qJD(3), -t33 * t100 - t32 * t103 - t19 * t106 + t138, t221 + t19 * t55 - t32 * t39 + t33 * t40 + (t97 + t192) * t123, 0, -t185, t189, t236, t241, t156, t30 * t100 + t28 * t103 + t16 * t106 + t138, t39 * qJD(3) - t41 * t100 - t26 * t106 - t34 * t98 - t50 * t86, t40 * qJD(3) - t41 * t103 - t26 * t107 - t143 * t34 - t50 * t87, -t16 * t55 + t26 * t50 + t28 * t39 - t30 * t40 + t41 * t34 + t221, t74 * t203 + (-t131 * t44 + t74 * t186) * t106, t178 + t229, t217 - t230, -t72 * t202 + (t72 * t187 - t38) * t106, t144 + t155, t103 * t240 + t61, -t17 * t202 - t153 * t103 + t2 * t107 + t24 * t72 + t4 * t240 + t43 * t45 + t9 * t87 + (-t11 * t133 + t17 * t187) * t106, t17 * t203 - t1 * t107 - t10 * t87 - t6 * t103 + t24 * t74 - t3 * t240 - t43 * t44 + (t11 * t131 + t17 * t186) * t106, -t10 * t45 - t3 * t72 - t4 * t74 + t9 * t44 + t158 * t100 + (-qJD(6) * t159 + t1 * t133 - t2 * t131) * t106, t1 * t10 + t11 * t43 - t153 * t4 + t17 * t24 + t2 * t9 + t6 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t198, -t197, 0, t82 * t132 - t83 * t134 + (-t132 * t88 + t134 * t89) * qJD(3), 0, 0, 0, 0, 0, 0, -t189, -t185, t242, -t32 * t100 + t33 * t103 + t19 * t107 + t208, 0, 0, 0, 0, 0, 0, t242, t189, t185, t28 * t100 - t30 * t103 - t16 * t107 + t208, 0, 0, 0, 0, 0, 0, t144 - t155, t217 + t230, t178 - t229, t100 * t159 + t17 * t103 + t106 * t139 + t11 * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t175, t193 * t136, 0, t175, 0, 0, -t110 * t191, -t110 * t190 - t121 + (t88 + t199) * qJD(3), 0, 0, t205, t237, t56, -t205, 0, 0, -t122 * t98 - t143 * t97 - t150, t36 * qJD(3) - t122 * t143 + t97 * t98 - t19 (-t32 + t36) * t98 + (t33 - t35) * t143 + (-t127 * t86 - t129 * t87) * pkin(3), t32 * t35 - t33 * t36 + (t127 * t19 - t129 * t18 - t97 * t191) * pkin(3), 0, -t56, 0, t205, t237, -t205, -t116 * t86 + t119 * t87 + (-t30 - t35) * t143 + (t28 + t195) * t98, t48 * t98 + t150 + t207, t48 * t143 - t41 * t98 + (0.2e1 * qJD(5) - t36) * qJD(3) + t19, -t16 * t116 + t18 * t119 + t195 * t30 - t28 * t35 - t41 * t48, -t167 * t74 - t209, -t38 - t244 + (t44 + t169) * t131, t149 + t219, t235 * t72 + t206, t142 - t218, t240 * t98, t116 * t45 + t162 * t131 + t133 * t231 - t153 * t98 + t194 * t72 - t240 * t7, -t116 * t44 - t131 * t231 + t162 * t133 + t194 * t74 + t240 * t8 - t6 * t98, t7 * t74 + t8 * t72 + (-t143 * t6 + t114 * t44 - t2 + (-t114 * t72 - t6) * qJD(6)) * t133 + (-t143 * t153 - t114 * t45 - t1 + (t114 * t74 - t153) * qJD(6)) * t131, t11 * t116 + t114 * t139 + t153 * t7 + t17 * t194 - t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t234 * qJD(3), t228, t238, t143 * t32 + t33 * t98 + t115, 0, 0, 0, 0, 0, 0, t238, -0.2e1 * t143 * qJD(3), -t228, t225 - t30 * t98 + (-qJD(5) - t28) * t143 + t171, 0, 0, 0, 0, 0, 0, t142 + t218, t219 - t149, -t243 * t131 + t244 - t38, -t131 * t232 + t157 * t133 + t17 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t205, -t96 - t135, t30 * qJD(3) + t18 + t207, 0, 0, 0, 0, 0, 0, -qJD(3) * t72 + t149, -qJD(3) * t74 + t142, t245 * t131 + t243 * t133, -t17 * qJD(3) + t157 * t131 + t133 * t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t220, -t72 ^ 2 + t74 ^ 2, -t243, -t220, t245, t87, -t17 * t74 + t232, t17 * t72 - t157, 0, 0;];
tauc_reg  = t5;
