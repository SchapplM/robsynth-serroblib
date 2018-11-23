% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 14:53
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 14:52:10
% EndTime: 2018-11-16 14:52:19
% DurationCPUTime: 3.79s
% Computational Cost: add. (5912->326), mult. (15303->478), div. (0->0), fcn. (12213->8), ass. (0->174)
t118 = cos(qJ(5));
t190 = qJD(5) * t118;
t115 = sin(qJ(4));
t226 = cos(qJ(4));
t116 = sin(qJ(3));
t119 = cos(qJ(2));
t227 = cos(qJ(3));
t176 = t227 * t119;
t117 = sin(qJ(2));
t194 = qJD(1) * t117;
t92 = -qJD(1) * t176 + t116 * t194;
t94 = t116 * t119 + t227 * t117;
t93 = qJD(1) * t94;
t61 = t115 * t93 + t226 * t92;
t237 = t118 * t61;
t242 = t190 - t237;
t114 = sin(qJ(5));
t139 = t94 * qJD(3);
t127 = t94 * qJD(2) + t139;
t126 = t127 * qJD(1);
t171 = t226 * qJD(4);
t192 = qJD(4) * t115;
t146 = t116 * t117 - t176;
t141 = t146 * qJD(2);
t129 = t146 * qJD(3) + t141;
t66 = t129 * qJD(1);
t147 = t115 * t126 + t92 * t171 + t93 * t192 + t226 * t66;
t148 = -t115 * t92 + t226 * t93;
t110 = qJD(2) + qJD(3);
t168 = qJD(4) + t110;
t149 = t118 * t168;
t191 = qJD(5) * t114;
t18 = -qJD(5) * t149 - t118 * t147 - t148 * t191;
t47 = t114 * t168 - t118 * t148;
t201 = qJD(5) * t47;
t19 = t114 * t147 + t201;
t235 = qJD(5) - t61;
t241 = t114 * t235;
t45 = -t114 * t148 - t149;
t1 = -t114 * t19 - t18 * t118 - t241 * t47 - t242 * t45;
t124 = t226 * t127;
t169 = -qJD(1) * t124 + t115 * t66;
t31 = -t148 * qJD(4) + t169;
t27 = t118 * t31;
t5 = -t148 * t45 - t235 * t241 + t27;
t107 = t119 * pkin(2) + pkin(1);
t99 = qJD(1) * t107;
t76 = -pkin(3) * t92 + t99;
t34 = -pkin(4) * t61 + pkin(6) * t148 + t76;
t175 = t226 * t116;
t103 = pkin(2) * t175;
t173 = t227 * qJD(2);
t161 = pkin(2) * t173;
t95 = t110 * pkin(3) + t161;
t78 = qJD(2) * t103 + t115 * t95;
t75 = t168 * pkin(6) + t78;
t150 = t114 * t75 - t118 * t34;
t240 = t235 * t150;
t24 = t114 * t34 + t118 * t75;
t63 = -pkin(3) * t139 + (-t117 * pkin(2) - t94 * pkin(3)) * qJD(2);
t52 = t63 * qJD(1);
t10 = t31 * pkin(4) - pkin(6) * t147 + t52;
t158 = t226 * t227;
t210 = qJD(2) * pkin(2);
t182 = qJD(3) * t210;
t160 = t116 * t182;
t185 = t116 * t210;
t163 = t115 * t185;
t133 = qJD(4) * t163 + t115 * t160 - t158 * t182 - t95 * t171;
t4 = -qJD(5) * t24 + t118 * t10 + t114 * t133;
t239 = -t235 * t24 - t4;
t15 = t18 * t114;
t7 = t242 * t47 - t15;
t213 = t114 * t31 + t190 * t235;
t6 = t148 * t47 - t235 * t237 + t213;
t77 = t226 * t95 - t163;
t74 = -t168 * pkin(4) - t77;
t219 = t61 * t74;
t221 = t235 * t148;
t236 = t61 * t148;
t22 = t148 ^ 2 - t61 ^ 2;
t130 = -t76 * t61 + t133;
t181 = t95 * t192;
t144 = t227 * t115 + t175;
t229 = t144 * qJD(3) + t116 * t171;
t122 = t148 * t76 - t229 * t210 - t181;
t20 = -t61 * t168 + t147;
t174 = -t148 * t150 + t74 * t191;
t132 = t229 * pkin(2);
t51 = qJD(2) * t132 + t181;
t159 = t51 * t114 - t24 * t148 + t74 * t190;
t21 = -t148 * t110 - t169;
t40 = -pkin(4) * t148 - pkin(6) * t61;
t189 = qJD(1) * qJD(2);
t232 = -0.2e1 * t189;
t231 = -0.2e1 * t99;
t230 = t114 * t150 + t118 * t24;
t3 = -t150 * qJD(5) + t10 * t114 - t118 * t133;
t228 = pkin(3) * t93;
t2 = t3 * t118;
t136 = t226 * t146;
t72 = -t115 * t94 - t136;
t225 = t31 * t72;
t142 = t115 * t146;
t73 = -t226 * t94 + t142;
t224 = t31 * t73;
t223 = t47 * t45;
t222 = t51 * t73;
t216 = t92 * t93;
t215 = t93 * t99;
t214 = t99 * t92;
t208 = t114 * t24;
t206 = t114 * t45;
t17 = t19 * t118;
t200 = qJD(5) * t235;
t199 = t115 * t116;
t121 = qJD(1) ^ 2;
t198 = t119 * t121;
t186 = t227 * pkin(2);
t106 = t186 + pkin(3);
t90 = t115 * t106 + t103;
t196 = t117 ^ 2 - t119 ^ 2;
t193 = qJD(3) * t116;
t188 = -t150 * t237 + t61 * t208 + t2;
t184 = t117 * t210;
t183 = pkin(2) * t194;
t180 = t73 * t191;
t179 = t73 * t190;
t178 = t117 * t198;
t177 = t227 * t110;
t172 = t227 * qJD(3);
t170 = t117 * t189;
t38 = -t228 + t40;
t35 = t38 - t183;
t85 = pkin(6) + t90;
t167 = qJD(5) * t85 + t35;
t164 = pkin(1) * t232;
t162 = pkin(3) * t171;
t156 = t119 * t170;
t87 = t144 * t210;
t155 = pkin(3) * t192 - t87;
t36 = -qJD(4) * t136 - t115 * t127 - t226 * t129 - t94 * t192;
t154 = -t235 * t36 + t224;
t153 = -t148 * t78 + t77 * t61;
t151 = -t118 * t150 + t208;
t143 = t158 - t199;
t64 = t106 * t171 + (t143 * qJD(3) - t116 * t192) * pkin(2);
t145 = -t235 * t64 - t31 * t85 - t219;
t89 = -pkin(2) * t199 + t226 * t106;
t37 = qJD(4) * t142 + t115 * t129 - t94 * t171 - t124;
t11 = t37 * pkin(4) + t36 * pkin(6) + t63;
t81 = -t146 * pkin(3) + t107;
t39 = t72 * pkin(4) - t73 * pkin(6) + t81;
t138 = qJD(5) * t73 * t74 + t11 * t235 + t31 * t39;
t137 = -t39 * t200 - t36 * t74 + t222;
t104 = t115 * pkin(3) + pkin(6);
t135 = -t104 * t31 - t162 * t235 - t219;
t134 = -t151 * qJD(5) - t4 * t114 + t2;
t125 = pkin(2) * t116 * t127;
t120 = qJD(2) ^ 2;
t105 = -t226 * pkin(3) - pkin(4);
t88 = t143 * t210;
t84 = -pkin(4) - t89;
t79 = -t183 - t228;
t65 = t106 * t192 + t132;
t49 = -t92 ^ 2 + t93 ^ 2;
t42 = -t93 * t110 + t126;
t41 = -t92 * t110 + t66;
t33 = t114 * t40 + t118 * t77;
t32 = -t114 * t77 + t118 * t40;
t29 = t114 * t38 + t118 * t88;
t28 = -t114 * t88 + t118 * t38;
t8 = t241 * t45 - t17;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t156, t196 * t232, -t120 * t119, -0.2e1 * t156, t120 * t117, 0, t117 * t164, t119 * t164, 0, 0, -t93 * t129 - t66 * t94, -t94 * t126 + t66 * t146 + t110 * (t146 * t92 - t93 * t94) t129 * t110, t146 * t126 + t92 * t127, t127 * t110, 0, t231 * t110 * t94 + t141 * t183 + t92 * t184, pkin(2) * t94 * t170 + t107 * t66 + t99 * t129 + t93 * t184, qJD(2) * t125 - t141 * t161 - t94 * t160, t231 * t184, t147 * t73 + t148 * t36, -t147 * t72 + t148 * t37 - t36 * t61 - t224, -t36 * t168, -t37 * t61 + t225, -t37 * t168, 0, t31 * t81 + t37 * t76 + t52 * t72 - t61 * t63, t147 * t81 - t148 * t63 - t36 * t76 + t52 * t73, t133 * t72 + t36 * t77 - t37 * t78 + t222, t52 * t81 + t63 * t76, -t47 * t180 + (-t18 * t73 - t36 * t47) * t118 (t114 * t47 + t118 * t45) * t36 + (t15 - t17 + (-t118 * t47 + t206) * qJD(5)) * t73, t118 * t154 - t18 * t72 - t180 * t235 + t37 * t47, t45 * t179 + (t19 * t73 - t36 * t45) * t114, -t114 * t154 - t179 * t235 - t19 * t72 - t37 * t45, t235 * t37 + t225, t114 * t137 + t118 * t138 - t150 * t37 + t4 * t72, -t114 * t138 + t118 * t137 - t24 * t37 - t3 * t72 (-t11 * t47 + t18 * t39 - t150 * t36 - t4 * t73 + (-t24 * t73 - t39 * t45) * qJD(5)) * t118 + (-t11 * t45 - t19 * t39 + t24 * t36 - t3 * t73 + (-t150 * t73 + t39 * t47) * qJD(5)) * t114, t151 * t11 + (t230 * qJD(5) + t114 * t3 + t118 * t4) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t178, t196 * t121, 0, t178, 0, 0, t121 * pkin(1) * t117, pkin(1) * t198, 0, 0, t216, t49, t41, -t216, t42, 0, t215 + (-t92 * t194 + (-qJD(2) - t110) * t193) * pkin(2), -t214 + (-t93 * t194 + (-t173 - t177) * qJD(3)) * pkin(2), qJD(1) * t125 + t92 * t161 - t93 * t185 - t66 * t186 + (t172 * t92 - t193 * t93) * pkin(2), t99 * t183, t236, t22, t20, -t236, t21, 0, -t65 * t168 + t61 * t79 + t122, t148 * t79 - t168 * t64 + t130, -t147 * t89 - t148 * t65 - t31 * t90 + t61 * t64 + t153, -t133 * t90 - t51 * t89 + t64 * t78 - t65 * t77 - t76 * t79, t7, t1, t6, t8, t5, t221, t19 * t84 + t45 * t65 + (-t167 * t235 - t51) * t118 + t145 * t114 + t174, t118 * t145 + t167 * t241 - t18 * t84 + t47 * t65 + t159 (-t19 * t85 + t35 * t47 - t45 * t64 + (t47 * t85 + t150) * qJD(5)) * t118 + (-t18 * t85 + t35 * t45 + t47 * t64 - t4 + (t45 * t85 - t24) * qJD(5)) * t114 + t188, t51 * t84 + t65 * t74 + (t150 * t167 + t24 * t64 + t3 * t85) * t118 + (t150 * t64 - t167 * t24 - t4 * t85) * t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t216, t49, t41, -t216, t42, 0, t215 + (-qJD(3) + t110) * t185, -t214 + (-t172 + t177) * t210, 0, 0, t236, t22, t20, -t236, t21, 0, t87 * t168 + (-t168 * t192 - t61 * t93) * pkin(3) + t122, t88 * t168 + (-t148 * t93 - t168 * t171) * pkin(3) + t130, -t88 * t61 + t87 * t148 + (-t226 * t147 - t115 * t31 + (-t115 * t148 + t226 * t61) * qJD(4)) * pkin(3) + t153, t77 * t87 - t78 * t88 + (-t226 * t51 - t115 * t133 + t76 * t93 + (-t115 * t77 + t226 * t78) * qJD(4)) * pkin(3), t7, t1, t6, t8, t5, t221, t105 * t19 - t28 * t235 + t155 * t45 + (-t104 * t200 - t51) * t118 + t135 * t114 + t174, -t105 * t18 + (t104 * t191 + t29) * t235 + t155 * t47 + t135 * t118 + t159, t28 * t47 + t29 * t45 + (-t45 * t162 - t104 * t19 + (t104 * t47 + t150) * qJD(5)) * t118 + (t47 * t162 - t104 * t18 - t4 + (t104 * t45 - t24) * qJD(5)) * t114 + t188, t51 * t105 + t150 * t28 - t24 * t29 - t74 * t87 + (t115 * t74 + t230 * t226) * qJD(4) * pkin(3) + t134 * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t236, t22, t20, -t236, t21, 0, t168 * t78 + t122, t168 * t77 + t130, 0, 0, t7, t1, t6, t206 * t235 - t17, t5, t221, -pkin(4) * t19 - t213 * pkin(6) - t114 * t219 - t51 * t118 - t235 * t32 - t45 * t78 + t174, -t74 * t237 + pkin(4) * t18 + t33 * t235 - t47 * t78 + (t191 * t235 - t27) * pkin(6) + t159, t32 * t47 + t33 * t45 + t2 + (t240 + (-t19 + t201) * pkin(6)) * t118 + ((qJD(5) * t45 - t18) * pkin(6) + t239) * t114, -pkin(4) * t51 + pkin(6) * t134 + t150 * t32 - t24 * t33 - t74 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t223, -t45 ^ 2 + t47 ^ 2, t235 * t45 - t18, -t223, t235 * t47 - t19, t31, -t47 * t74 - t239, t45 * t74 - t240 - t3, 0, 0;];
tauc_reg  = t9;
