% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x31]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:22:27
% EndTime: 2019-12-31 22:22:35
% DurationCPUTime: 2.81s
% Computational Cost: add. (4538->286), mult. (11678->396), div. (0->0), fcn. (8834->8), ass. (0->179)
t146 = cos(qJ(5));
t196 = qJD(5) * t146;
t148 = cos(qJ(3));
t149 = cos(qJ(2));
t201 = qJD(1) * t149;
t190 = t148 * t201;
t144 = sin(qJ(3));
t145 = sin(qJ(2));
t202 = qJD(1) * t145;
t191 = t144 * t202;
t102 = t190 - t191;
t103 = -t144 * t201 - t148 * t202;
t143 = sin(qJ(4));
t147 = cos(qJ(4));
t77 = t147 * t102 + t143 * t103;
t248 = t146 * t77;
t251 = t196 - t248;
t135 = -t149 * pkin(2) - pkin(1);
t121 = t135 * qJD(1);
t89 = -t102 * pkin(3) + t121;
t226 = t89 * t77;
t199 = qJD(4) * t143;
t237 = pkin(6) + pkin(7);
t123 = t237 * t149;
t118 = qJD(1) * t123;
t108 = t148 * t118;
t122 = t237 * t145;
t116 = qJD(1) * t122;
t221 = qJD(2) * pkin(2);
t110 = -t116 + t221;
t163 = -t144 * t110 - t108;
t192 = qJD(2) * t237;
t168 = qJD(1) * t192;
t111 = t145 * t168;
t112 = t149 * t168;
t178 = t144 * t111 - t148 * t112;
t154 = t163 * qJD(3) + t178;
t139 = qJD(2) + qJD(3);
t195 = qJD(1) * qJD(2);
t187 = t149 * t195;
t82 = qJD(3) * t190 - t139 * t191 + t148 * t187;
t36 = -t82 * pkin(8) + t154;
t235 = t102 * pkin(8);
t59 = -t163 + t235;
t183 = t143 * t36 - t59 * t199;
t200 = qJD(3) * t144;
t172 = -t144 * t112 - t118 * t200;
t242 = (qJD(3) * t110 - t111) * t148;
t115 = t144 * t149 + t148 * t145;
t88 = t139 * t115;
t83 = t88 * qJD(1);
t35 = -t83 * pkin(8) + t172 + t242;
t104 = t144 * t118;
t177 = t148 * t110 - t104;
t99 = t103 * pkin(8);
t58 = t177 + t99;
t52 = t139 * pkin(3) + t58;
t3 = (qJD(4) * t52 + t35) * t147 + t183;
t250 = -t226 - t3;
t205 = -qJD(5) + t77;
t249 = qJD(5) + t205;
t142 = sin(qJ(5));
t138 = qJD(4) + t139;
t164 = t143 * t102 - t147 * t103;
t197 = qJD(5) * t142;
t198 = qJD(4) * t147;
t30 = t102 * t198 + t103 * t199 - t143 * t83 + t147 * t82;
t18 = t138 * t196 + t146 * t30 - t164 * t197;
t16 = t18 * t142;
t62 = t142 * t138 + t146 * t164;
t7 = t251 * t62 + t16;
t31 = t164 * qJD(4) + t143 * t82 + t147 * t83;
t67 = t205 * t196;
t224 = t142 * t31 - t67;
t6 = -t164 * t62 + t205 * t248 + t224;
t247 = t142 * t205;
t27 = t146 * t31;
t60 = -t146 * t138 + t142 * t164;
t5 = t164 * t60 - t205 * t247 + t27;
t19 = qJD(5) * t62 + t142 * t30;
t1 = -t142 * t19 + t18 * t146 + t247 * t62 - t251 * t60;
t218 = t143 * t59;
t32 = t147 * t52 - t218;
t28 = -t138 * pkin(4) - t32;
t233 = t28 * t77;
t229 = t164 * t77;
t227 = t89 * t164;
t184 = t143 * t35 - t147 * t36;
t215 = t147 * t59;
t33 = t143 * t52 + t215;
t4 = t33 * qJD(4) + t184;
t246 = -t4 - t227;
t22 = t164 ^ 2 - t77 ^ 2;
t51 = pkin(4) * t164 - t77 * pkin(9);
t20 = -t77 * t138 + t30;
t244 = -0.2e1 * t195;
t230 = t205 * t164;
t241 = qJD(1) * t115;
t29 = t138 * pkin(9) + t33;
t39 = -pkin(4) * t77 - pkin(9) * t164 + t89;
t165 = t142 * t29 - t146 * t39;
t188 = t164 * t165 + t28 * t197;
t14 = t142 * t39 + t146 * t29;
t169 = t14 * t164 + t4 * t142 + t28 * t196;
t21 = t138 * t164 - t31;
t117 = t145 * t192;
t119 = t149 * t192;
t209 = t148 * t122;
t159 = -qJD(3) * t209 - t148 * t117 - t144 * t119 - t123 * t200;
t46 = -t88 * pkin(8) + t159;
t162 = t144 * t122 - t148 * t123;
t152 = t162 * qJD(3) + t144 * t117 - t148 * t119;
t114 = t144 * t145 - t148 * t149;
t87 = t139 * t114;
t47 = t87 * pkin(8) + t152;
t68 = -t115 * pkin(8) - t144 * t123 - t209;
t69 = -t114 * pkin(8) - t162;
t49 = t143 * t69 - t147 * t68;
t10 = -t49 * qJD(4) + t143 * t47 + t147 * t46;
t85 = t147 * t114 + t143 * t115;
t43 = -t85 * qJD(4) - t143 * t88 - t147 * t87;
t86 = -t143 * t114 + t147 * t115;
t92 = t114 * pkin(3) + t135;
t48 = t85 * pkin(4) - t86 * pkin(9) + t92;
t50 = t143 * t68 + t147 * t69;
t238 = (qJD(5) * t48 + t10) * t205 - (qJD(5) * t39 + t3) * t85 + t28 * t43 - t50 * t31 + t4 * t86;
t234 = t103 * pkin(3);
t232 = t28 * t86;
t231 = t48 * t31;
t228 = t86 * t31;
t134 = t148 * pkin(2) + pkin(3);
t211 = t143 * t144;
t171 = t144 * t116 - t108;
t64 = t171 - t235;
t204 = -t148 * t116 - t104;
t65 = t99 + t204;
t223 = t143 * t64 + t147 * t65 - t134 * t198 - (-t144 * t199 + (t147 * t148 - t211) * qJD(3)) * pkin(2);
t210 = t144 * t147;
t222 = -t143 * t65 + t147 * t64 + t134 * t199 + (t144 * t198 + (t143 * t148 + t210) * qJD(3)) * pkin(2);
t213 = t103 * t102;
t212 = t121 * t103;
t151 = qJD(1) ^ 2;
t208 = t149 * t151;
t150 = qJD(2) ^ 2;
t207 = t150 * t145;
t206 = t150 * t149;
t203 = t145 ^ 2 - t149 ^ 2;
t137 = t145 * t221;
t136 = pkin(2) * t202;
t193 = t86 * t197;
t66 = t83 * pkin(3) + qJD(2) * t136;
t79 = t88 * pkin(3) + t137;
t189 = -pkin(3) * t138 - t52;
t186 = -pkin(2) * t139 - t110;
t45 = -t234 + t51;
t98 = pkin(2) * t210 + t143 * t134 + pkin(9);
t179 = qJD(5) * t98 + t136 + t45;
t174 = pkin(1) * t244;
t132 = t143 * pkin(3) + pkin(9);
t173 = qJD(5) * t132 + t45;
t37 = t143 * t58 + t215;
t167 = pkin(3) * t199 - t37;
t166 = -t205 * t43 + t228;
t161 = -t121 * t102 - t172;
t158 = -t205 * t223 - t98 * t31 - t233;
t38 = t147 * t58 - t218;
t153 = -t132 * t31 - t233 - (-pkin(3) * t198 + t38) * t205;
t133 = -t147 * pkin(3) - pkin(4);
t97 = pkin(2) * t211 - t147 * t134 - pkin(4);
t90 = t136 - t234;
t63 = -t102 ^ 2 + t103 ^ 2;
t55 = (-t103 - t241) * t139;
t54 = -t102 * t139 + t82;
t44 = t86 * qJD(4) - t143 * t87 + t147 * t88;
t12 = t44 * pkin(4) - t43 * pkin(9) + t79;
t11 = t50 * qJD(4) + t143 * t46 - t147 * t47;
t9 = t31 * pkin(4) - t30 * pkin(9) + t66;
t8 = t146 * t9;
t2 = [0, 0, 0, 0.2e1 * t145 * t187, t203 * t244, t206, -t207, 0, -pkin(6) * t206 + t145 * t174, pkin(6) * t207 + t149 * t174, t103 * t87 + t82 * t115, -t87 * t102 + t103 * t88 - t82 * t114 - t115 * t83, -t87 * t139, -t88 * t139, 0, t135 * t83 + t121 * t88 + t152 * t139 + (qJD(1) * t114 - t102) * t137, t135 * t82 - t121 * t87 - t159 * t139 + (-t103 + t241) * t137, t164 * t43 + t30 * t86, -t164 * t44 - t30 * t85 + t43 * t77 - t228, t43 * t138, -t44 * t138, 0, -t11 * t138 + t92 * t31 + t89 * t44 + t66 * t85 - t77 * t79, -t10 * t138 + t164 * t79 + t92 * t30 + t89 * t43 + t66 * t86, -t62 * t193 + (t18 * t86 + t43 * t62) * t146, (-t142 * t62 - t146 * t60) * t43 + (-t16 - t146 * t19 + (t142 * t60 - t146 * t62) * qJD(5)) * t86, t146 * t166 + t18 * t85 + t193 * t205 + t62 * t44, -t142 * t166 - t19 * t85 - t60 * t44 + t67 * t86, -t205 * t44 + t31 * t85, t11 * t60 - t165 * t44 + t49 * t19 + t8 * t85 + (-t12 * t205 + t231 + (t205 * t50 - t29 * t85 + t232) * qJD(5)) * t146 + t238 * t142, t11 * t62 - t14 * t44 + t49 * t18 + ((-qJD(5) * t50 + t12) * t205 - t231 - (-qJD(5) * t29 + t9) * t85 - qJD(5) * t232) * t142 + t238 * t146; 0, 0, 0, -t145 * t208, t203 * t151, 0, 0, 0, t151 * pkin(1) * t145, pkin(1) * t208, t213, t63, t54, t55, 0, t102 * t136 + t212 - t171 * t139 + (t186 * t144 - t108) * qJD(3) + t178, t103 * t136 + t204 * t139 + (t186 * qJD(3) + t111) * t148 + t161, -t229, t22, t20, t21, 0, -t138 * t222 + t77 * t90 + t246, t138 * t223 - t164 * t90 + t250, t7, t1, t6, t5, t230, t97 * t19 + t222 * t60 + (t179 * t205 - t4) * t146 + t158 * t142 + t188, t146 * t158 - t179 * t247 + t97 * t18 + t222 * t62 + t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, t63, t54, t55, 0, -t139 * t163 + t154 + t212, t177 * t139 + t161 - t242, -t229, t22, t20, t21, 0, -t77 * t234 + t37 * t138 - t227 + (t143 * t189 - t215) * qJD(4) - t184, t164 * t234 + t38 * t138 - t226 + (qJD(4) * t189 - t35) * t147 - t183, t7, t1, t6, t5, t230, t133 * t19 + t167 * t60 + (t173 * t205 - t4) * t146 + t153 * t142 + t188, t133 * t18 + t146 * t153 + t167 * t62 - t173 * t247 + t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t229, t22, t20, t21, 0, t33 * t138 + t246, t32 * t138 + t250, t7, t1, t6, t5, t230, -pkin(4) * t19 - t4 * t146 + (-t142 * t32 + t146 * t51) * t205 - t33 * t60 - t142 * t233 - t224 * pkin(9) + t188, -pkin(4) * t18 - (t142 * t51 + t146 * t32) * t205 - t33 * t62 - t28 * t248 + (-t197 * t205 - t27) * pkin(9) + t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t60, -t60 ^ 2 + t62 ^ 2, -t205 * t60 + t18, -t205 * t62 - t19, t31, -t249 * t14 - t142 * t3 - t28 * t62 + t8, -t142 * t9 - t146 * t3 + t249 * t165 + t28 * t60;];
tauc_reg = t2;
