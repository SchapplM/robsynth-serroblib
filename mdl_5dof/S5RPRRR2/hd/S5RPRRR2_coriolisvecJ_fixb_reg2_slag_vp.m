% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRR2
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
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:12:14
% EndTime: 2019-12-05 18:12:27
% DurationCPUTime: 3.30s
% Computational Cost: add. (6778->299), mult. (18466->396), div. (0->0), fcn. (14551->8), ass. (0->163)
t165 = cos(qJ(5));
t162 = sin(qJ(5));
t192 = qJD(5) * t162;
t161 = cos(pkin(9));
t167 = cos(qJ(3));
t198 = t167 * t161;
t187 = qJD(1) * t198;
t160 = sin(pkin(9));
t164 = sin(qJ(3));
t199 = t164 * t160;
t188 = qJD(1) * t199;
t132 = -t187 + t188;
t143 = t167 * t160 + t164 * t161;
t134 = t143 * qJD(1);
t163 = sin(qJ(4));
t166 = cos(qJ(4));
t100 = -t166 * t132 - t163 * t134;
t232 = t100 * pkin(8);
t215 = pkin(6) + qJ(2);
t147 = t215 * t160;
t144 = qJD(1) * t147;
t148 = t215 * t161;
t145 = qJD(1) * t148;
t109 = -t164 * t144 + t167 * t145;
t84 = -t132 * pkin(7) + t109;
t81 = t166 * t84;
t200 = t164 * t145;
t108 = -t167 * t144 - t200;
t83 = -t134 * pkin(7) + t108;
t82 = qJD(3) * pkin(3) + t83;
t35 = t163 * t82 + t81;
t29 = t35 + t232;
t137 = t143 * qJD(3);
t124 = qJD(1) * t137;
t190 = qJD(1) * qJD(2);
t184 = t160 * t190;
t151 = qJD(2) * t198;
t195 = qJD(3) * t167;
t197 = qJD(1) * t151 - t144 * t195;
t77 = (-qJD(3) * t145 - t184) * t164 + t197;
t64 = -t124 * pkin(7) + t77;
t150 = qJD(3) * t187;
t123 = qJD(3) * t188 - t150;
t173 = t143 * qJD(2);
t172 = qJD(1) * t173;
t78 = -t109 * qJD(3) - t172;
t65 = t123 * pkin(7) + t78;
t181 = -t163 * t64 + t166 * t65;
t17 = -t35 * qJD(4) + t181;
t193 = qJD(4) * t166;
t194 = qJD(4) * t163;
t174 = -t166 * t123 - t163 * t124 - t132 * t193 - t134 * t194;
t7 = -pkin(8) * t174 + t17;
t185 = t162 * t7 - t29 * t192;
t159 = qJD(3) + qJD(4);
t176 = t163 * t132 - t166 * t134;
t233 = pkin(8) * t176;
t79 = t163 * t84;
t34 = t166 * t82 - t79;
t28 = t34 + t233;
t26 = t159 * pkin(4) + t28;
t169 = t176 * qJD(4) + t163 * t123 - t166 * t124;
t180 = -t163 * t65 - t166 * t64 - t82 * t193 + t84 * t194;
t6 = pkin(8) * t169 - t180;
t1 = (qJD(5) * t26 + t6) * t165 + t185;
t52 = t165 * t100 + t162 * t176;
t154 = -t161 * pkin(2) - pkin(1);
t146 = t154 * qJD(1) + qJD(2);
t110 = t132 * pkin(3) + t146;
t66 = -pkin(4) * t100 + t110;
t231 = t66 * t52;
t240 = -t1 - t231;
t191 = qJD(5) * t165;
t18 = -t100 * t191 - t162 * t169 - t165 * t174 - t176 * t192;
t156 = qJD(5) + t159;
t208 = t52 * t156;
t228 = -t18 - t208;
t226 = t162 * t100 - t165 * t176;
t217 = t226 ^ 2;
t218 = t52 ^ 2;
t229 = t217 - t218;
t189 = -t162 * t6 + t165 * t7;
t210 = t165 * t29;
t9 = t162 * t26 + t210;
t2 = -t9 * qJD(5) + t189;
t224 = t66 * t226;
t239 = t2 - t224;
t216 = t226 * t52;
t211 = t162 * t29;
t8 = t165 * t26 - t211;
t238 = t9 * t226 + t8 * t52;
t19 = qJD(5) * t226 + t162 * t174 - t165 * t169;
t209 = t226 * t156;
t223 = -t19 + t209;
t207 = t176 * t159;
t236 = t169 - t207;
t204 = t100 * t159;
t235 = t174 - t204;
t213 = t176 ^ 2;
t214 = t100 ^ 2;
t234 = t213 - t214;
t221 = pkin(4) * t176;
t212 = t100 * t176;
t230 = t110 * t176;
t227 = -t110 * t100 + t180;
t112 = -t164 * t147 + t167 * t148;
t222 = t134 ^ 2;
t220 = t124 * pkin(3);
t219 = t134 * pkin(3);
t38 = t166 * t83 - t79;
t111 = -t167 * t147 - t164 * t148;
t95 = -t143 * pkin(7) + t111;
t142 = t198 - t199;
t96 = t142 * pkin(7) + t112;
t41 = t163 * t95 + t166 * t96;
t155 = t166 * pkin(3) + pkin(4);
t202 = t162 * t163;
t37 = -t163 * t83 - t81;
t30 = t37 - t232;
t31 = t38 + t233;
t206 = -t162 * t30 - t165 * t31 + t155 * t191 + (-t163 * t192 + (t165 * t166 - t202) * qJD(4)) * pkin(3);
t201 = t163 * t165;
t205 = t162 * t31 - t165 * t30 - t155 * t192 + (-t163 * t191 + (-t162 * t166 - t201) * qJD(4)) * pkin(3);
t203 = t134 * t132;
t196 = t160 ^ 2 + t161 ^ 2;
t186 = -pkin(4) * t156 - t26;
t40 = -t163 * t96 + t166 * t95;
t179 = t196 * qJD(1) ^ 2;
t36 = -pkin(4) * t169 + t220;
t107 = t163 * t142 + t166 * t143;
t32 = -t107 * pkin(8) + t40;
t106 = -t166 * t142 + t163 * t143;
t33 = -t106 * pkin(8) + t41;
t20 = -t162 * t33 + t165 * t32;
t21 = t162 * t32 + t165 * t33;
t63 = -t162 * t106 + t165 * t107;
t116 = -t142 * pkin(3) + t154;
t175 = 0.2e1 * t196 * t190;
t86 = -t147 * t195 + t151 + (-qJD(2) * t160 - qJD(3) * t148) * t164;
t70 = -t137 * pkin(7) + t86;
t136 = qJD(3) * t199 - t161 * t195;
t87 = -t112 * qJD(3) - t173;
t71 = t136 * pkin(7) + t87;
t24 = t163 * t71 + t166 * t70 + t95 * t193 - t96 * t194;
t25 = -t41 * qJD(4) - t163 * t70 + t166 * t71;
t129 = pkin(3) * t201 + t162 * t155;
t128 = -pkin(3) * t202 + t165 * t155;
t127 = t132 ^ 2;
t75 = t106 * pkin(4) + t116;
t72 = t219 - t221;
t62 = t165 * t106 + t162 * t107;
t61 = t107 * qJD(4) - t163 * t136 + t166 * t137;
t60 = t166 * t136 + t163 * t137 - t142 * t193 + t143 * t194;
t39 = t137 * pkin(3) + t61 * pkin(4);
t23 = t63 * qJD(5) - t162 * t60 + t165 * t61;
t22 = t106 * t191 + t107 * t192 + t162 * t61 + t165 * t60;
t15 = t60 * pkin(8) + t25;
t14 = -t61 * pkin(8) + t24;
t11 = t165 * t28 - t211;
t10 = -t162 * t28 - t210;
t4 = -t21 * qJD(5) - t162 * t14 + t165 * t15;
t3 = t20 * qJD(5) + t165 * t14 + t162 * t15;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175, qJ(2) * t175, -t123 * t143 - t134 * t136, -t123 * t142 - t143 * t124 + t136 * t132 - t134 * t137, -t136 * qJD(3), -t124 * t142 + t132 * t137, -t137 * qJD(3), 0, t87 * qJD(3) + t154 * t124 + t146 * t137, -t86 * qJD(3) - t154 * t123 - t146 * t136, t108 * t136 - t109 * t137 + t111 * t123 - t112 * t124 - t86 * t132 - t87 * t134 + t77 * t142 - t78 * t143, t108 * t87 + t109 * t86 + t78 * t111 + t77 * t112, t107 * t174 + t176 * t60, -t100 * t60 - t106 * t174 + t107 * t169 + t176 * t61, -t60 * t159, -t100 * t61 - t106 * t169, -t61 * t159, 0, t110 * t61 - t116 * t169 + t25 * t159 + (-t100 * t137 + t106 * t124) * pkin(3), -t110 * t60 + t116 * t174 - t24 * t159 + (t107 * t124 - t137 * t176) * pkin(3), t100 * t24 + t106 * t180 - t17 * t107 + t169 * t41 - t174 * t40 + t176 * t25 + t34 * t60 - t35 * t61, -t180 * t41 + t17 * t40 + t35 * t24 + t34 * t25 + (t110 * t137 + t116 * t124) * pkin(3), -t18 * t63 - t22 * t226, t18 * t62 - t63 * t19 - t22 * t52 - t226 * t23, -t22 * t156, t19 * t62 - t23 * t52, -t23 * t156, 0, t4 * t156 + t75 * t19 + t66 * t23 + t36 * t62 - t39 * t52, -t3 * t156 - t75 * t18 - t66 * t22 + t226 * t39 + t36 * t63, -t1 * t62 + t20 * t18 - t21 * t19 - t2 * t63 + t8 * t22 - t226 * t4 - t9 * t23 + t3 * t52, t1 * t21 + t2 * t20 + t9 * t3 + t36 * t75 + t66 * t39 + t8 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t179, -qJ(2) * t179, 0, 0, 0, 0, 0, 0, 0.2e1 * t134 * qJD(3), t150 + (-t132 - t188) * qJD(3), -t127 - t222, t108 * t134 + t109 * t132, 0, 0, 0, 0, 0, 0, -t169 - t207, t174 + t204, -t213 - t214, -t35 * t100 - t176 * t34 + t220, 0, 0, 0, 0, 0, 0, t19 + t209, -t18 + t208, -t217 - t218, t226 * t8 - t9 * t52 + t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t203, -t127 + t222, t150 + (t132 - t188) * qJD(3), -t203, 0, 0, -t146 * t134 - t172, t164 * t184 + t146 * t132 + (t108 + t200) * qJD(3) - t197, 0, 0, t212, t234, t235, -t212, t236, 0, t100 * t219 + t230 - t37 * t159 + (-t81 + (-pkin(3) * t159 - t82) * t163) * qJD(4) + t181, t38 * t159 + (t134 * t176 - t159 * t193) * pkin(3) + t227, t34 * t100 - t37 * t176 - t35 * t176 - t38 * t100 + (t163 * t169 - t166 * t174 + (t100 * t166 - t163 * t176) * qJD(4)) * pkin(3), -t34 * t37 - t35 * t38 + (-t110 * t134 - t180 * t163 + t166 * t17 + (-t163 * t34 + t166 * t35) * qJD(4)) * pkin(3), -t216, t229, t228, t216, t223, 0, t205 * t156 + t52 * t72 + t239, -t206 * t156 - t226 * t72 + t240, t128 * t18 - t129 * t19 - t205 * t226 + t206 * t52 + t238, t1 * t129 + t2 * t128 + t205 * t8 + t206 * t9 - t66 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t212, t234, t235, -t212, t236, 0, t35 * t159 + t17 + t230, t34 * t159 + t227, 0, 0, -t216, t229, t228, t216, t223, 0, -t52 * t221 - t10 * t156 - t224 + (t162 * t186 - t210) * qJD(5) + t189, t226 * t221 + t11 * t156 - t231 + (qJD(5) * t186 - t6) * t165 - t185, t10 * t226 - t11 * t52 + (-t162 * t19 + t165 * t18 + (t162 * t226 + t165 * t52) * qJD(5)) * pkin(4) + t238, -t8 * t10 - t9 * t11 + (t1 * t162 + t176 * t66 + t165 * t2 + (-t162 * t8 + t165 * t9) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t216, t229, t228, t216, t223, 0, t9 * t156 + t239, t8 * t156 + t240, 0, 0;];
tauc_reg = t5;
