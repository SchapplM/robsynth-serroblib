% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tau_reg [5x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:14:48
% EndTime: 2021-01-15 22:14:56
% DurationCPUTime: 1.77s
% Computational Cost: add. (2914->301), mult. (4369->346), div. (0->0), fcn. (2695->12), ass. (0->177)
t244 = qJ(4) + pkin(7);
t140 = sin(pkin(8));
t141 = cos(pkin(8));
t143 = sin(qJ(3));
t146 = cos(qJ(3));
t185 = qJD(3) * t244;
t164 = -t143 * qJD(4) - t146 * t185;
t147 = cos(qJ(2));
t203 = qJD(1) * t147;
t193 = pkin(1) * t203;
t127 = t146 * qJD(4);
t80 = -t143 * t185 + t127;
t212 = t141 * t146;
t213 = t140 * t143;
t86 = -t212 + t213;
t223 = t140 * t164 + t141 * t80 + t86 * t193;
t139 = qJ(1) + qJ(2);
t129 = cos(t139);
t116 = g(2) * t129;
t144 = sin(qJ(2));
t228 = t144 * pkin(1);
t194 = qJD(1) * t228;
t225 = t147 * pkin(1);
t205 = -qJD(2) * t194 + qJDD(1) * t225;
t133 = qJDD(1) + qJDD(2);
t230 = t133 * pkin(2);
t243 = -t205 - t230 + t116;
t226 = t146 * pkin(3);
t120 = pkin(2) + t226;
t134 = qJD(1) + qJD(2);
t71 = -t120 * t134 + qJD(4) - t193;
t192 = t134 * t212;
t72 = t134 * t213 - t192;
t87 = t140 * t146 + t141 * t143;
t74 = t87 * t134;
t30 = t72 * pkin(4) - t74 * qJ(5) + t71;
t210 = t146 * t133;
t211 = t143 * t133;
t173 = t140 * t211 - t141 * t210;
t81 = t87 * qJD(3);
t43 = t134 * t81 + t173;
t200 = qJD(3) * t143;
t189 = t140 * t200;
t162 = t87 * t133 - t134 * t189;
t199 = qJD(3) * t146;
t188 = t141 * t199;
t44 = t134 * t188 + t162;
t190 = t134 * t200;
t48 = pkin(3) * t190 - t120 * t133 + qJDD(4) - t205;
t151 = t43 * pkin(4) - t44 * qJ(5) + t48;
t7 = -t74 * qJD(5) + t151;
t242 = t30 * t81 + t7 * t86;
t82 = t188 - t189;
t241 = -t30 * t82 - t7 * t87;
t240 = t48 * t86 + t71 * t81;
t239 = t48 * t87 + t71 * t82;
t69 = t74 ^ 2;
t238 = -t72 ^ 2 - t69;
t128 = sin(t139);
t237 = g(1) * t129 + g(2) * t128;
t117 = g(1) * t128;
t236 = t117 - t116;
t184 = t244 * t134 + t194;
t64 = t184 * t146;
t221 = t140 * t64;
t63 = t184 * t143;
t62 = qJD(3) * pkin(3) - t63;
t35 = t141 * t62 - t221;
t33 = -qJD(3) * pkin(4) + qJD(5) - t35;
t58 = t141 * t64;
t36 = t140 * t62 + t58;
t34 = qJD(3) * qJ(5) + t36;
t198 = qJDD(1) * t144;
t201 = qJD(2) * t147;
t77 = t133 * pkin(7) + (qJD(1) * t201 + t198) * pkin(1);
t168 = qJ(4) * t133 + qJD(4) * t134 + t77;
t171 = qJD(3) * t184;
t25 = qJDD(3) * pkin(3) - t168 * t143 - t146 * t171;
t29 = -t143 * t171 + t168 * t146;
t11 = t140 * t25 + t141 * t29;
t136 = qJDD(3) * qJ(5);
t8 = qJD(3) * qJD(5) + t11 + t136;
t10 = -t140 * t29 + t141 * t25;
t219 = qJDD(3) * pkin(4);
t9 = qJDD(5) - t219 - t10;
t235 = t33 * t82 - t34 * t81 - t8 * t86 + t9 * t87;
t234 = -t10 * t87 - t11 * t86 - t35 * t82 - t36 * t81;
t135 = qJ(3) + pkin(8);
t125 = sin(t135);
t126 = cos(t135);
t233 = g(3) * t125 + t126 * t237 - t11;
t232 = pkin(3) * t143;
t231 = g(3) * t146;
t229 = t134 * pkin(2);
t145 = sin(qJ(1));
t227 = t145 * pkin(1);
t224 = t140 * t80 - t141 * t164 - t87 * t193;
t217 = t125 * t129;
t218 = t125 * t128;
t222 = g(1) * t218 - g(2) * t217;
t91 = -t193 - t229;
t220 = t146 * t117 + t91 * t200;
t216 = t126 * t129;
t215 = t129 * t244;
t214 = t134 * t143;
t38 = -t140 * t63 + t58;
t209 = t38 * qJD(3);
t119 = pkin(7) + t228;
t208 = -qJ(4) - t119;
t39 = -t141 * t63 - t221;
t207 = qJD(5) - t39;
t137 = t143 ^ 2;
t204 = -t146 ^ 2 + t137;
t202 = qJD(2) * t144;
t100 = t129 * t120;
t197 = pkin(4) * t216 + qJ(5) * t217 + t100;
t196 = t243 * t143 + t91 * t199;
t195 = pkin(1) * t201;
t123 = pkin(3) * t200;
t191 = t134 * t202;
t186 = t244 * t143;
t183 = t208 * t143;
t182 = t128 * t244 + t100;
t181 = qJD(3) * t208;
t180 = t134 * t194;
t179 = t205 + t236;
t177 = -g(2) * t216 + t126 * t117;
t37 = t81 * pkin(4) - t82 * qJ(5) - t87 * qJD(5) + t123;
t176 = -t37 + t194;
t130 = t146 * qJ(4);
t106 = t146 * pkin(7) + t130;
t61 = t141 * t106 - t140 * t186;
t175 = t61 * qJDD(3) + t222;
t172 = -t126 * pkin(4) - t125 * qJ(5);
t170 = -t128 * t120 + t215;
t169 = g(1) * t217 + g(2) * t218 - g(3) * t126 + t10;
t60 = t140 * t106 + t141 * t186;
t167 = -t60 * qJDD(3) + t177;
t153 = (-qJD(4) - t195) * t143 + t146 * t181;
t55 = t143 * t181 + t146 * t195 + t127;
t27 = t140 * t153 + t141 * t55;
t85 = t146 * t119 + t130;
t52 = t140 * t183 + t141 * t85;
t166 = -t27 * qJD(3) - t52 * qJDD(3) - t222;
t53 = t86 * pkin(4) - t87 * qJ(5) - t120;
t163 = -t134 * t91 + t237 - t77;
t26 = t140 * t55 - t141 * t153;
t51 = t140 * t85 - t141 * t183;
t161 = t26 * t74 - t27 * t72 - t52 * t43 + t51 * t44 - t237;
t149 = qJD(3) ^ 2;
t160 = pkin(7) * t149 - t180 - t230;
t159 = -t30 * t74 - qJDD(5) + t169;
t158 = -t26 * qJD(3) - t51 * qJDD(3) + t177;
t121 = -pkin(2) - t225;
t157 = pkin(1) * t191 + t119 * t149 + t121 * t133;
t156 = -pkin(7) * qJDD(3) + (t193 - t229) * qJD(3);
t155 = -qJDD(3) * t119 + (t121 * t134 - t195) * qJD(3);
t154 = (-g(1) * (-t120 + t172) - g(2) * t244) * t128;
t152 = -t223 * t72 + t224 * t74 - t61 * t43 + t60 * t44 - t237;
t150 = 0.2e1 * t74 * qJD(3) + t173;
t148 = cos(qJ(1));
t132 = t134 ^ 2;
t131 = t148 * pkin(1);
t124 = pkin(1) * t202;
t114 = -t141 * pkin(3) - pkin(4);
t112 = t140 * pkin(3) + qJ(5);
t104 = -t120 - t225;
t102 = qJDD(3) * t146 - t149 * t143;
t101 = qJDD(3) * t143 + t149 * t146;
t89 = t124 + t123;
t78 = t137 * t133 + 0.2e1 * t146 * t190;
t57 = -0.2e1 * t204 * t134 * qJD(3) + 0.2e1 * t143 * t210;
t47 = t53 - t225;
t40 = pkin(3) * t214 + t74 * pkin(4) + t72 * qJ(5);
t31 = t124 + t37;
t28 = (-t72 + t192) * qJD(3) + t162;
t1 = [qJDD(1), g(1) * t145 - g(2) * t148, g(1) * t148 + g(2) * t145, t133, (t133 * t147 - t191) * pkin(1) + t179, ((-qJDD(1) - t133) * t144 + (-qJD(1) - t134) * t201) * pkin(1) + t237, t78, t57, t101, t102, 0, t155 * t143 + (-t157 - t243) * t146 + t220, t155 * t146 + (t157 - t117) * t143 + t196, t104 * t43 + t89 * t72 + t158 + t240, t104 * t44 + t89 * t74 + t166 + t239, t161 + t234, t11 * t52 + t36 * t27 - t10 * t51 - t35 * t26 + t48 * t104 + t71 * t89 - g(1) * (t170 - t227) - g(2) * (t131 + t182), t31 * t72 + t47 * t43 + t158 + t242, t161 + t235, -t31 * t74 - t47 * t44 - t166 + t241, t8 * t52 + t34 * t27 + t7 * t47 + t30 * t31 + t9 * t51 + t33 * t26 - g(1) * (t215 - t227) - g(2) * (t131 + t197) + t154; 0, 0, 0, t133, t179 + t180, (-t198 + (-qJD(2) + t134) * t203) * pkin(1) + t237, t78, t57, t101, t102, 0, t156 * t143 + (-t160 - t243) * t146 + t220, t156 * t146 + (t160 - t117) * t143 + t196, -t72 * t194 - t120 * t43 + (t72 * t232 - t224) * qJD(3) + t167 + t240, -t74 * t194 - t120 * t44 + (t74 * t232 - t223) * qJD(3) - t175 + t239, t152 + t234, t11 * t61 - t10 * t60 - t48 * t120 - g(1) * t170 - g(2) * t182 + (-t194 + t123) * t71 + t223 * t36 - t224 * t35, -qJD(3) * t224 - t176 * t72 + t53 * t43 + t167 + t242, t152 + t235, qJD(3) * t223 + t176 * t74 - t53 * t44 + t175 + t241, -g(1) * t215 - g(2) * t197 - t176 * t30 + t223 * t34 + t224 * t33 + t7 * t53 + t9 * t60 + t8 * t61 + t154; 0, 0, 0, 0, 0, 0, -t143 * t132 * t146, t204 * t132, t211, t210, qJDD(3), t143 * t163 - t231, g(3) * t143 + t146 * t163, t209 - t71 * t74 + (qJDD(3) * t141 - t214 * t72) * pkin(3) + t169, t39 * qJD(3) + t71 * t72 + (-qJDD(3) * t140 - t214 * t74) * pkin(3) + t233, (t36 - t38) * t74 + (-t35 + t39) * t72 + (-t140 * t43 - t141 * t44) * pkin(3), t35 * t38 - t36 * t39 + (-t231 + t10 * t141 + t11 * t140 + (-t134 * t71 + t237) * t143) * pkin(3), t209 - t40 * t72 + (pkin(4) - t114) * qJDD(3) + t159, -t112 * t43 + t114 * t44 + (t34 - t38) * t74 + (t33 - t207) * t72, t112 * qJDD(3) - t30 * t72 + t40 * t74 + t136 + (0.2e1 * qJD(5) - t39) * qJD(3) - t233, t8 * t112 + t9 * t114 - t30 * t40 - t33 * t38 - g(3) * (-t172 + t226) + t207 * t34 + t237 * (pkin(4) * t125 - qJ(5) * t126 + t232); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, t28, t238, t35 * t74 + t36 * t72 - t236 + t48, t150, t238, -t28, t34 * t72 + (-qJD(5) - t33) * t74 + t151 - t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74 * t72 - qJDD(3), (t72 + t192) * qJD(3) + t162, -t69 - t149, -t34 * qJD(3) - t159 - t219;];
tau_reg = t1;
