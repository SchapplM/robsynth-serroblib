% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPR1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:51
% EndTime: 2019-12-05 17:47:55
% DurationCPUTime: 2.05s
% Computational Cost: add. (3345->309), mult. (6566->375), div. (0->0), fcn. (4485->12), ass. (0->170)
t157 = sin(qJ(5));
t160 = cos(qJ(5));
t205 = qJD(5) * t160;
t206 = qJD(5) * t157;
t154 = sin(pkin(8));
t158 = sin(qJ(3));
t203 = qJD(1) * qJD(3);
t193 = t158 * t203;
t155 = cos(pkin(8));
t161 = cos(qJ(3));
t199 = t161 * qJDD(1);
t192 = t161 * t203;
t200 = t158 * qJDD(1);
t247 = t192 + t200;
t196 = t154 * t199 + t247 * t155;
t51 = t154 * t193 - t196;
t179 = -t154 * t200 + t155 * t199;
t102 = t154 * t161 + t155 * t158;
t97 = t102 * qJD(3);
t52 = qJD(1) * t97 - t179;
t92 = t102 * qJD(1);
t210 = qJD(1) * t158;
t195 = t154 * t210;
t209 = qJD(1) * t161;
t95 = t155 * t209 - t195;
t12 = -t157 * t51 + t160 * t52 + t92 * t205 + t95 * t206;
t146 = qJD(3) + qJD(5);
t44 = t157 * t95 + t160 * t92;
t220 = t44 * t146;
t250 = -t12 + t220;
t178 = -t157 * t92 + t160 * t95;
t228 = t178 ^ 2;
t229 = t44 ^ 2;
t249 = t228 - t229;
t227 = t44 * t178;
t159 = sin(qJ(1));
t162 = cos(qJ(1));
t242 = g(1) * t159 - g(2) * t162;
t13 = qJD(5) * t178 - t157 * t52 - t160 * t51;
t221 = t178 * t146;
t248 = -t13 + t221;
t165 = qJD(1) ^ 2;
t246 = -t165 * qJ(2) - t242;
t149 = qJDD(1) * qJ(2);
t232 = t95 * pkin(7);
t163 = -pkin(1) - pkin(6);
t115 = t163 * qJD(1) + qJD(2);
t85 = -qJ(4) * t210 + t158 * t115;
t68 = t154 * t85;
t86 = -qJ(4) * t209 + t161 * t115;
t72 = qJD(3) * pkin(3) + t86;
t34 = t155 * t72 - t68;
t23 = qJD(3) * pkin(4) - t232 + t34;
t233 = t92 * pkin(7);
t222 = t155 * t85;
t35 = t154 * t72 + t222;
t24 = t35 - t233;
t114 = t163 * qJDD(1) + qJDD(2);
t104 = t161 * t114;
t202 = qJD(1) * qJD(4);
t208 = qJD(3) * t158;
t33 = -t161 * t202 - t115 * t208 + qJDD(3) * pkin(3) + t104 + (t193 - t199) * qJ(4);
t207 = qJD(3) * t161;
t38 = (-qJ(4) * qJD(1) + t115) * t207 + (-qJ(4) * qJDD(1) + t114 - t202) * t158;
t14 = -t154 * t38 + t155 * t33;
t6 = qJDD(3) * pkin(4) + t52 * pkin(7) + t14;
t15 = t154 * t33 + t155 * t38;
t9 = t51 * pkin(7) + t15;
t1 = (qJD(5) * t23 + t9) * t160 + t157 * t6 - t24 * t206;
t147 = qJ(3) + pkin(8);
t137 = qJ(5) + t147;
t125 = sin(t137);
t126 = cos(t137);
t108 = pkin(3) * t210 + qJD(1) * qJ(2) + qJD(4);
t59 = t92 * pkin(4) + t108;
t245 = g(3) * t126 + t242 * t125 + t59 * t44 - t1;
t103 = -t154 * t158 + t155 * t161;
t176 = -t157 * t102 + t160 * t103;
t244 = t12 * t176;
t144 = qJDD(3) + qJDD(5);
t243 = t144 * t176;
t152 = t158 ^ 2;
t153 = t161 ^ 2;
t212 = t152 + t153;
t188 = t212 * t114;
t184 = g(1) * t162 + g(2) * t159;
t150 = qJD(1) * qJD(2);
t198 = 0.2e1 * t150;
t241 = 0.2e1 * t149 + t198 - t184;
t8 = t157 * t23 + t160 * t24;
t2 = -qJD(5) * t8 - t157 * t9 + t160 * t6;
t240 = g(3) * t125 - t242 * t126 - t178 * t59 + t2;
t94 = t154 * t208 - t155 * t207;
t169 = qJD(5) * t176 - t157 * t97 - t160 * t94;
t177 = t160 * t102 + t157 * t103;
t239 = t13 * t177 + t169 * t44;
t238 = -t144 * t177 - t146 * t169;
t237 = t1 * t177 + t169 * t8 + t176 * t2 - t242;
t236 = t95 ^ 2;
t234 = t51 * pkin(4);
t231 = pkin(3) * t154;
t230 = g(3) * t158;
t140 = t158 * pkin(3);
t226 = t95 * t92;
t156 = -qJ(4) - pkin(6);
t41 = -t154 * t86 - t222;
t27 = t41 + t233;
t42 = t155 * t86 - t68;
t28 = t42 - t232;
t127 = t155 * pkin(3) + pkin(4);
t89 = t157 * t127 + t160 * t231;
t225 = t89 * qJD(5) - t157 * t28 + t160 * t27;
t88 = t160 * t127 - t157 * t231;
t224 = -t88 * qJD(5) + t157 * t27 + t160 * t28;
t217 = qJ(4) - t163;
t82 = -t161 * qJD(4) + t217 * t208;
t110 = t217 * t161;
t83 = -qJD(3) * t110 - t158 * qJD(4);
t37 = t154 * t82 + t155 * t83;
t223 = -t97 * qJD(3) + t103 * qJDD(3);
t109 = t217 * t158;
t58 = -t155 * t109 - t154 * t110;
t219 = pkin(1) * qJDD(1);
t128 = qJ(2) + t140;
t216 = (t198 + t149) * qJ(2);
t215 = t162 * pkin(1) + t159 * qJ(2);
t213 = t152 - t153;
t164 = qJD(3) ^ 2;
t211 = -t164 - t165;
t204 = t108 * qJD(1);
t119 = pkin(3) * t207 + qJD(2);
t201 = qJDD(3) * t158;
t197 = t161 * t165 * t158;
t36 = -t154 * t83 + t155 * t82;
t189 = t157 * t94 - t160 * t97;
t57 = t154 * t109 - t155 * t110;
t187 = t212 * qJDD(1);
t186 = qJDD(2) - t219;
t185 = t158 * t192;
t181 = -t102 * t51 - t92 * t94;
t180 = -t103 * t52 - t95 * t97;
t39 = -t103 * pkin(7) + t57;
t40 = -t102 * pkin(7) + t58;
t16 = -t157 * t40 + t160 * t39;
t17 = t157 * t39 + t160 * t40;
t67 = t247 * pkin(3) + qJDD(4) + t149 + t150;
t175 = t94 * qJD(3) - t102 * qJDD(3);
t174 = 0.2e1 * qJ(2) * t203 + qJDD(3) * t163;
t172 = -t184 + t67;
t170 = t15 * t102 + t14 * t103 - t34 * t97 - t35 * t94 - t242;
t168 = -t163 * t164 + t241;
t145 = -pkin(7) + t156;
t139 = t162 * qJ(2);
t136 = qJDD(3) * t161;
t135 = cos(t147);
t134 = sin(t147);
t105 = pkin(4) * t134 + t140;
t90 = t92 ^ 2;
t73 = t102 * pkin(4) + t128;
t63 = pkin(3) * t209 + t95 * pkin(4);
t60 = -t94 * pkin(4) + t119;
t29 = t67 - t234;
t26 = t94 * pkin(7) + t37;
t25 = t97 * pkin(7) + t36;
t20 = -qJD(5) * t177 + t189;
t19 = t102 * t205 + t103 * t206 - t189;
t7 = -t157 * t24 + t160 * t23;
t4 = -qJD(5) * t17 - t157 * t26 + t160 * t25;
t3 = qJD(5) * t16 + t157 * t25 + t160 * t26;
t5 = [0, 0, 0, 0, 0, qJDD(1), t242, t184, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - t242 - 0.2e1 * t219, t241, -t186 * pkin(1) - g(1) * (-t159 * pkin(1) + t139) - g(2) * t215 + t216, t153 * qJDD(1) - 0.2e1 * t185, -0.2e1 * t158 * t199 + 0.2e1 * t213 * t203, -t164 * t158 + t136, t152 * qJDD(1) + 0.2e1 * t185, -t164 * t161 - t201, 0, t158 * t168 + t161 * t174, -t158 * t174 + t161 * t168, -t163 * t187 - t188 + t242, -g(1) * (t163 * t159 + t139) - g(2) * (t162 * pkin(6) + t215) + t163 * t188 + t216, t180, t52 * t102 + t103 * t51 + t97 * t92 + t95 * t94, t223, t181, t175, 0, t36 * qJD(3) + t57 * qJDD(3) + t67 * t102 - t108 * t94 + t119 * t92 - t128 * t51 - t134 * t184, -t37 * qJD(3) - t58 * qJDD(3) + t67 * t103 - t108 * t97 + t119 * t95 - t128 * t52 - t135 * t184, -t36 * t95 - t37 * t92 + t58 * t51 + t57 * t52 - t170, t15 * t58 + t35 * t37 + t14 * t57 + t34 * t36 + t67 * t128 + t108 * t119 - g(1) * (t162 * t140 + t139 + (-pkin(1) + t156) * t159) - g(2) * (t140 * t159 - t162 * t156 + t215), -t178 * t19 - t244, t12 * t177 - t13 * t176 - t169 * t178 + t19 * t44, -t19 * t146 + t243, t239, t238, 0, -t125 * t184 + t73 * t13 + t16 * t144 + t4 * t146 + t169 * t59 + t177 * t29 + t60 * t44, -t73 * t12 - t126 * t184 - t17 * t144 - t3 * t146 + t176 * t29 + t178 * t60 - t59 * t19, t16 * t12 - t17 * t13 - t178 * t4 + t7 * t19 - t3 * t44 - t237, t1 * t17 + t8 * t3 + t2 * t16 + t7 * t4 + t29 * t73 + t59 * t60 - g(1) * (t162 * t105 + t139 + (-pkin(1) + t145) * t159) - g(2) * (t159 * t105 - t162 * t145 + t215); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t165, t246 + t186, 0, 0, 0, 0, 0, 0, t211 * t158 + t136, t211 * t161 - t201, -t187, t188 + t246, 0, 0, 0, 0, 0, 0, -qJD(1) * t92 + t223, -qJD(1) * t95 + t175, -t180 - t181, t170 - t204, 0, 0, 0, 0, 0, 0, -qJD(1) * t44 + t20 * t146 + t243, -qJD(1) * t178 + t238, -t178 * t20 - t239 + t244, -t59 * qJD(1) + t7 * t20 + t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, -t213 * t165, t199, -t197, -t200, qJDD(3), t161 * t246 + t104 + t230, g(3) * t161 + (-t114 - t246) * t158, 0, 0, t226, -t90 + t236, t179, -t226, (t95 + t195) * qJD(3) - t196, qJDD(3), g(3) * t134 - t41 * qJD(3) - t108 * t95 - t242 * t135 + (qJDD(3) * t155 - t92 * t209) * pkin(3) + t14, g(3) * t135 + t42 * qJD(3) + t108 * t92 + t242 * t134 + (-qJDD(3) * t154 - t209 * t95) * pkin(3) - t15, (t35 + t41) * t95 + (-t34 + t42) * t92 + (t154 * t51 + t155 * t52) * pkin(3), -t34 * t41 - t35 * t42 + (t230 + t14 * t155 + t15 * t154 + (-t242 - t204) * t161) * pkin(3), t227, t249, t250, -t227, t248, t144, t88 * t144 - t146 * t225 - t63 * t44 + t240, -t89 * t144 + t146 * t224 - t178 * t63 + t245, t88 * t12 - t89 * t13 + (t225 + t8) * t178 + (t224 - t7) * t44, g(3) * t105 + t1 * t89 + t2 * t88 - t59 * t63 - t224 * t8 - t225 * t7 - t242 * (t161 * pkin(3) + pkin(4) * t135); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t95 - t195) * qJD(3) + t196, -0.2e1 * t92 * qJD(3) + t179, -t90 - t236, t34 * t95 + t35 * t92 + t172, 0, 0, 0, 0, 0, 0, t13 + t221, -t12 - t220, -t228 - t229, t178 * t7 + t44 * t8 + t172 - t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, t249, t250, -t227, t248, t144, t8 * t146 + t240, t7 * t146 + t245, 0, 0;];
tau_reg = t5;
