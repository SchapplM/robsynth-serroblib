% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:48
% EndTime: 2019-12-05 16:19:52
% DurationCPUTime: 1.95s
% Computational Cost: add. (2924->296), mult. (6621->371), div. (0->0), fcn. (4825->12), ass. (0->164)
t152 = sin(qJ(5));
t220 = cos(qJ(5));
t149 = sin(pkin(9));
t150 = cos(pkin(9));
t153 = sin(qJ(3));
t154 = cos(qJ(3));
t102 = t149 * t154 + t150 * t153;
t226 = t102 * qJD(2);
t196 = qJD(2) * t154;
t185 = t150 * t196;
t197 = qJD(2) * t153;
t96 = t149 * t197 - t185;
t170 = t152 * t96 - t220 * t226;
t190 = t154 * qJDD(2);
t191 = t153 * qJDD(2);
t172 = t149 * t191 - t150 * t190;
t230 = t226 * qJD(3);
t53 = t172 + t230;
t192 = qJD(2) * qJD(3);
t183 = t153 * t192;
t164 = t102 * qJDD(2) - t149 * t183;
t182 = t154 * t192;
t54 = t150 * t182 + t164;
t160 = t170 * qJD(5) - t152 * t54 - t220 * t53;
t144 = qJD(3) + qJD(5);
t202 = t170 * t144;
t233 = t160 - t202;
t184 = qJD(5) * t220;
t195 = qJD(5) * t152;
t168 = -t152 * t53 - t96 * t184 - t195 * t226 + t220 * t54;
t48 = -t152 * t226 - t220 * t96;
t201 = t48 * t144;
t232 = t168 - t201;
t211 = t48 ^ 2;
t212 = t170 ^ 2;
t231 = -t211 + t212;
t210 = t48 * t170;
t143 = pkin(8) + qJ(2);
t133 = sin(t143);
t135 = cos(t143);
t175 = g(1) * t135 + g(2) * t133;
t222 = t226 * pkin(7);
t208 = qJ(4) + pkin(6);
t118 = t208 * t154;
t194 = t153 * qJD(1);
t92 = qJD(2) * t118 + t194;
t71 = t149 * t92;
t204 = qJD(3) * pkin(3);
t117 = t208 * t153;
t140 = t154 * qJD(1);
t90 = -qJD(2) * t117 + t140;
t80 = t90 + t204;
t36 = t150 * t80 - t71;
t24 = qJD(3) * pkin(4) - t222 + t36;
t223 = t96 * pkin(7);
t203 = t150 * t92;
t37 = t149 * t80 + t203;
t25 = t37 - t223;
t138 = t154 * qJDD(1);
t178 = qJD(3) * t208;
t176 = t154 * t178;
t193 = qJD(1) * qJD(3);
t35 = qJDD(3) * pkin(3) + t138 - qJD(2) * t176 + (-qJD(2) * qJD(4) - t208 * qJDD(2) - t193) * t153;
t186 = pkin(6) * t190 + t153 * qJDD(1) + t154 * t193;
t91 = t154 * qJD(4) - t153 * t178;
t38 = qJ(4) * t190 + t91 * qJD(2) + t186;
t17 = -t149 * t38 + t150 * t35;
t7 = qJDD(3) * pkin(4) - t54 * pkin(7) + t17;
t18 = t149 * t35 + t150 * t38;
t8 = -t53 * pkin(7) + t18;
t1 = t152 * t7 + t24 * t184 - t25 * t195 + t220 * t8;
t145 = qJ(3) + pkin(9);
t139 = qJ(5) + t145;
t127 = sin(t139);
t128 = cos(t139);
t213 = t154 * pkin(3);
t131 = pkin(2) + t213;
t114 = -t131 * qJD(2) + qJD(4);
t60 = t96 * pkin(4) + t114;
t229 = g(3) * t127 + t175 * t128 - t60 * t48 - t1;
t10 = t152 * t24 + t220 * t25;
t2 = -t10 * qJD(5) - t152 * t8 + t220 * t7;
t228 = -g(3) * t128 + t175 * t127 + t170 * t60 + t2;
t227 = g(1) * t133 - g(2) * t135;
t225 = t226 ^ 2;
t155 = qJD(3) ^ 2;
t224 = t53 * pkin(4);
t101 = t149 * t153 - t150 * t154;
t100 = t101 * qJD(3);
t167 = t102 * qJD(3);
t21 = t220 * t100 + t101 * t184 + t102 * t195 + t152 * t167;
t58 = -t152 * t101 + t220 * t102;
t221 = t160 * t58 - t21 * t48;
t219 = pkin(3) * t149;
t215 = g(3) * t154;
t214 = t153 * pkin(3);
t209 = t226 * t96;
t207 = t100 * t96 - t102 * t53;
t41 = t150 * t90 - t71;
t93 = -t153 * qJD(4) - t176;
t42 = t149 * t93 + t150 * t91;
t39 = -t149 * t90 - t203;
t26 = t39 + t223;
t28 = t41 - t222;
t129 = t150 * pkin(3) + pkin(4);
t88 = t220 * t129 - t152 * t219;
t206 = t88 * qJD(5) - t152 * t26 - t220 * t28;
t89 = t152 * t129 + t220 * t219;
t205 = -t89 * qJD(5) + t152 * t28 - t220 * t26;
t200 = pkin(6) * qJDD(2);
t199 = qJDD(2) * pkin(2);
t62 = -t149 * t117 + t150 * t118;
t146 = t153 ^ 2;
t147 = t154 ^ 2;
t198 = t146 - t147;
t189 = pkin(6) * t197;
t188 = pkin(6) * t196;
t156 = qJD(2) ^ 2;
t187 = t153 * t156 * t154;
t136 = cos(t145);
t181 = pkin(4) * t136 + t213;
t40 = -t149 * t91 + t150 * t93;
t61 = -t150 * t117 - t149 * t118;
t177 = t153 * t182;
t22 = t58 * qJD(5) - t152 * t100 + t220 * t167;
t57 = t220 * t101 + t152 * t102;
t173 = t168 * t57 - t170 * t22;
t141 = qJDD(3) + qJDD(5);
t171 = t58 * t141 - t21 * t144;
t43 = -t102 * pkin(7) + t61;
t44 = -t101 * pkin(7) + t62;
t19 = -t152 * t44 + t220 * t43;
t20 = t152 * t43 + t220 * t44;
t169 = -0.2e1 * pkin(2) * t192 - pkin(6) * qJDD(3);
t165 = pkin(2) * t156 + t175;
t86 = pkin(3) * t183 - t131 * qJDD(2) + qJDD(4);
t163 = -pkin(6) * t155 + 0.2e1 * t199 + t227;
t162 = -t54 * t101 - t167 * t226;
t161 = -t227 + t86;
t112 = t140 - t189;
t113 = t188 + t194;
t63 = -pkin(6) * t183 + t186;
t64 = -t153 * t193 + t138 + (-t182 - t191) * pkin(6);
t157 = -t64 * t153 + t63 * t154 + (-t112 * t154 - t113 * t153) * qJD(3) - t175;
t148 = qJDD(1) - g(3);
t142 = -pkin(7) - t208;
t134 = sin(t145);
t116 = qJDD(3) * t154 - t155 * t153;
t115 = qJDD(3) * t153 + t155 * t154;
t110 = pkin(2) + t181;
t94 = t96 ^ 2;
t70 = t101 * pkin(4) - t131;
t66 = (t102 * pkin(4) + t214) * qJD(3);
t65 = pkin(3) * t197 + pkin(4) * t226;
t56 = -t100 * qJD(3) + t102 * qJDD(3);
t55 = -t101 * qJDD(3) - t102 * t155;
t30 = t86 + t224;
t29 = -pkin(7) * t167 + t42;
t27 = t100 * pkin(7) + t40;
t13 = -t57 * t141 - t22 * t144;
t9 = -t152 * t25 + t220 * t24;
t4 = -t20 * qJD(5) - t152 * t29 + t220 * t27;
t3 = t19 * qJD(5) + t152 * t27 + t220 * t29;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t148, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, 0, 0, 0, 0, 0, 0, t116, -t115, 0, t63 * t153 + t64 * t154 - g(3) + (-t112 * t153 + t113 * t154) * qJD(3), 0, 0, 0, 0, 0, 0, t55, -t56, -t162 + t207, -t37 * t100 - t17 * t101 + t18 * t102 - t36 * t167 - g(3), 0, 0, 0, 0, 0, 0, t13, -t171, t173 + t221, t1 * t58 - t10 * t21 - t2 * t57 - t9 * t22 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t227, t175, 0, 0, t146 * qJDD(2) + 0.2e1 * t177, 0.2e1 * t153 * t190 - 0.2e1 * t198 * t192, t115, t147 * qJDD(2) - 0.2e1 * t177, t116, 0, t169 * t153 + t163 * t154, -t163 * t153 + t169 * t154, (t146 + t147) * t200 + t157, (t227 + t199) * pkin(2) + t157 * pkin(6), -t100 * t226 + t54 * t102, t162 + t207, t56, t53 * t101 + t96 * t167, t55, 0, t61 * qJDD(3) + t86 * t101 - t131 * t53 + t227 * t136 + (t114 * t102 + t96 * t214 + t40) * qJD(3), -t62 * qJDD(3) - t114 * t100 + t86 * t102 - t131 * t54 - t227 * t134 + (t214 * t226 - t42) * qJD(3), t36 * t100 - t18 * t101 - t17 * t102 - t37 * t167 - t226 * t40 - t42 * t96 - t62 * t53 - t61 * t54 - t175, t18 * t62 + t37 * t42 + t17 * t61 + t36 * t40 - t86 * t131 + t114 * t153 * t204 - g(1) * (-t133 * t131 + t135 * t208) - g(2) * (t135 * t131 + t133 * t208), t168 * t58 + t170 * t21, -t173 + t221, t171, -t160 * t57 - t22 * t48, t13, 0, t128 * t227 + t19 * t141 + t4 * t144 - t160 * t70 + t60 * t22 + t30 * t57 - t48 * t66, -t127 * t227 - t20 * t141 - t3 * t144 + t168 * t70 - t170 * t66 - t60 * t21 + t30 * t58, -t1 * t57 - t10 * t22 + t160 * t20 - t168 * t19 + t170 * t4 - t2 * t58 + t9 * t21 + t3 * t48 - t175, t1 * t20 + t10 * t3 + t2 * t19 + t9 * t4 + t30 * t70 + t60 * t66 - g(1) * (-t133 * t110 - t135 * t142) - g(2) * (t135 * t110 - t133 * t142); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t187, t198 * t156, t191, t187, t190, qJDD(3), -t215 + t138 + (t113 - t188) * qJD(3) + (t165 - t193 - t200) * t153, g(3) * t153 + (t112 + t189) * qJD(3) + t165 * t154 - t186, 0, 0, t209, -t94 + t225, (t96 + t185) * qJD(3) + t164, -t209, -t172, qJDD(3), -g(3) * t136 - t39 * qJD(3) - t114 * t226 + t175 * t134 + (qJDD(3) * t150 - t96 * t197) * pkin(3) + t17, g(3) * t134 + t41 * qJD(3) + t114 * t96 + t175 * t136 + (-qJDD(3) * t149 - t197 * t226) * pkin(3) - t18, (t37 + t39) * t226 + (-t36 + t41) * t96 + (-t149 * t53 - t150 * t54) * pkin(3), -t36 * t39 - t37 * t41 + (-t215 + t149 * t18 + t150 * t17 + (-qJD(2) * t114 + t175) * t153) * pkin(3), t210, t231, t232, -t210, t233, t141, t88 * t141 + t205 * t144 + t48 * t65 + t228, -t89 * t141 - t206 * t144 + t170 * t65 + t229, t89 * t160 - t168 * t88 + (t206 + t9) * t48 + (-t10 + t205) * t170, t1 * t89 + t2 * t88 - t60 * t65 - g(3) * t181 + t205 * t9 - t175 * (-pkin(4) * t134 - t214) + t206 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172 + 0.2e1 * t230, (-t96 + t185) * qJD(3) + t164, -t94 - t225, t226 * t36 + t37 * t96 + t161, 0, 0, 0, 0, 0, 0, -t160 - t202, t168 + t201, -t211 - t212, -t10 * t48 - t170 * t9 + t161 + t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, t231, t232, -t210, t233, t141, t10 * t144 + t228, t9 * t144 + t229, 0, 0;];
tau_reg = t5;
