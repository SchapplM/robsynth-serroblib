% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPR16
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR16_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR16_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR16_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:39:37
% EndTime: 2019-12-31 18:39:42
% DurationCPUTime: 2.25s
% Computational Cost: add. (1978->336), mult. (3607->400), div. (0->0), fcn. (1904->6), ass. (0->183)
t107 = cos(qJ(3));
t110 = -pkin(1) - pkin(6);
t64 = t110 * qJD(1) + qJD(2);
t29 = (pkin(4) * qJD(1) - t64) * t107;
t180 = qJD(4) + t29;
t166 = qJD(1) * qJD(3);
t151 = t107 * t166;
t104 = sin(qJ(3));
t162 = t104 * qJDD(1);
t230 = t151 + t162;
t105 = sin(qJ(1));
t108 = cos(qJ(1));
t95 = g(2) * t108;
t226 = g(1) * t105 - t95;
t175 = qJD(1) * t107;
t177 = qJD(1) * t104;
t196 = pkin(3) * t177 + qJD(1) * qJ(2);
t32 = -qJ(4) * t175 + t196;
t229 = -qJD(1) * t32 - t226;
t112 = qJD(1) ^ 2;
t181 = t112 * qJ(2);
t228 = t226 + t181;
t98 = qJDD(1) * qJ(2);
t211 = pkin(4) - t110;
t103 = sin(qJ(5));
t106 = cos(qJ(5));
t109 = -pkin(3) - pkin(7);
t20 = t109 * qJD(3) + t180;
t215 = pkin(7) * t104;
t89 = t107 * qJ(4);
t128 = -t89 + t215;
t23 = t128 * qJD(1) + t196;
t6 = t103 * t20 + t106 * t23;
t152 = t104 * t166;
t99 = qJD(1) * qJD(2);
t138 = t230 * pkin(3) + qJ(4) * t152 + t98 + t99;
t145 = qJD(3) * pkin(7) - qJD(4);
t167 = qJ(4) * qJDD(1);
t7 = pkin(7) * t162 + (t145 * qJD(1) - t167) * t107 + t138;
t63 = t110 * qJDD(1) + qJDD(2);
t142 = -t107 * t63 + qJDD(4);
t172 = qJD(3) * t104;
t49 = t64 * t172;
t130 = t142 + t49;
t86 = t107 * qJDD(1);
t223 = -t152 + t86;
t9 = t223 * pkin(4) + t109 * qJDD(3) + t130;
t2 = -qJD(5) * t6 - t103 * t7 + t106 * t9;
t72 = qJD(5) + t175;
t227 = t6 * t72 + t2;
t101 = t104 ^ 2;
t102 = t107 ^ 2;
t178 = t101 + t102;
t148 = t178 * t63;
t176 = qJD(1) * t106;
t153 = t104 * t176;
t173 = qJD(3) * t103;
t46 = -t153 + t173;
t125 = t46 * t72;
t169 = qJD(5) * t103;
t14 = qJD(3) * t169 - qJD(5) * t153 - t106 * qJDD(3) - t230 * t103;
t225 = t14 - t125;
t171 = qJD(3) * t106;
t48 = t103 * t177 + t171;
t189 = qJD(5) * t48;
t15 = t103 * qJDD(3) - t230 * t106 + t189;
t212 = t48 * t72;
t224 = -t15 + t212;
t45 = -qJDD(5) - t223;
t30 = t106 * t45;
t121 = -t72 * t169 - t30;
t159 = 0.2e1 * t99;
t222 = t159 + 0.2e1 * t98;
t174 = qJD(3) * qJ(4);
t52 = t104 * t64;
t28 = -pkin(4) * t177 + t52;
t24 = t28 + t174;
t221 = -t109 * t45 + t24 * t72;
t165 = qJDD(3) * qJ(4);
t51 = t104 * t63;
t143 = -t51 - t165;
t199 = t107 * t64;
t18 = (-qJD(4) - t199) * qJD(3) + t143;
t188 = qJDD(3) * pkin(3);
t19 = t130 - t188;
t144 = -qJD(4) + t199;
t205 = qJD(3) * pkin(3);
t31 = -t144 - t205;
t33 = -t52 - t174;
t116 = -t18 * t104 - t19 * t107 + (t104 * t31 - t107 * t33) * qJD(3);
t163 = qJDD(3) * t110;
t208 = t104 * pkin(3) - t89;
t53 = qJ(2) + t208;
t220 = (qJD(1) * t53 + t32) * qJD(3) + t163;
t5 = -t103 * t23 + t106 * t20;
t1 = qJD(5) * t5 + t103 * t9 + t106 * t7;
t132 = t103 * t5 - t106 * t6;
t219 = -qJD(5) * t132 + t1 * t103 + t2 * t106;
t217 = t5 * t72;
t94 = g(3) * t104;
t214 = g(3) * t107;
t213 = t48 * t46;
t184 = t105 * t107;
t187 = t104 * t105;
t210 = pkin(3) * t184 + qJ(4) * t187;
t50 = pkin(3) * t175 + qJ(4) * t177;
t209 = (t159 + t98) * qJ(2);
t207 = t108 * pkin(1) + t105 * qJ(2);
t204 = t103 * t45;
t203 = t103 * t72;
t202 = t104 * t48;
t201 = t106 * t14;
t200 = t106 * t46;
t197 = t15 * t103;
t195 = pkin(1) * qJDD(1);
t194 = qJ(4) * t104;
t192 = qJD(3) * t24;
t191 = qJD(3) * t46;
t190 = qJD(3) * t48;
t186 = t104 * t108;
t185 = t104 * t112;
t183 = t106 * t107;
t182 = t107 * t108;
t179 = t101 - t102;
t170 = qJD(3) * t107;
t168 = qJD(5) * t106;
t164 = qJDD(3) * t104;
t161 = g(1) * t184 - g(2) * t182 - t94;
t160 = t108 * pkin(6) + t207;
t158 = pkin(3) * t170 + qJ(4) * t172 + qJD(2);
t157 = t110 * t105;
t156 = t72 * t173;
t155 = t72 * t171;
t149 = qJD(3) * t211;
t147 = pkin(3) * t187 + t160;
t111 = qJD(3) ^ 2;
t58 = qJDD(3) * t107 - t111 * t104;
t56 = t178 * qJDD(1);
t141 = qJDD(2) - t195;
t140 = -t208 - t215;
t139 = t104 * t151;
t136 = g(1) * t108 + g(2) * t105;
t133 = t103 * t6 + t106 * t5;
t129 = -pkin(3) * t107 - t194;
t40 = qJ(2) - t140;
t55 = t211 * t107;
t17 = t103 * t55 + t106 * t40;
t16 = -t103 * t40 + t106 * t55;
t90 = t108 * qJ(2);
t124 = pkin(3) * t186 - qJ(4) * t182 + t90;
t122 = t72 * t168 - t204;
t120 = 0.2e1 * qJ(2) * t166 + t163;
t119 = -t110 * t56 + t226;
t118 = -t110 * t111 - t136;
t117 = t32 * t175 + t142 + t161;
t115 = t118 + t222;
t13 = (-qJD(1) * qJD(4) - t167) * t107 + t138;
t26 = -qJD(4) * t107 + t158;
t114 = -qJD(1) * t26 - qJDD(1) * t53 - t118 - t13;
t10 = -pkin(4) * t162 + (qJD(4) - t29) * qJD(3) - t143;
t113 = -qJD(5) * t109 * t72 - t104 * t226 + t10 - t214;
t69 = t107 * t185;
t59 = t179 * t112;
t57 = t107 * t111 + t164;
t54 = t211 * t104;
t44 = t107 * t149;
t43 = t104 * t149;
t42 = qJDD(1) * t102 - 0.2e1 * t139;
t41 = qJDD(1) * t101 + 0.2e1 * t139;
t39 = t58 - t185;
t38 = t164 + (t111 + t112) * t107;
t37 = -t103 * t184 + t106 * t108;
t36 = -t103 * t108 - t105 * t183;
t35 = -t103 * t182 - t105 * t106;
t34 = t103 * t105 - t106 * t182;
t27 = pkin(7) * t175 + t50;
t22 = -0.2e1 * t104 * t86 + 0.2e1 * t179 * t166;
t21 = t145 * t107 + t158;
t12 = t103 * t28 + t106 * t27;
t11 = -t103 * t27 + t106 * t28;
t4 = -t17 * qJD(5) - t103 * t21 - t106 * t43;
t3 = t16 * qJD(5) - t103 * t43 + t106 * t21;
t8 = [0, 0, 0, 0, 0, qJDD(1), t226, t136, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - 0.2e1 * t195 - t226, -t136 + t222, -t141 * pkin(1) - g(1) * (-t105 * pkin(1) + t90) - g(2) * t207 + t209, t42, t22, t58, t41, -t57, 0, t115 * t104 + t107 * t120, -t120 * t104 + t107 * t115, t119 - t148, -g(1) * (t90 + t157) - g(2) * t160 + t110 * t148 + t209, 0, -t58, t57, t42, t22, t41, -t116 + t119, t114 * t104 - t220 * t107, t220 * t104 + t114 * t107, t13 * t53 + t32 * t26 - g(1) * (t157 + t124) - g(2) * (-t105 * t89 + t147) + t116 * t110, t168 * t202 + (-t104 * t14 + t170 * t48) * t103, (-t103 * t46 + t106 * t48) * t170 + (-t197 - t201 + (-t103 * t48 - t200) * qJD(5)) * t104, (-t14 + t156) * t107 + (t122 - t190) * t104, -t170 * t200 + (-t106 * t15 + t169 * t46) * t104, (-t15 + t155) * t107 + (t121 + t191) * t104, -t107 * t45 - t172 * t72, -g(1) * t35 - g(2) * t37 - t15 * t54 - t16 * t45 + t4 * t72 - t44 * t46 + (-t171 * t24 + t2) * t107 + (-qJD(3) * t5 - t10 * t106 + t169 * t24) * t104, -g(1) * t34 - g(2) * t36 + t14 * t54 + t17 * t45 - t3 * t72 - t44 * t48 + (t173 * t24 - t1) * t107 + (qJD(3) * t6 + t10 * t103 + t168 * t24) * t104, t14 * t16 - t15 * t17 - t3 * t46 - t4 * t48 - t132 * t170 + (-qJD(5) * t133 + t1 * t106 - t103 * t2 - t136) * t104, t1 * t17 + t6 * t3 + t2 * t16 + t5 * t4 - t10 * t54 - t24 * t44 - g(1) * (pkin(7) * t186 + t124) - g(2) * (pkin(4) * t108 + t147) + (g(1) * t211 - g(2) * t128) * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t112, -t228 + t141, 0, 0, 0, 0, 0, 0, t39, -t38, -t56, t148 - t228, 0, 0, 0, 0, 0, 0, -t56, -t39, t38, t116 + t229, 0, 0, 0, 0, 0, 0, qJD(1) * t203 + (t15 + t155) * t104 + (-t121 + t191) * t107, t72 * t176 + (-t14 - t156) * t104 + (t122 + t190) * t107, (-t48 * t172 + qJD(1) * t46 + (qJD(5) * t46 - t14) * t107) * t106 + (-t46 * t172 - qJD(1) * t48 + (t15 - t189) * t107) * t103, t132 * qJD(1) + (qJD(3) * t133 + t10) * t104 + (t192 - t219) * t107 - t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t59, t86, -t69, -t162, qJDD(3), (t63 - t181) * t107 - t161, t228 * t104 + t214 - t51, 0, 0, qJDD(3), -t86, t162, t69, -t59, -t69, t129 * qJDD(1) + ((-t33 - t174) * t107 + (-qJD(4) + t31 + t205) * t104) * qJD(1), t177 * t50 + t117 - 0.2e1 * t188, 0.2e1 * t165 + 0.2e1 * qJD(3) * qJD(4) + t51 + (qJD(1) * t50 - g(3)) * t107 + t229 * t104, -t19 * pkin(3) - g(1) * t210 + g(3) * t208 - t18 * qJ(4) - t129 * t95 + t144 * t33 - t31 * t52 - t32 * t50, -t203 * t48 - t201, (-t15 - t212) * t106 + (t14 + t125) * t103, (-t107 * t203 + t202) * qJD(1) + t121, t106 * t125 + t197, (-t104 * t46 - t183 * t72) * qJD(1) - t122, t72 * t177, qJ(4) * t15 + t113 * t103 + t221 * t106 - t11 * t72 + t5 * t177 + t180 * t46, -qJ(4) * t14 - t221 * t103 + t113 * t106 + t12 * t72 - t6 * t177 + t180 * t48, t11 * t48 + t12 * t46 + (-t6 * t175 + t109 * t14 - t2 + (-t109 * t46 - t6) * qJD(5)) * t106 + (t5 * t175 - t109 * t15 - t1 + (t109 * t48 + t5) * qJD(5)) * t103 - t161, t10 * qJ(4) - t6 * t12 - t5 * t11 - g(1) * (pkin(7) * t184 + t210) - g(3) * t140 + t180 * t24 - (t107 * t109 - t194) * t95 + t219 * t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, qJDD(3) - t69, -t102 * t112 - t111, qJD(3) * t33 + t117 - t188 + t49, 0, 0, 0, 0, 0, 0, -t203 * t72 - t191 - t30, -t106 * t72 ^ 2 - t190 + t204, t224 * t103 + t225 * t106, -t192 + t227 * t106 + (t1 - t217) * t103 + t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, -t46 ^ 2 + t48 ^ 2, -t225, -t213, t224, -t45, -g(1) * t36 + g(2) * t34 - t106 * t94 - t24 * t48 + t227, g(1) * t37 - g(2) * t35 + t24 * t46 + t217 + (-qJD(5) * t20 - t7) * t106 + (qJD(5) * t23 - t9 + t94) * t103, 0, 0;];
tau_reg = t8;
