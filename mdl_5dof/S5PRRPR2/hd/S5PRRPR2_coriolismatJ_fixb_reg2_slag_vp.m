% Calculate inertial parameters regressor of coriolis matrix for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRPR2_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:41
% EndTime: 2019-12-05 16:17:44
% DurationCPUTime: 1.46s
% Computational Cost: add. (1572->152), mult. (3529->229), div. (0->0), fcn. (2809->6), ass. (0->148)
t175 = qJD(2) + qJD(3);
t134 = sin(pkin(9));
t135 = cos(pkin(9));
t138 = cos(qJ(5));
t218 = t138 ^ 2;
t136 = sin(qJ(5));
t219 = t136 ^ 2;
t73 = (t218 / 0.2e1 + t219 / 0.2e1 - 0.1e1 / 0.2e1) * t135 * t134;
t225 = t175 * t73;
t132 = t134 ^ 2;
t101 = (-t218 + t219) * t132;
t224 = t175 * t101;
t133 = t135 ^ 2;
t124 = t132 + t133;
t102 = t124 * t136;
t58 = t175 * t102;
t103 = t124 * t138;
t59 = t175 * t103;
t223 = t175 * t124;
t156 = t175 * t135;
t137 = sin(qJ(3));
t214 = t137 * pkin(2);
t130 = qJ(4) + t214;
t222 = qJ(4) + t130;
t180 = t73 * qJD(1);
t192 = t132 * qJ(4);
t146 = -t135 * pkin(4) - t134 * pkin(7) - pkin(3);
t187 = t135 * t136;
t80 = qJ(4) * t187 - t138 * t146;
t186 = t135 * t138;
t81 = qJ(4) * t186 + t136 * t146;
t20 = t192 + (t136 * t80 + t81 * t138) * t135;
t139 = cos(qJ(3));
t216 = pkin(2) * t139;
t141 = t146 - t216;
t61 = t130 * t186 + t136 * t141;
t169 = t81 / 0.2e1 + t61 / 0.2e1;
t60 = t130 * t187 - t138 * t141;
t170 = t80 / 0.2e1 + t60 / 0.2e1;
t191 = t132 * t130;
t179 = t191 / 0.2e1 + t192 / 0.2e1;
t140 = (t170 * t136 + t169 * t138) * t135 + t179;
t182 = t138 * t139;
t185 = t136 * t137;
t92 = (t135 * t182 + t185) * pkin(2);
t203 = t92 * t136;
t183 = t138 * t137;
t184 = t136 * t139;
t91 = (-t135 * t184 + t183) * pkin(2);
t204 = t91 * t138;
t144 = -t204 / 0.2e1 - t203 / 0.2e1;
t4 = t140 + t144;
t221 = qJD(2) * t4 + qJD(3) * t20 + t180;
t15 = t191 + (t60 * t136 + t61 * t138) * t135;
t220 = qJD(2) * t15 + t180;
t217 = t130 / 0.2e1;
t215 = pkin(3) * t137;
t174 = t132 * t216;
t55 = t92 * t135 + t138 * t174;
t96 = t103 * qJD(4);
t212 = t55 * qJD(3) + t96;
t54 = t135 * t91 - t136 * t174;
t95 = t102 * qJD(4);
t211 = -t54 * qJD(3) + t95;
t210 = pkin(2) * qJD(2);
t209 = pkin(2) * qJD(3);
t208 = t60 * t135;
t207 = t61 * t135;
t206 = t80 * t135;
t205 = t81 * t135;
t116 = t124 * qJD(4);
t97 = t124 * t216;
t202 = t97 * qJD(3) + t116;
t22 = (t138 * t92 / 0.2e1 - t136 * t91 / 0.2e1 - t135 * t216 / 0.2e1) * t134;
t200 = qJD(2) * t22;
t190 = t132 * t136;
t32 = -t130 * t190 - t208;
t199 = qJD(2) * t32;
t189 = t132 * t138;
t104 = t130 * t189;
t33 = -t104 - t207;
t198 = qJD(2) * t33;
t34 = (t203 + t204) * t134;
t197 = qJD(2) * t34;
t196 = qJD(2) * t54;
t195 = qJD(2) * t55;
t89 = t124 * t130;
t194 = qJD(2) * t89;
t193 = qJD(2) * t97;
t188 = t133 * qJ(4);
t45 = (-t215 + (t89 - t214) * t139) * pkin(2);
t181 = t45 * qJD(2);
t178 = qJD(4) * t135;
t177 = qJD(5) * t136;
t176 = qJD(5) * t138;
t173 = t137 * t209;
t172 = t137 * t210;
t171 = t214 / 0.2e1;
t168 = t136 * t189;
t167 = t136 * t178;
t166 = t138 * t178;
t165 = t135 * t177;
t164 = t135 * t176;
t162 = -t184 / 0.2e1;
t161 = t183 / 0.2e1;
t160 = -t182 / 0.2e1;
t126 = qJ(4) * t189;
t159 = -t104 / 0.2e1 - t126 / 0.2e1;
t157 = pkin(2) * t175;
t155 = t134 * t172;
t154 = t137 * t157;
t153 = t136 * t156;
t152 = t138 * t156;
t9 = t130 * t174 - t60 * t91 + t61 * t92;
t151 = t22 * qJD(1) + t9 * qJD(2);
t10 = pkin(2) * t161 + (pkin(2) * t162 - t169) * t135 + t159;
t51 = -t126 - t205;
t149 = -qJD(2) * t10 - qJD(3) * t51;
t11 = (-t214 / 0.2e1 + (qJ(4) / 0.2e1 + t217) * t132) * t136 + (pkin(2) * t160 + t170) * t135;
t50 = -qJ(4) * t190 - t206;
t148 = -qJD(2) * t11 + qJD(3) * t50;
t114 = t124 * qJ(4);
t43 = t171 + t222 * (-t133 / 0.2e1 - t132 / 0.2e1);
t147 = qJD(2) * t43 - qJD(3) * t114;
t145 = -qJD(5) + t156;
t143 = t145 * t136;
t142 = t145 * t138;
t121 = qJ(4) * t174;
t120 = t134 * t173;
t115 = qJD(5) * t168;
t113 = t134 * t164;
t112 = t134 * t165;
t94 = t101 * qJD(5);
t86 = t175 * t168;
t83 = t134 * t152;
t82 = t134 * t153;
t72 = t134 * t142;
t71 = t134 * t143;
t62 = t73 * qJD(4);
t44 = t188 / 0.2e1 + t133 * t217 + t171 + t179;
t36 = t164 - t59;
t35 = t165 - t58;
t31 = t34 * qJD(3);
t13 = t205 / 0.2e1 + t207 / 0.2e1 + (t135 * t162 + t161) * pkin(2) - t159;
t12 = -t206 / 0.2e1 - t208 / 0.2e1 + (t135 * t160 - t185 / 0.2e1) * pkin(2) - t222 * t190 / 0.2e1;
t5 = qJD(3) * t22 + t62;
t3 = t140 - t144;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t200 + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t225; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134 * t176, t134 * t177, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t173, -t139 * t209, 0, 0, 0, 0, 0, 0, 0, 0, -t135 * t173, t120, t202, qJD(3) * t45 + qJD(4) * t89, -t115, t94, t112, t115, t113, 0, -qJD(5) * t33 + t211, qJD(5) * t32 + t212, -t31, qJD(3) * t9 + qJD(4) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154, -t139 * t157, 0, 0, 0, 0, 0, 0, 0, 0, -t135 * t154, t120 + t155, t193 + t202, t181 + (t121 + (t139 * t188 - t215) * pkin(2)) * qJD(3) + t44 * qJD(4), -t115, t94, t112, t115, t113, 0, qJD(5) * t13 - t196 + t211, qJD(5) * t12 + t195 + t212, -t31 - t197, (-t80 * t91 + t81 * t92 + t121) * qJD(3) + t3 * qJD(4) + t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t223, qJD(3) * t44 + t194, 0, 0, 0, 0, 0, 0, t58, t59, 0, qJD(3) * t3 + t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, t224, t71, t86, t72, 0, qJD(3) * t13 - qJD(5) * t61 - t198, qJD(3) * t12 + qJD(5) * t60 + t199, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t200 + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, t139 * t210, 0, 0, 0, 0, 0, 0, 0, 0, t135 * t172, -t155, t116 - t193, -qJD(4) * t43 - t181, -t115, t94, t112, t115, t113, 0, -qJD(5) * t10 + t196 + t95, -qJD(5) * t11 - t195 + t96, t197, qJD(4) * t4 - t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, t114 * qJD(4), -t115, t94, t112, t115, t113, 0, -qJD(5) * t51 + t95, qJD(5) * t50 + t96, 0, t20 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t223, -t147, 0, 0, 0, 0, 0, 0, t58, t59, 0, t221; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, t224, t71, t86, t72, 0, -qJD(5) * t81 + t149, qJD(5) * t80 + t148, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t225; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t223, qJD(3) * t43 - t194, 0, 0, 0, 0, 0, 0, t35, t36, 0, -qJD(3) * t4 - t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t223, t147, 0, 0, 0, 0, 0, 0, t35, t36, 0, -t221; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, t142, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, -t224, -t82, -t86, -t83, 0, qJD(3) * t10 - t167 + t198, qJD(3) * t11 - t166 - t199, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, -t224, -t82, -t86, -t83, 0, -t149 - t167, -t148 - t166, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t153, -t152, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
