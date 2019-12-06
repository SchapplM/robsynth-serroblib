% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPRR6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:57:57
% EndTime: 2019-12-05 15:58:01
% DurationCPUTime: 1.62s
% Computational Cost: add. (1474->184), mult. (3798->309), div. (0->0), fcn. (4332->10), ass. (0->166)
t133 = cos(pkin(10));
t229 = cos(qJ(4));
t180 = t229 * t133;
t131 = sin(pkin(10));
t135 = sin(qJ(4));
t204 = t135 * t131;
t150 = t180 - t204;
t106 = t150 ^ 2;
t181 = t229 * t131;
t203 = t135 * t133;
t111 = t181 + t203;
t107 = t111 ^ 2;
t232 = -t107 - t106;
t186 = t107 - t106;
t134 = sin(qJ(5));
t129 = t134 ^ 2;
t137 = cos(qJ(5));
t130 = t137 ^ 2;
t124 = t130 - t129;
t205 = t134 * t137;
t168 = 0.2e1 * t111 * t205;
t142 = qJD(2) * t168 - qJD(4) * t124;
t127 = t131 ^ 2;
t227 = t111 * pkin(4);
t228 = t150 * pkin(8);
t75 = t227 - t228;
t231 = t75 / 0.2e1;
t230 = -t111 / 0.2e1;
t226 = pkin(7) + qJ(3);
t215 = cos(pkin(5));
t170 = t215 * t131;
t132 = sin(pkin(5));
t136 = sin(qJ(2));
t207 = t132 * t136;
t103 = t133 * t207 + t170;
t141 = -t131 * t207 + t215 * t133;
t54 = t229 * t103 + t135 * t141;
t225 = t134 * t54;
t138 = cos(qJ(2));
t206 = t132 * t138;
t83 = t150 * t206;
t224 = t134 * t83;
t223 = t137 * t54;
t119 = t226 * t133;
t171 = t226 * t131;
t81 = t229 * t119 - t135 * t171;
t222 = t137 * t81;
t221 = t137 * t83;
t36 = t137 * t206 + t225;
t220 = t36 * t150;
t37 = -t134 * t206 + t223;
t219 = t37 * t150;
t218 = t81 * t134;
t82 = t111 * t206;
t217 = t82 * t134;
t216 = t82 * t137;
t105 = t181 / 0.2e1 + t203 / 0.2e1;
t43 = (t111 / 0.2e1 - t105) * t206;
t214 = qJD(1) * t43;
t143 = -t180 / 0.2e1 + t204 / 0.2e1;
t44 = (t150 / 0.2e1 + t143) * t206;
t213 = qJD(1) * t44;
t40 = t186 * t134;
t212 = qJD(2) * t40;
t41 = t232 * t134;
t211 = qJD(2) * t41;
t42 = t186 * t137;
t210 = qJD(2) * t42;
t70 = t232 * t137;
t209 = qJD(2) * t70;
t208 = t111 * t137;
t61 = t134 * t150;
t63 = t134 * t111;
t66 = t137 * t150;
t27 = ((-0.1e1 + t127) * t207 + (-t170 + t103) * t133) * t206;
t202 = t27 * qJD(1);
t201 = t186 * qJD(2);
t200 = t61 * qJD(2);
t199 = t63 * qJD(2);
t198 = t66 * qJD(2);
t120 = t133 ^ 2 + t127;
t197 = qJD(2) * t132;
t196 = qJD(3) * t137;
t195 = qJD(4) * t134;
t194 = qJD(4) * t137;
t193 = qJD(5) * t134;
t192 = qJD(5) * t137;
t191 = t105 * qJD(2);
t190 = t150 * qJD(2);
t104 = t150 * qJD(4);
t189 = t111 * qJD(2);
t188 = t111 * qJD(4);
t187 = t120 * qJD(2);
t185 = t134 * t207;
t184 = -t224 / 0.2e1;
t183 = -t221 / 0.2e1;
t53 = t103 * t135 - t229 * t141;
t182 = t53 * t230;
t126 = -pkin(3) * t133 - pkin(2);
t179 = t134 * t194;
t178 = t150 * t192;
t177 = t150 * t189;
t176 = t150 * t188;
t175 = t136 * t197;
t174 = t134 * t192;
t173 = t137 * t189;
t172 = t207 / 0.2e1;
t169 = t120 * t138;
t167 = qJD(2) * t126 + qJD(3);
t166 = -qJD(5) + t190;
t165 = -pkin(4) * t150 - t111 * pkin(8);
t164 = qJD(4) * t168;
t155 = (-t36 / 0.2e1 + t225 / 0.2e1) * t111;
t1 = t216 / 0.2e1 + t155;
t140 = t126 + t165;
t21 = -t137 * t140 + t218;
t5 = (-t21 + t218) * t111 - t75 * t66;
t162 = t1 * qJD(1) + t5 * qJD(2);
t154 = (-t37 / 0.2e1 + t223 / 0.2e1) * t111;
t4 = -t217 / 0.2e1 + t154;
t22 = t134 * t140 + t222;
t6 = (-t22 + t222) * t111 + t75 * t61;
t161 = t4 * qJD(1) + t6 * qJD(2);
t80 = t119 * t135 + t229 * t171;
t15 = -t150 * t21 - t63 * t80;
t151 = t172 + t182;
t8 = t183 + t220 / 0.2e1 - t151 * t134;
t160 = -qJD(1) * t8 + qJD(2) * t15;
t16 = t150 * t22 + t208 * t80;
t7 = t184 - t219 / 0.2e1 + t151 * t137;
t159 = qJD(1) * t7 - qJD(2) * t16;
t116 = t120 * qJ(3);
t139 = t103 * t133 / 0.2e1 - t141 * t131 / 0.2e1;
t38 = t172 - t139;
t158 = qJD(1) * t38 - qJD(2) * t116;
t157 = t166 * t137;
t156 = -t228 / 0.2e1 + t227 / 0.2e1;
t146 = t231 + t156;
t11 = t146 * t134;
t153 = pkin(4) * t194 - qJD(2) * t11;
t13 = t146 * t137;
t152 = pkin(4) * t195 + qJD(2) * t13;
t60 = (t129 / 0.2e1 - t130 / 0.2e1) * t111;
t149 = -qJD(2) * t60 + t179;
t148 = t111 * t157;
t147 = qJD(5) * t105 - t177;
t145 = qJD(2) * t107 * t205 + qJD(4) * t60;
t69 = t124 * t107;
t144 = qJD(2) * t69 + t164;
t102 = t105 * qJD(4);
t101 = t137 * t188;
t57 = t61 * qJD(5);
t56 = t60 * qJD(5);
t51 = -t193 + t200;
t46 = (-t105 + t230) * t206;
t45 = (-t150 / 0.2e1 + t143) * t206;
t39 = t172 + t139;
t20 = t53 * t137;
t19 = t53 * t134;
t14 = t80 * t134 + (-t156 + t231) * t137;
t12 = t80 * t137 + (-t75 / 0.2e1 + t156) * t134;
t10 = t219 / 0.2e1 + t53 * t208 / 0.2e1 + t184 + t137 * t172;
t9 = -t220 / 0.2e1 + t134 * t182 + t183 - t185 / 0.2e1;
t3 = t217 / 0.2e1 + t154;
t2 = -t216 / 0.2e1 + t155;
t17 = [0, 0, 0, 0, 0, 0, 0, t27 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t175, -t138 * t197, -t133 * t175, t131 * t175, t169 * t197, t202 + t39 * qJD(3) + (-pkin(2) * t136 + qJ(3) * t169) * t197, 0, 0, 0, 0, 0, qJD(4) * t46 - t150 * t175, qJD(4) * t45 + t111 * t175, 0, 0, 0, 0, 0, (-(t137 * t207 - t224) * t150 + t82 * t63) * qJD(2) + t2 * qJD(4) + t10 * qJD(5), ((t185 + t221) * t150 + t82 * t208) * qJD(2) + t3 * qJD(4) + t9 * qJD(5); 0, 0, 0, 0, 0, 0, 0, t39 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t46 - qJD(4) * t54, qJD(2) * t45 + qJD(4) * t53, 0, 0, 0, 0, 0, qJD(2) * t2 + qJD(5) * t19 - t194 * t54, qJD(2) * t3 + qJD(5) * t20 + t195 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t10 + qJD(4) * t19 - qJD(5) * t37, qJD(2) * t9 + qJD(4) * t20 + qJD(5) * t36; 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t38 - t202, 0, 0, 0, 0, 0, -t43 * qJD(4), -t44 * qJD(4), 0, 0, 0, 0, 0, qJD(4) * t1 - qJD(5) * t7, qJD(4) * t4 - qJD(5) * t8; 0, 0, 0, 0, 0, 0, t120 * qJD(3), t116 * qJD(3), t176, -t186 * qJD(4), 0, 0, 0, t126 * t188, t126 * t104, -t107 * t174 + t130 * t176, -qJD(5) * t69 - t150 * t164, t111 * t150 * t193 + qJD(4) * t42, -qJD(4) * t40 + t111 * t178, -t176, -qJD(3) * t41 + qJD(4) * t5 + qJD(5) * t16, -qJD(3) * t70 + qJD(4) * t6 + qJD(5) * t15; 0, 0, 0, 0, 0, 0, t187, -t158, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t211, -t209; 0, 0, 0, 0, 0, 0, 0, 0, t177, -t201, t104, -t188, 0, -qJD(4) * t81 + t126 * t189 - t214, qJD(4) * t80 + t126 * t190 - t213, -t56 - (-t130 * t189 - t179) * t150, -0.2e1 * t111 * t174 - t142 * t150, t134 * t188 + t210, t101 - t212, t147, (t134 * t165 - t222) * qJD(4) + t14 * qJD(5) + t162, (t137 * t165 + t218) * qJD(4) + t12 * qJD(5) + t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, -t144, t166 * t63, t148, t102, qJD(4) * t14 - qJD(5) * t22 - t159, qJD(4) * t12 + qJD(5) * t21 + t160; 0, 0, 0, 0, 0, 0, 0, t38 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t187, t158, 0, 0, 0, 0, 0, t188, t104, 0, 0, 0, 0, 0, t101 + t57 + t211, -qJD(4) * t63 + t178 + t209; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t190, 0, 0, 0, 0, 0, t173, -t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43 * qJD(2), t44 * qJD(2), 0, 0, 0, 0, 0, -qJD(2) * t1, -qJD(2) * t4; 0, 0, 0, 0, 0, 0, 0, 0, -t177, t201, 0, 0, 0, -t111 * t167 + t214, -t150 * t167 + t213, -t130 * t177 - t56, 0.2e1 * t134 * t148, -qJD(5) * t66 - t210, t57 + t212, -t147, -qJD(5) * t13 - t111 * t196 - t162, qJD(3) * t63 + qJD(5) * t11 - t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t189, -t190, 0, 0, 0, 0, 0, -t173, t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, t124 * qJD(5), 0, 0, 0, -pkin(4) * t193, -pkin(4) * t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, -t142, t192 - t198, t51, -t191, -pkin(8) * t192 - t152, pkin(8) * t193 - t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t7, qJD(2) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t145, t144, qJD(4) * t66 - t134 * t177, -qJD(4) * t61 - t150 * t173, t102, -qJD(3) * t61 + qJD(4) * t13 + t159, -qJD(4) * t11 - t150 * t196 - t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t200, -t137 * t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, t142, t198, -t200, t191, t152, t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t17;
