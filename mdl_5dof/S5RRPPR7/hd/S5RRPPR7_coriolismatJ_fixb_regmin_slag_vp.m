% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x25]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:00
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPPR7_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:59:04
% EndTime: 2021-01-15 19:59:11
% DurationCPUTime: 1.94s
% Computational Cost: add. (1703->200), mult. (3467->289), div. (0->0), fcn. (3725->6), ass. (0->170)
t116 = cos(qJ(5));
t113 = sin(pkin(8));
t115 = sin(qJ(2));
t218 = -qJ(3) - pkin(6);
t98 = t218 * t115;
t212 = t113 * t98;
t202 = cos(pkin(8));
t117 = cos(qJ(2));
t99 = t218 * t117;
t97 = t202 * t99;
t231 = -t97 + t212;
t106 = t202 * t117;
t199 = t113 * t115;
t93 = -t106 + t199;
t235 = -t93 * pkin(4) + t231;
t207 = t116 * t235;
t114 = sin(qJ(5));
t164 = qJD(5) * t114;
t95 = t113 * t117 + t202 * t115;
t48 = t114 * t95;
t185 = t48 * qJD(1);
t236 = t185 + t164;
t225 = t95 ^ 2;
t90 = t93 ^ 2;
t232 = t90 + t225;
t234 = qJD(3) * t232;
t233 = t232 * qJD(1);
t60 = t113 * t99 + t202 * t98;
t134 = -t231 * t93 - t60 * t95;
t229 = qJD(3) * t134;
t227 = t134 * qJD(1);
t37 = -t95 * pkin(4) + t60;
t226 = t37 * t93;
t111 = t114 ^ 2;
t112 = t116 ^ 2;
t102 = t111 - t112;
t51 = t116 * t93;
t141 = 0.2e1 * t114 * t51;
t121 = qJD(1) * t141 - t102 * qJD(2);
t224 = -t95 / 0.2e1;
t138 = -t97 / 0.2e1;
t223 = pkin(3) + pkin(7);
t221 = t95 * pkin(3);
t220 = -t114 / 0.2e1;
t110 = t115 * pkin(2);
t216 = qJD(2) * pkin(2);
t108 = -t117 * pkin(2) - pkin(1);
t130 = -t95 * qJ(4) + t108;
t24 = t223 * t93 + t130;
t13 = t114 * t24 + t37 * t116;
t203 = t93 * qJ(4);
t140 = t110 + t203;
t25 = t223 * t95 + t140;
t1 = -t226 * t116 + t13 * t93 - t25 * t48;
t215 = t1 * qJD(1);
t105 = t113 * pkin(2) + qJ(4);
t214 = t105 * t93;
t213 = t105 * t95;
t211 = t116 * t25;
t209 = t37 * t114;
t14 = t116 * t24 - t209;
t2 = t226 * t114 + t14 * t93 - t95 * t211;
t210 = t2 * qJD(1);
t208 = t235 * t114;
t47 = t93 * pkin(3) + t130;
t53 = t140 + t221;
t7 = t47 * t53;
t206 = t7 * qJD(1);
t8 = t13 * t95 + t93 * t207;
t205 = t8 * qJD(1);
t9 = -t14 * t95 + t93 * t208;
t204 = t9 * qJD(1);
t12 = t108 * t110;
t198 = t12 * qJD(1);
t15 = -t47 * t95 - t53 * t93;
t197 = t15 * qJD(1);
t16 = t47 * t93 - t53 * t95;
t196 = t16 * qJD(1);
t107 = -t202 * pkin(2) - pkin(3);
t109 = t110 / 0.2e1;
t19 = t109 + (pkin(3) / 0.2e1 - t107 / 0.2e1) * t95 + (qJ(4) / 0.2e1 + t105 / 0.2e1) * t93;
t193 = t19 * qJD(1);
t21 = t232 * t114;
t192 = t21 * qJD(1);
t157 = t225 - t90;
t26 = t157 * t114;
t191 = t26 * qJD(1);
t27 = t157 * t116;
t190 = t27 * qJD(1);
t28 = t232 * t116;
t189 = t28 * qJD(1);
t118 = (-t113 * t93 / 0.2e1 + t202 * t224) * pkin(2);
t35 = -t110 / 0.2e1 + t118;
t188 = t35 * qJD(1);
t40 = t108 * t95 + t93 * t110;
t187 = t40 * qJD(1);
t41 = -t108 * t93 + t95 * t110;
t186 = t41 * qJD(1);
t184 = t48 * qJD(5);
t183 = t51 * qJD(1);
t182 = t51 * qJD(2);
t179 = t231 * qJD(2);
t178 = t60 * qJD(2);
t87 = t199 / 0.2e1 - t106 / 0.2e1;
t177 = t87 * qJD(1);
t176 = t225 * qJD(1);
t175 = t93 * qJD(1);
t82 = t93 * qJD(2);
t174 = t93 * qJD(3);
t173 = t93 * qJD(4);
t172 = t95 * qJD(1);
t171 = t95 * qJD(2);
t84 = t95 * qJD(3);
t170 = t95 * qJD(4);
t169 = qJD(1) * t114;
t168 = qJD(1) * t116;
t167 = qJD(1) * t117;
t166 = qJD(4) * t114;
t165 = qJD(4) * t116;
t163 = qJD(5) * t116;
t103 = -t115 ^ 2 + t117 ^ 2;
t162 = t103 * qJD(1);
t161 = t114 * qJD(2);
t160 = t115 * qJD(2);
t159 = t116 * qJD(2);
t158 = t117 * qJD(2);
t156 = pkin(1) * t115 * qJD(1);
t155 = pkin(1) * t167;
t154 = t47 * t172;
t153 = t93 * t172;
t152 = t93 * t171;
t150 = t111 * t175;
t149 = t95 * t164;
t148 = t95 * t163;
t147 = t93 * t169;
t146 = t93 * t161;
t145 = t225 * t169;
t77 = t95 * t168;
t144 = t115 * t167;
t143 = t114 * t163;
t142 = t114 * t159;
t139 = qJD(5) + t172;
t136 = qJD(2) * t141;
t104 = -pkin(7) + t107;
t133 = t104 * t93 + t213;
t132 = t139 * t114;
t57 = t138 + t97 / 0.2e1;
t131 = t57 * qJD(1) + t105 * qJD(2);
t129 = t104 * t224 + t214 / 0.2e1;
t128 = t87 * qJD(5) + t153;
t127 = t93 * t132;
t119 = t25 / 0.2e1 + t129;
t5 = t119 * t116;
t126 = -t5 * qJD(1) + t105 * t161;
t3 = t119 * t114;
t125 = -t3 * qJD(1) - t105 * t159;
t50 = (t112 / 0.2e1 - t111 / 0.2e1) * t93;
t124 = -t50 * qJD(1) + t142;
t123 = t114 * t90 * t168 + t50 * qJD(2);
t54 = t102 * t90;
t122 = t54 * qJD(1) + t136;
t78 = t87 * qJD(2);
t66 = -t77 - t163;
t45 = t50 * qJD(5);
t39 = 0.2e1 * t138 + t212;
t34 = t109 + t118;
t20 = -t214 / 0.2e1 + t107 * t95 / 0.2e1 + t109 + t203 / 0.2e1 + t221 / 0.2e1;
t6 = -t208 / 0.2e1 - t211 / 0.2e1 + t235 * t220 + t129 * t116;
t4 = t129 * t114 + t25 * t220 + t207;
t10 = [0, 0, 0, t115 * t158, t103 * qJD(2), 0, 0, 0, -pkin(1) * t160, -pkin(1) * t158, t40 * qJD(2), t41 * qJD(2), t234, t12 * qJD(2) + t229, t234, t15 * qJD(2) + t93 * t170, t16 * qJD(2) + qJD(4) * t225, t7 * qJD(2) - t47 * t170 + t229, t111 * t152 + t143 * t90, -t54 * qJD(5) + t136 * t95, t26 * qJD(2) + t148 * t93, t27 * qJD(2) - t149 * t93, -t152, t1 * qJD(2) + t28 * qJD(3) + t9 * qJD(5) + t166 * t225, t2 * qJD(2) - t21 * qJD(3) + t8 * qJD(5) + t165 * t225; 0, 0, 0, t144, t162, t158, -t160, 0, -pkin(6) * t158 - t156, pkin(6) * t160 - t155, -t179 + t187, -t178 + t186, (-t113 * t95 + t202 * t93) * t216, t198 + (t113 * t60 - t202 * t231) * t216 + t34 * qJD(3), (-t107 * t93 - t213) * qJD(2) - t173, t179 + t197, t178 + t196, t206 + (t60 * t105 + t107 * t231) * qJD(2) + t20 * qJD(3) + t39 * qJD(4), t45 + (t142 + t150) * t95, t121 * t95 - 0.2e1 * t93 * t143, -t159 * t93 + t191, t146 + t190, -t128, t215 + (-t116 * t133 + t209) * qJD(2) - t51 * qJD(4) + t4 * qJD(5), t37 * t159 + t210 + t6 * qJD(5) + (qJD(2) * t133 + t173) * t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t233, t34 * qJD(2) + t227, t233, 0, 0, t20 * qJD(2) + t227, 0, 0, 0, 0, 0, t189, -t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, t153, t176, t39 * qJD(2) - t154, 0, 0, 0, 0, 0, t145 - t182, t168 * t225 + t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, -t122, t139 * t51, -t127, -t78, t4 * qJD(2) - t14 * qJD(5) + t204, t6 * qJD(2) + t13 * qJD(5) + t205; 0, 0, 0, -t144, -t162, 0, 0, 0, t156, t155, -t84 - t187, t174 - t186, 0, t35 * qJD(3) - t198, 0, t84 - t197, -t174 - t196, -t19 * qJD(3) + t57 * qJD(4) - t206, -t95 * t150 + t45, -0.2e1 * t116 * t127, -t149 - t191, -t148 - t190, t128, t3 * qJD(5) - t114 * t174 - t215, -t51 * qJD(3) + t5 * qJD(5) - t210; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t105 * qJD(4), -t143, t102 * qJD(5), 0, 0, 0, t105 * t163 + t166, -t105 * t164 + t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172, t175, 0, t188, 0, t172, -t175, -t193, 0, 0, 0, 0, 0, -t147, -t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t131, 0, 0, 0, 0, 0, t161, t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, -t121, -t132, t66, t177, -t104 * t164 - t125, -t104 * t163 - t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, -t82, -t233, -t35 * qJD(2) - t227, -t233, -t171, t82, t19 * qJD(2) - t170 - t227, 0, 0, 0, 0, 0, t146 - t148 - t189, t182 + t184 + t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, -t175, 0, -t188, 0, -t172, t175, t193, 0, 0, 0, 0, 0, t147, t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t153, -t176, -t57 * qJD(2) + t154 + t84, 0, 0, 0, 0, 0, -t145 - t184, (-qJD(5) * t95 - t176) * t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t131, 0, 0, 0, 0, 0, -t161, -t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t236, -t139 * t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123, t122, (-t168 * t93 + t161) * t95, (t147 + t159) * t95, -t78, -t3 * qJD(2) + t48 * qJD(4) + t116 * t84 - t204, -t5 * qJD(2) - t48 * qJD(3) + t165 * t95 - t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, t121, t95 * t169, t77, -t177, t125, t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185, t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t10;
