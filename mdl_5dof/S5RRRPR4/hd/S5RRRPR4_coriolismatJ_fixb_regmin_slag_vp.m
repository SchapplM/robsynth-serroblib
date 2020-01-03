% Calculate minimal parameter regressor of coriolis matrix for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x24]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRPR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:42
% EndTime: 2019-12-31 21:11:48
% DurationCPUTime: 1.85s
% Computational Cost: add. (1466->226), mult. (2832->263), div. (0->0), fcn. (2533->6), ass. (0->184)
t151 = cos(qJ(3));
t226 = t151 * qJ(4);
t148 = sin(qJ(3));
t246 = t148 * pkin(3);
t117 = -t226 + t246;
t228 = t117 * t151;
t152 = cos(qJ(2));
t245 = t152 * pkin(1);
t135 = -pkin(2) - t245;
t204 = -t151 * pkin(3) - t148 * qJ(4);
t86 = t135 + t204;
t76 = t86 * t148;
t191 = pkin(2) - t204;
t96 = t191 * t148;
t239 = t76 / 0.2e1 - t96 / 0.2e1;
t264 = t239 - t228;
t192 = qJD(3) - qJD(5);
t193 = qJD(1) + qJD(2);
t147 = sin(qJ(5));
t150 = cos(qJ(5));
t97 = t148 * t147 + t151 * t150;
t99 = -t151 * t147 + t148 * t150;
t33 = t97 ^ 2 - t99 ^ 2;
t263 = t193 * t33;
t194 = t151 * qJD(4);
t262 = t204 * qJD(3) + t194;
t143 = t151 * pkin(4);
t71 = t143 - t86;
t87 = t143 + t191;
t181 = t87 / 0.2e1 + t71 / 0.2e1;
t261 = t181 * t97;
t260 = t181 * t99;
t259 = t192 * t97;
t258 = t192 * t99;
t182 = pkin(2) / 0.2e1 - t135 / 0.2e1;
t255 = t182 * t148;
t145 = t148 ^ 2;
t146 = t151 ^ 2;
t130 = t146 - t145;
t254 = t193 * t130;
t253 = t193 * t97 * t99;
t142 = t148 * pkin(8);
t120 = t148 * pkin(7) - t142;
t121 = (pkin(7) - pkin(8)) * t151;
t252 = t192 * (t147 * t120 + t150 * t121);
t251 = t192 * (t150 * t120 - t147 * t121);
t149 = sin(qJ(2));
t134 = t149 * pkin(1) + pkin(7);
t89 = t148 * t134 - t142;
t90 = (-pkin(8) + t134) * t151;
t250 = t192 * (t147 * t89 + t150 * t90);
t249 = t192 * (-t147 * t90 + t150 * t89);
t248 = -pkin(3) - pkin(4);
t247 = pkin(2) * t151;
t244 = t71 * t97;
t243 = t71 * t99;
t242 = t87 * t97;
t241 = t87 * t99;
t237 = pkin(1) * qJD(1);
t236 = pkin(1) * qJD(2);
t235 = t86 * t117;
t234 = t86 * t151;
t233 = qJD(1) * t71;
t232 = qJD(2) * t87;
t231 = t191 * t151;
t230 = t117 * t191;
t229 = t117 * t148;
t227 = t135 * t151;
t95 = t248 * t148 + t226;
t58 = t95 * t97;
t16 = t58 - t243;
t225 = t16 * qJD(1);
t59 = t95 * t99;
t17 = t59 + t244;
t224 = t17 * qJD(1);
t172 = (t145 + t146) * t152;
t29 = (t134 * t172 + t149 * t86) * pkin(1);
t219 = t29 * qJD(1);
t44 = t229 + t234;
t214 = t44 * qJD(1);
t45 = -t76 + t228;
t213 = t45 * qJD(1);
t88 = pkin(1) * t172;
t212 = t88 * qJD(1);
t211 = t97 * qJD(5);
t210 = t99 * qJD(5);
t183 = t245 / 0.2e1;
t123 = t148 * t183;
t125 = t151 * t183;
t209 = t147 * t123 + t150 * t125;
t184 = -t245 / 0.2e1;
t124 = t148 * t184;
t208 = t150 * t124 + t147 * t125;
t126 = t151 * t184;
t207 = t150 * t123 + t147 * t126;
t206 = t147 * t124 + t150 * t126;
t190 = t149 * t236;
t129 = t148 * t190;
t139 = t145 * qJD(4);
t205 = t139 - t129;
t203 = qJD(1) * t148;
t202 = qJD(2) * t148;
t201 = qJD(3) * qJ(4);
t200 = qJD(4) * t148;
t199 = t147 * qJD(3);
t198 = t147 * qJD(4);
t197 = t148 * qJD(3);
t196 = t150 * qJD(3);
t195 = t150 * qJD(4);
t140 = t151 * qJD(3);
t189 = pkin(7) * t197;
t188 = t149 * t237;
t187 = pkin(7) * t140;
t186 = t97 * t233;
t185 = t99 * t233;
t180 = qJD(1) * t235;
t179 = t86 * t203;
t178 = -t191 / 0.2e1 + t86 / 0.2e1;
t177 = t135 * t203;
t176 = qJD(1) * t227;
t175 = t134 * t197;
t174 = t134 * t140;
t173 = pkin(1) * t193;
t171 = t148 * t193;
t170 = t151 * t190;
t169 = t149 * t173;
t1 = t208 - t58 + t260;
t18 = t58 - t241;
t168 = t1 * qJD(1) - t18 * qJD(2);
t19 = t59 + t242;
t2 = t209 - t59 - t261;
t167 = t2 * qJD(1) - t19 * qJD(2);
t25 = t124 - t264;
t55 = t96 + t228;
t166 = t25 * qJD(1) + t55 * qJD(2);
t26 = t178 * t151 + t125 + t229;
t54 = t229 - t231;
t165 = t26 * qJD(1) + t54 * qJD(2);
t164 = qJD(3) * t117 - t200;
t131 = t148 * t194;
t163 = t131 - t170;
t162 = -t244 / 0.2e1 - t242 / 0.2e1;
t161 = t243 / 0.2e1 + t241 / 0.2e1;
t61 = t124 + t255;
t160 = pkin(2) * t202 + t61 * qJD(1);
t62 = t182 * t151 + t126;
t159 = t62 * qJD(1) + qJD(2) * t247;
t10 = t207 - t260;
t158 = t10 * qJD(1) - t99 * t232;
t11 = t206 + t261;
t157 = t11 * qJD(1) + t97 * t232;
t154 = (t226 / 0.2e1 - t246 / 0.2e1) * t245;
t14 = -t178 * t117 + t154;
t156 = t14 * qJD(1) + qJD(2) * t230;
t30 = t123 + t239;
t155 = t30 * qJD(1) - t191 * t202;
t132 = t148 * t140;
t128 = t151 * t188;
t127 = t148 * t188;
t119 = t130 * qJD(3);
t116 = t192 * t150;
t115 = t193 * t145;
t114 = t192 * t147;
t102 = t150 * qJ(4) + t147 * t248;
t101 = t147 * qJ(4) - t150 * t248;
t82 = t151 * t171;
t81 = t88 * qJD(2);
t78 = t99 * t200;
t77 = t97 * t200;
t64 = -t247 / 0.2e1 + t227 / 0.2e1 + t126;
t63 = t124 - t255;
t39 = t99 * t171;
t38 = t97 * t171;
t31 = t123 - t239;
t28 = t124 + t264;
t27 = -t229 - t234 / 0.2e1 + t125 + t231 / 0.2e1;
t20 = t97 * t258;
t15 = -t230 / 0.2e1 + t235 / 0.2e1 + t154;
t13 = t161 + t207;
t12 = t162 + t206;
t5 = t192 * t33;
t4 = t59 - t162 + t209;
t3 = t58 - t161 + t208;
t6 = [0, 0, 0, 0, -t190, -t152 * t236, t132, t119, 0, 0, 0, t135 * t197 - t170, t135 * t140 + t129, -t45 * qJD(3) + t163, t81, -t44 * qJD(3) + t205, t29 * qJD(2) + t164 * t86, t20, -t5, 0, 0, 0, t16 * qJD(3) - t190 * t97 + t210 * t71 + t77, t17 * qJD(3) - t190 * t99 - t211 * t71 + t78; 0, 0, 0, 0, -t169, -t152 * t173, t132, t119, 0, 0, 0, t63 * qJD(3) - t128 - t170, t64 * qJD(3) + t127 + t129, t28 * qJD(3) - t128 + t163, t81 + t212, t27 * qJD(3) - t127 + t205, t219 + t15 * qJD(3) + t31 * qJD(4) + (pkin(7) * t172 - t149 * t191) * t236, t20, -t5, 0, 0, 0, t3 * qJD(3) + t13 * qJD(5) - t169 * t97 + t77, t4 * qJD(3) + t12 * qJD(5) - t169 * t99 + t78; 0, 0, 0, 0, 0, 0, t82, t254, t140, -t197, 0, t63 * qJD(2) - t174 + t177, t64 * qJD(2) + t175 + t176, t28 * qJD(2) - t174 - t213, t262, t27 * qJD(2) - t175 - t214, t15 * qJD(2) + t134 * t262 + t180, t253, -t263, -t259, -t258, 0, t3 * qJD(2) + t225 - t250, t4 * qJD(2) + t224 - t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t140, t115, t31 * qJD(2) + t174 - t179, 0, 0, 0, 0, 0, t38, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t253, t263, t259, t258, 0, t13 * qJD(2) + t185 + t250, t12 * qJD(2) - t186 + t249; 0, 0, 0, 0, t188, t152 * t237, t132, t119, 0, 0, 0, -t61 * qJD(3) + t128, -t62 * qJD(3) - t127, -t25 * qJD(3) + t128 + t131, -t212, -t26 * qJD(3) + t127 + t139, -t14 * qJD(3) - t30 * qJD(4) - t219, t20, -t5, 0, 0, 0, -t1 * qJD(3) - t10 * qJD(5) + t188 * t97 + t77, -t2 * qJD(3) - t11 * qJD(5) + t188 * t99 + t78; 0, 0, 0, 0, 0, 0, t132, t119, 0, 0, 0, -pkin(2) * t197, -pkin(2) * t140, -t55 * qJD(3) + t131, 0, -t54 * qJD(3) + t139, -t164 * t191, t20, -t5, 0, 0, 0, t18 * qJD(3) + t87 * t210 + t77, t19 * qJD(3) - t87 * t211 + t78; 0, 0, 0, 0, 0, 0, t82, t254, t140, -t197, 0, -t160 - t187, -t159 + t189, -t166 - t187, t262, -t165 - t189, pkin(7) * t262 - t156, t253, -t263, -t259, -t258, 0, -t168 - t252, -t167 - t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t140, t115, -t155 + t187, 0, 0, 0, 0, 0, t38, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t253, t263, t259, t258, 0, -t158 + t252, -t157 + t251; 0, 0, 0, 0, 0, 0, -t82, -t254, 0, 0, 0, t61 * qJD(2) - t177, t62 * qJD(2) - t176, t25 * qJD(2) + t213, 0, t26 * qJD(2) + t214, t14 * qJD(2) - t180, -t253, t263, 0, 0, 0, t1 * qJD(2) - t225, t2 * qJD(2) - t224; 0, 0, 0, 0, 0, 0, -t82, -t254, 0, 0, 0, t160, t159, t166, 0, t165, t156, -t253, t263, 0, 0, 0, t168, t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4), 0, 0, 0, 0, 0, t102 * qJD(5) + t198, -t101 * qJD(5) + t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t201, 0, 0, 0, 0, 0, t199, t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192 * t102, -t192 * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, 0, -t115, t30 * qJD(2) + t179, 0, 0, 0, 0, 0, -t38, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, 0, -t115, t155, 0, 0, 0, 0, 0, -t38, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t201, 0, 0, 0, 0, 0, -t114, -t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, -t263, 0, 0, 0, t10 * qJD(2) - t185, t11 * qJD(2) + t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, -t263, 0, 0, 0, t158, t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102 * qJD(3) - t198, t101 * qJD(3) - t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199, -t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
