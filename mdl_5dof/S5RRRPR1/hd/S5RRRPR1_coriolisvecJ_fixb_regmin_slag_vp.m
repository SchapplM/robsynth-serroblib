% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:50
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:48:58
% EndTime: 2021-01-15 22:49:09
% DurationCPUTime: 2.52s
% Computational Cost: add. (3662->265), mult. (9959->364), div. (0->0), fcn. (7348->8), ass. (0->166)
t160 = sin(qJ(5));
t163 = cos(qJ(5));
t164 = cos(qJ(3));
t165 = cos(qJ(2));
t200 = qJD(1) * t165;
t193 = t164 * t200;
t161 = sin(qJ(3));
t162 = sin(qJ(2));
t201 = qJD(1) * t162;
t194 = t161 * t201;
t120 = -t193 + t194;
t122 = -t161 * t200 - t164 * t201;
t158 = sin(pkin(9));
t159 = cos(pkin(9));
t177 = -t158 * t120 - t159 * t122;
t198 = qJD(5) * t160;
t95 = -t159 * t120 + t158 * t122;
t219 = t163 * t95;
t155 = qJD(2) + qJD(3);
t197 = qJD(1) * qJD(2);
t191 = t165 * t197;
t98 = qJD(3) * t193 - t155 * t194 + t164 * t191;
t134 = t161 * t165 + t164 * t162;
t104 = t155 * t134;
t99 = t104 * qJD(1);
t54 = t158 * t98 + t159 * t99;
t55 = -t158 * t99 + t159 * t98;
t11 = qJD(5) * t219 - t160 * t54 + t163 * t55 - t177 * t198;
t154 = qJD(5) + t155;
t49 = -t160 * t177 + t219;
t217 = t49 * t154;
t4 = t11 - t217;
t46 = t160 * t95 + t163 * t177;
t223 = t46 * t49;
t170 = -qJD(5) * t46 - t160 * t55 - t163 * t54;
t218 = t46 * t154;
t5 = t170 + t218;
t6 = t46 ^ 2 - t49 ^ 2;
t226 = t95 * pkin(8);
t227 = pkin(6) + pkin(7);
t142 = t227 * t165;
t137 = qJD(1) * t142;
t127 = t164 * t137;
t141 = t227 * t162;
t135 = qJD(1) * t141;
t221 = qJD(2) * pkin(2);
t129 = -t135 + t221;
t176 = -t161 * t129 - t127;
t212 = t120 * qJ(4);
t78 = -t176 - t212;
t220 = t159 * t78;
t117 = t122 * qJ(4);
t123 = t161 * t137;
t187 = t164 * t129 - t123;
t77 = t117 + t187;
t66 = t155 * pkin(3) + t77;
t33 = t158 * t66 + t220;
t19 = t33 + t226;
t151 = -t165 * pkin(2) - pkin(1);
t140 = t151 * qJD(1);
t105 = t120 * pkin(3) + qJD(4) + t140;
t61 = -pkin(4) * t95 + t105;
t192 = t19 * t198 - t61 * t49;
t195 = qJD(2) * t227;
t181 = qJD(1) * t195;
t131 = t165 * t181;
t199 = qJD(3) * t161;
t185 = -t161 * t131 - t137 * t199;
t130 = t162 * t181;
t230 = (qJD(3) * t129 - t130) * t164;
t27 = -t99 * qJ(4) - t120 * qJD(4) + t185 + t230;
t186 = t161 * t130 - t164 * t131;
t169 = t176 * qJD(3) + t186;
t28 = -t98 * qJ(4) + t122 * qJD(4) + t169;
t9 = -t158 * t27 + t159 * t28;
t2 = -t55 * pkin(8) + t9;
t10 = t158 * t28 + t159 * t27;
t3 = -t54 * pkin(8) + t10;
t174 = -t160 * t3 + t163 * t2 - t61 * t46;
t233 = -0.2e1 * t197;
t89 = t177 * pkin(8);
t208 = t159 * t161;
t222 = pkin(2) * qJD(3);
t184 = t161 * t135 - t127;
t79 = t184 + t212;
t203 = -t164 * t135 - t123;
t80 = t117 + t203;
t215 = t158 * t80 - t159 * t79 + (-t158 * t164 - t208) * t222;
t209 = t158 * t161;
t213 = -t158 * t79 - t159 * t80 + (t159 * t164 - t209) * t222;
t229 = qJD(1) * t134;
t228 = qJD(5) - t154;
t225 = pkin(3) * t158;
t224 = t122 * pkin(3);
t133 = t161 * t162 - t164 * t165;
t136 = t162 * t195;
t138 = t165 * t195;
t207 = t164 * t141;
t171 = -qJD(3) * t207 - t164 * t136 - t161 * t138 - t142 * t199;
t40 = -t104 * qJ(4) - t133 * qJD(4) + t171;
t103 = t155 * t133;
t175 = t161 * t141 - t164 * t142;
t168 = t175 * qJD(3) + t161 * t136 - t164 * t138;
t41 = t103 * qJ(4) - t134 * qJD(4) + t168;
t16 = t158 * t41 + t159 * t40;
t69 = t158 * t78;
t35 = t159 * t77 - t69;
t90 = -t134 * qJ(4) - t161 * t142 - t207;
t91 = -t133 * qJ(4) - t175;
t45 = t158 * t90 + t159 * t91;
t216 = t226 + t215;
t214 = -t89 - t213;
t211 = t122 * t120;
t210 = t140 * t122;
t167 = qJD(1) ^ 2;
t206 = t165 * t167;
t166 = qJD(2) ^ 2;
t205 = t166 * t162;
t204 = t166 * t165;
t202 = t162 ^ 2 - t165 ^ 2;
t153 = t162 * t221;
t152 = pkin(2) * t201;
t82 = t99 * pkin(3) + qJD(2) * t152;
t190 = -pkin(2) * t155 - t129;
t97 = t104 * pkin(3) + t153;
t15 = -t158 * t40 + t159 * t41;
t32 = t159 * t66 - t69;
t34 = -t158 * t77 - t220;
t44 = -t158 * t91 + t159 * t90;
t188 = pkin(1) * t233;
t182 = -t105 * t95 - t10;
t150 = t164 * pkin(2) + pkin(3);
t115 = -pkin(2) * t209 + t159 * t150;
t65 = pkin(4) * t177 - t224;
t180 = t177 * t33 + t32 * t95;
t17 = t155 * pkin(4) + t32 - t89;
t179 = -t160 * t17 - t163 * t19;
t101 = t159 * t133 + t158 * t134;
t102 = -t158 * t133 + t159 * t134;
t56 = t163 * t101 + t160 * t102;
t57 = -t160 * t101 + t163 * t102;
t107 = t133 * pkin(3) + t151;
t173 = -t105 * t177 + t9;
t172 = t140 * t120 - t185;
t149 = t159 * pkin(3) + pkin(4);
t116 = pkin(2) * t208 + t158 * t150;
t110 = pkin(4) + t115;
t106 = t152 - t224;
t81 = -t120 ^ 2 + t122 ^ 2;
t73 = t101 * pkin(4) + t107;
t68 = (-t122 - t229) * t155;
t67 = t120 * t155 + t98;
t62 = t152 + t65;
t60 = -t159 * t103 - t158 * t104;
t59 = -t158 * t103 + t159 * t104;
t39 = t59 * pkin(4) + t97;
t31 = -t101 * pkin(8) + t45;
t30 = -t102 * pkin(8) + t44;
t29 = t54 * pkin(4) + t82;
t21 = -t89 + t35;
t20 = t34 - t226;
t14 = t57 * qJD(5) + t160 * t60 + t163 * t59;
t13 = -t56 * qJD(5) - t160 * t59 + t163 * t60;
t8 = -t59 * pkin(8) + t16;
t7 = -t60 * pkin(8) + t15;
t1 = [0, 0, 0, 0.2e1 * t162 * t191, t202 * t233, t204, -t205, 0, -pkin(6) * t204 + t162 * t188, pkin(6) * t205 + t165 * t188, t122 * t103 + t98 * t134, t103 * t120 + t122 * t104 - t98 * t133 - t134 * t99, -t103 * t155, -t104 * t155, 0, t151 * t99 + t140 * t104 + t168 * t155 + (qJD(1) * t133 + t120) * t153, t151 * t98 - t140 * t103 - t171 * t155 + (-t122 + t229) * t153, t82 * t101 + t105 * t59 + t107 * t54 + t15 * t155 - t95 * t97, t82 * t102 + t105 * t60 + t107 * t55 - t16 * t155 + t177 * t97, -t10 * t101 - t9 * t102 - t15 * t177 + t16 * t95 - t32 * t60 - t33 * t59 - t44 * t55 - t45 * t54, t10 * t45 + t105 * t97 + t82 * t107 + t32 * t15 + t33 * t16 + t9 * t44, t11 * t57 + t13 * t46, -t11 * t56 + t13 * t49 - t14 * t46 + t170 * t57, t13 * t154, -t14 * t154, 0, -t39 * t49 - t73 * t170 + t29 * t56 + t61 * t14 + (-t160 * t8 + t163 * t7 + (-t160 * t30 - t163 * t31) * qJD(5)) * t154, t39 * t46 + t73 * t11 + t29 * t57 + t61 * t13 - (t160 * t7 + t163 * t8 + (-t160 * t31 + t163 * t30) * qJD(5)) * t154; 0, 0, 0, -t162 * t206, t202 * t167, 0, 0, 0, t167 * pkin(1) * t162, pkin(1) * t206, -t211, t81, t67, t68, 0, -t120 * t152 + t210 - t184 * t155 + (t161 * t190 - t127) * qJD(3) + t186, t122 * t152 + t203 * t155 + (qJD(3) * t190 + t130) * t164 + t172, t106 * t95 + t215 * t155 + t173, -t106 * t177 - t213 * t155 + t182, -t115 * t55 - t116 * t54 - t177 * t215 + t213 * t95 + t180, t10 * t116 - t105 * t106 + t9 * t115 + t213 * t33 + t215 * t32, -t223, t6, t4, t5, 0, t62 * t49 + (t214 * t160 + t216 * t163) * t154 + ((-t110 * t160 - t116 * t163) * t154 + t179) * qJD(5) + t174, -t62 * t46 + (-t2 + (qJD(5) * t116 - t216) * t154) * t160 + (-qJD(5) * t17 - t3 + (-qJD(5) * t110 + t214) * t154) * t163 + t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t211, t81, t67, t68, 0, -t155 * t176 + t169 + t210, t155 * t187 + t172 - t230, -t34 * t155 - t224 * t95 + t173, t35 * t155 + t177 * t224 + t182, t34 * t177 - t35 * t95 + (-t158 * t54 - t159 * t55) * pkin(3) + t180, -t32 * t34 - t33 * t35 + (t10 * t158 + t105 * t122 + t159 * t9) * pkin(3), -t223, t6, t4, t5, 0, t65 * t49 - (-t160 * t21 + t163 * t20) * t154 + ((-t149 * t160 - t163 * t225) * t154 + t179) * qJD(5) + t174, -t163 * t3 - t160 * t2 - t65 * t46 + (t160 * t20 + t163 * t21) * t154 + (-(t149 * t163 - t160 * t225) * t154 - t163 * t17) * qJD(5) + t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155 * t177 + t54, t95 * t155 + t55, -t177 ^ 2 - t95 ^ 2, t177 * t32 - t33 * t95 + t82, 0, 0, 0, 0, 0, -t170 + t218, t11 + t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t223, t6, t4, t5, 0, t228 * t179 + t174, (-t19 * t154 - t2) * t160 + (-t228 * t17 - t3) * t163 + t192;];
tauc_reg = t1;
