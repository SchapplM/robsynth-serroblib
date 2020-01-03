% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:02:00
% EndTime: 2019-12-31 21:02:08
% DurationCPUTime: 2.56s
% Computational Cost: add. (3103->311), mult. (7795->437), div. (0->0), fcn. (5020->6), ass. (0->164)
t142 = sin(qJ(3));
t186 = qJD(3) * t142;
t145 = cos(qJ(2));
t190 = qJD(1) * t145;
t232 = -t142 * t190 + t186;
t140 = sin(pkin(8));
t141 = cos(pkin(8));
t144 = cos(qJ(3));
t99 = t140 * t144 + t141 * t142;
t89 = t99 * qJD(3);
t218 = t99 * t190 - t89;
t184 = qJD(3) * t144;
t204 = t141 * t144;
t217 = t232 * t140 - t141 * t184 + t190 * t204;
t125 = -qJD(3) + t190;
t143 = sin(qJ(2));
t191 = qJD(1) * t143;
t174 = t142 * t191;
t181 = t144 * qJD(2);
t103 = t174 - t181;
t189 = qJD(2) * t142;
t105 = t144 * t191 + t189;
t61 = t141 * t103 + t105 * t140;
t231 = t125 * t61;
t187 = qJD(2) * t145;
t173 = t142 * t187;
t230 = t143 * t184 + t173;
t155 = -t103 * t140 + t141 * t105;
t229 = t155 ^ 2;
t180 = qJD(1) * qJD(2);
t228 = -0.2e1 * t180;
t223 = -qJ(4) - pkin(7);
t163 = qJD(3) * t223;
t150 = -qJD(4) * t142 + t144 * t163;
t200 = t144 * t145;
t152 = pkin(3) * t143 - qJ(4) * t200;
t157 = pkin(2) * t143 - pkin(7) * t145;
t106 = t157 * qJD(1);
t211 = pkin(6) * t174 + t144 * t106;
t55 = qJD(1) * t152 + t211;
t201 = t143 * t144;
t202 = t142 * t145;
t91 = t142 * t106;
t65 = t91 + (-pkin(6) * t201 - qJ(4) * t202) * qJD(1);
t183 = qJD(4) * t144;
t86 = t142 * t163 + t183;
t219 = (t150 - t55) * t141 + (t65 - t86) * t140;
t113 = -qJD(2) * pkin(2) + pkin(6) * t191;
t74 = pkin(3) * t103 + qJD(4) + t113;
t23 = pkin(4) * t61 - qJ(5) * t155 + t74;
t227 = t23 * t155;
t134 = pkin(6) * t190;
t226 = t232 * pkin(3) - t134;
t179 = qJD(2) * qJD(3);
t76 = qJD(1) * t230 + t142 * t179;
t225 = pkin(6) * t142;
t114 = qJD(2) * pkin(7) + t134;
t110 = -pkin(2) * t145 - pkin(7) * t143 - pkin(1);
t95 = t110 * qJD(1);
t213 = t142 * t95;
t69 = t114 * t144 + t213;
t49 = -qJ(4) * t103 + t69;
t46 = t141 * t49;
t68 = -t114 * t142 + t144 * t95;
t48 = -qJ(4) * t105 + t68;
t21 = t140 * t48 + t46;
t224 = t21 * t155;
t164 = t143 * t180;
t159 = pkin(6) * t164;
t107 = t157 * qJD(2);
t96 = qJD(1) * t107;
t212 = -t142 * t159 - t144 * t96;
t149 = -qJD(3) * t69 - t212;
t166 = t145 * t180;
t185 = qJD(3) * t143;
t172 = t142 * t185;
t75 = -qJD(1) * t172 + (t166 + t179) * t144;
t15 = pkin(3) * t164 - qJ(4) * t75 - qJD(4) * t105 + t149;
t153 = -t114 * t186 + t142 * t96 + t95 * t184;
t148 = -t144 * t159 + t153;
t20 = -qJ(4) * t76 - qJD(4) * t103 + t148;
t3 = -t140 * t20 + t141 * t15;
t4 = t140 * t15 + t141 * t20;
t30 = t140 * t55 + t141 * t65;
t25 = qJ(5) * t191 + t30;
t53 = t140 * t150 + t141 * t86;
t222 = t25 - t53;
t221 = pkin(4) * t191 - t219;
t128 = pkin(6) * t200;
t188 = qJD(2) * t143;
t210 = t144 * t107 + t188 * t225;
t31 = -t143 * t183 + t152 * qJD(2) + (-t128 + (qJ(4) * t143 - t110) * t142) * qJD(3) + t210;
t216 = t142 * t107 + t110 * t184;
t36 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t201 + (-qJD(4) * t143 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t145) * t142 + t216;
t9 = t140 * t31 + t141 * t36;
t220 = t218 * pkin(4) - t217 * qJ(5) + qJD(5) * t99 - t226;
t44 = -pkin(3) * t125 + t48;
t17 = t140 * t44 + t46;
t101 = t144 * t110;
t66 = -qJ(4) * t201 + t101 + (-pkin(3) - t225) * t145;
t194 = t142 * t110 + t128;
t203 = t142 * t143;
t70 = -qJ(4) * t203 + t194;
t38 = t140 * t66 + t141 * t70;
t215 = t140 * t49;
t214 = t142 * t75;
t209 = t103 * t125;
t208 = t105 * t125;
t207 = t113 * t142;
t206 = t113 * t144;
t205 = t125 * t144;
t147 = qJD(1) ^ 2;
t199 = t145 * t147;
t146 = qJD(2) ^ 2;
t198 = t146 * t143;
t197 = t146 * t145;
t22 = t141 * t48 - t215;
t196 = qJD(5) - t22;
t193 = pkin(3) * t203 + t143 * pkin(6);
t138 = t143 ^ 2;
t192 = -t145 ^ 2 + t138;
t182 = t113 * qJD(3);
t178 = qJ(5) * t164 + t4;
t176 = t230 * pkin(3) + pkin(6) * t187;
t175 = -pkin(3) * t144 - pkin(2);
t170 = t125 * t184;
t168 = t145 * t181;
t59 = pkin(3) * t76 + pkin(6) * t166;
t167 = t223 * t142;
t42 = t140 * t75 + t141 * t76;
t162 = pkin(1) * t228;
t161 = t103 + t181;
t160 = -t105 + t189;
t43 = -t140 * t76 + t141 * t75;
t112 = t223 * t144;
t72 = -t112 * t140 - t141 * t167;
t73 = -t141 * t112 + t140 * t167;
t158 = -t73 * t42 + t43 * t72 - t53 * t61;
t156 = -t61 ^ 2 - t229;
t8 = -t140 * t36 + t141 * t31;
t16 = t141 * t44 - t215;
t37 = -t140 * t70 + t141 * t66;
t98 = t140 * t142 - t204;
t154 = qJD(1) * t138 - t125 * t145;
t2 = -pkin(4) * t164 - t3;
t5 = pkin(4) * t42 - qJ(5) * t43 - qJD(5) * t155 + t59;
t131 = -pkin(3) * t141 - pkin(4);
t129 = pkin(3) * t140 + qJ(5);
t82 = -t140 * t203 + t141 * t201;
t81 = t99 * t143;
t57 = pkin(4) * t98 - qJ(5) * t99 + t175;
t51 = t140 * t173 - t141 * t168 + t143 * t89;
t50 = t185 * t98 - t187 * t99;
t45 = pkin(4) * t81 - qJ(5) * t82 + t193;
t35 = pkin(4) * t145 - t37;
t34 = -qJ(5) * t145 + t38;
t24 = pkin(3) * t105 + pkin(4) * t155 + qJ(5) * t61;
t12 = -qJ(5) * t125 + t17;
t11 = pkin(4) * t125 + qJD(5) - t16;
t10 = -pkin(4) * t50 + qJ(5) * t51 - qJD(5) * t82 + t176;
t7 = -pkin(4) * t188 - t8;
t6 = qJ(5) * t188 - qJD(5) * t145 + t9;
t1 = -qJD(5) * t125 + t178;
t13 = [0, 0, 0, 0.2e1 * t145 * t164, t192 * t228, t197, -t198, 0, -pkin(6) * t197 + t143 * t162, pkin(6) * t198 + t145 * t162, t75 * t201 + (t168 - t172) * t105, (-t103 * t144 - t105 * t142) * t187 + (-t214 - t144 * t76 + (t103 * t142 - t105 * t144) * qJD(3)) * t143, t125 * t172 - t145 * t75 + (t105 * t143 + t144 * t154) * qJD(2), t143 * t170 + t145 * t76 + (-t103 * t143 - t142 * t154) * qJD(2), (-t125 - t190) * t188, -(-t110 * t186 + t210) * t125 + (t144 * t182 + pkin(6) * t76 + (qJD(1) * t101 + t68) * qJD(2)) * t143 + ((pkin(6) * t103 + t207) * qJD(2) + (t213 + (pkin(6) * t125 + t114) * t144) * qJD(3) + t212) * t145, (-pkin(6) * t145 * t186 + t216) * t125 + t153 * t145 + (pkin(6) * t75 - t142 * t182) * t143 + ((pkin(6) * t105 + t206) * t145 + (-pkin(6) * t205 - qJD(1) * t194 - t69) * t143) * qJD(2), -t155 * t8 + t16 * t51 + t17 * t50 - t3 * t82 - t37 * t43 - t38 * t42 - t4 * t81 - t61 * t9, t16 * t8 + t17 * t9 + t176 * t74 + t193 * t59 + t3 * t37 + t4 * t38, t10 * t61 + t125 * t7 + t145 * t2 - t23 * t50 + t42 * t45 + t5 * t81 + (-qJD(1) * t35 - t11) * t188, -t1 * t81 - t11 * t51 + t12 * t50 + t155 * t7 + t2 * t82 - t34 * t42 + t35 * t43 - t6 * t61, -t1 * t145 - t10 * t155 - t125 * t6 + t23 * t51 - t43 * t45 - t5 * t82 + (qJD(1) * t34 + t12) * t188, t1 * t34 + t10 * t23 + t11 * t7 + t12 * t6 + t2 * t35 + t45 * t5; 0, 0, 0, -t143 * t199, t192 * t147, 0, 0, 0, t147 * pkin(1) * t143, pkin(1) * t199, -t105 * t205 + t214, (t75 + t209) * t144 + (-t76 + t208) * t142, -t170 + (t125 * t200 + t143 * t160) * qJD(1), t125 * t186 + (-t125 * t202 + t143 * t161) * qJD(1), t125 * t191, -pkin(2) * t76 + t211 * t125 + (pkin(7) * t205 + t207) * qJD(3) + ((-pkin(7) * t189 - t68) * t143 + (-pkin(6) * t161 - t207) * t145) * qJD(1), -pkin(2) * t75 - t91 * t125 + (-pkin(7) * t125 * t142 + t206) * qJD(3) + (-t113 * t200 + (-pkin(7) * t181 + t69) * t143 + (t125 * t201 + t145 * t160) * pkin(6)) * qJD(1), -t155 * t219 + t217 * t16 + t218 * t17 - t3 * t99 + t30 * t61 - t4 * t98 + t158, t4 * t73 - t3 * t72 + t59 * t175 + t226 * t74 + (t53 - t30) * t17 + t219 * t16, t42 * t57 + t5 * t98 - t220 * t61 - t218 * t23 + t221 * t125 + (-qJD(2) * t72 + t11) * t191, -t1 * t98 - t217 * t11 + t218 * t12 + t155 * t221 + t2 * t99 + t25 * t61 + t158, -t43 * t57 - t5 * t99 + t220 * t155 + t217 * t23 + t222 * t125 + (qJD(2) * t73 - t12) * t191, t1 * t73 + t221 * t11 - t222 * t12 + t2 * t72 - t220 * t23 + t5 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105 * t103, -t103 ^ 2 + t105 ^ 2, t75 - t209, -t208 - t76, t164, -t105 * t113 - t125 * t69 + t149, t103 * t113 - t125 * t68 - t148, t17 * t155 - t224 + (-t140 * t42 - t141 * t43) * pkin(3) + (-t16 + t22) * t61, t16 * t21 - t17 * t22 + (-t105 * t74 + t140 * t4 + t141 * t3) * pkin(3), -t125 * t21 - t227 - t24 * t61 + (pkin(4) - t131) * t164 + t3, t12 * t155 - t129 * t42 + t131 * t43 - t224 + (t11 - t196) * t61, t129 * t164 - t23 * t61 + t24 * t155 + (-0.2e1 * qJD(5) + t22) * t125 + t178, t1 * t129 - t11 * t21 + t12 * t196 + t131 * t2 - t23 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, t155 * t16 + t17 * t61 + t59, -t125 * t155 + t42, t156, -t43 - t231, -t11 * t155 + t12 * t61 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155 * t61 - t164, t43 - t231, -t125 ^ 2 - t229, t12 * t125 + t2 + t227;];
tauc_reg = t13;
