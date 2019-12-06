% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRR1
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:25:33
% EndTime: 2019-12-05 18:25:42
% DurationCPUTime: 1.92s
% Computational Cost: add. (5919->263), mult. (13036->356), div. (0->0), fcn. (8355->8), ass. (0->183)
t154 = sin(qJ(5));
t155 = sin(qJ(4));
t160 = cos(qJ(2));
t156 = sin(qJ(2));
t159 = cos(qJ(4));
t197 = t156 * t159;
t119 = (t160 * t155 + t197) * qJD(1);
t151 = qJD(2) + qJD(4);
t158 = cos(qJ(5));
t101 = t158 * t119 + t154 * t151;
t99 = t154 * t119 - t158 * t151;
t81 = t101 * t99;
t146 = t156 * qJDD(1);
t192 = qJD(1) * qJD(2);
t182 = t160 * t192;
t128 = t146 + t182;
t147 = t160 * qJDD(1);
t183 = t156 * t192;
t129 = t147 - t183;
t179 = t155 * t128 - t159 * t129;
t83 = -t119 * qJD(4) - t179;
t82 = qJDD(5) - t83;
t225 = -t81 + t82;
t230 = t154 * t225;
t150 = qJDD(2) + qJDD(4);
t195 = qJD(1) * t156;
t117 = -t159 * t160 * qJD(1) + t155 * t195;
t95 = t119 * t117;
t223 = -t95 + t150;
t229 = t155 * t223;
t228 = t158 * t225;
t227 = t159 * t223;
t162 = qJD(1) ^ 2;
t215 = pkin(3) + qJ(3);
t218 = pkin(1) + pkin(2);
t226 = -t215 * t128 + (t218 * t162 * t156 + t215 * t192 - g(3)) * t160;
t109 = t151 * t117;
t84 = -t117 * qJD(4) + t159 * t128 + t155 * t129;
t224 = -t84 + t109;
t171 = qJD(2) * pkin(2) - pkin(3) * t195;
t133 = qJD(2) * pkin(1) - qJ(3) * t195;
t157 = sin(qJ(1));
t217 = cos(qJ(1));
t136 = t157 * g(1) - t217 * g(2);
t173 = -t133 * t195 - qJDD(3) + t136;
t153 = t160 ^ 2;
t200 = t153 * t162;
t68 = t218 * t129 - t171 * t195 + t215 * t200 + t173;
t163 = t224 * pkin(4) - t68;
t191 = qJD(1) * qJD(3);
t187 = -0.2e1 * t191;
t137 = t217 * g(1) + t157 * g(2);
t198 = t156 * t137;
t206 = qJDD(2) * pkin(1);
t184 = t160 * t191;
t143 = 0.2e1 * t184;
t113 = -t156 * g(3) - t160 * t137;
t186 = -t129 * qJ(3) - t113;
t64 = t129 * pkin(3) + t143 - t218 * t200 + (-t133 - t171) * qJD(2) - t186;
t35 = t159 * t64 + (qJDD(2) * pkin(2) + t156 * t187 + t198 + t206 + t226) * t155;
t92 = t95 + t150;
t26 = t92 * pkin(4) + t35;
t11 = t154 * t26 - t158 * t163;
t12 = t154 * t163 + t158 * t26;
t6 = t154 * t11 + t158 * t12;
t161 = qJD(2) ^ 2;
t222 = -t161 + t200;
t219 = t151 ^ 2;
t220 = t119 ^ 2;
t106 = -t219 - t220;
t115 = qJD(5) + t117;
t180 = -t158 * t150 + t154 * t84;
t50 = (qJD(5) - t115) * t101 + t180;
t177 = t137 + t187;
t196 = t160 * t162;
t221 = -(pkin(1) * t196 + t177) * t156 - t206;
t97 = t99 ^ 2;
t98 = t101 ^ 2;
t111 = t115 ^ 2;
t116 = t117 ^ 2;
t216 = t160 * g(3);
t174 = t177 * t156;
t34 = t155 * t64 - t159 * (t218 * qJDD(2) + t174 + t226);
t29 = t106 * pkin(4) + t34;
t27 = t154 * t29;
t57 = t81 + t82;
t209 = t158 * t57;
t78 = -t98 - t111;
t38 = -t154 * t78 - t209;
t214 = pkin(4) * t38 + t27;
t28 = t158 * t29;
t69 = -t111 - t97;
t33 = t158 * t69 - t230;
t213 = pkin(4) * t33 - t28;
t212 = t154 * t57;
t211 = t155 * t68;
t210 = t155 * t92;
t208 = t159 * t68;
t207 = t159 * t92;
t205 = t115 * t154;
t204 = t115 * t158;
t203 = t151 * t155;
t202 = t151 * t159;
t152 = t156 ^ 2;
t201 = t152 * t162;
t140 = t156 * t196;
t134 = qJDD(2) + t140;
t199 = t156 * t134;
t193 = qJD(5) + t115;
t175 = -t154 * t150 - t158 * t84;
t63 = -t99 * qJD(5) - t175;
t89 = t115 * t99;
t54 = t63 + t89;
t25 = t154 * t54 - t158 * t50;
t190 = pkin(4) * t25 + t6;
t189 = t155 * t81;
t188 = t159 * t81;
t185 = t156 * t215;
t181 = t155 * t34 + t159 * t35;
t178 = -t201 + t222;
t176 = -t128 + t182;
t5 = -t158 * t11 + t154 * t12;
t13 = t155 * t35 - t159 * t34;
t172 = qJD(2) * t133 + t186;
t169 = t143 - t172;
t167 = (-qJD(4) + t151) * t119 - t179;
t166 = -t129 * pkin(1) - t173;
t165 = qJ(3) * t176 - t216;
t135 = qJDD(2) - t140;
t132 = (-t152 - t153) * t162;
t131 = (t152 - t153) * t162;
t130 = t147 - 0.2e1 * t183;
t127 = t146 + 0.2e1 * t182;
t112 = -t198 + t216;
t108 = t219 - t220;
t107 = t116 - t219;
t105 = t199 + t160 * (t161 - t201);
t104 = t160 * t135 + t156 * t222;
t103 = (t128 + t182) * t156;
t102 = (t129 - t183) * t160;
t96 = t160 * t127 + t156 * t130;
t94 = -t116 + t220;
t90 = -t219 - t116;
t88 = -qJ(3) * t200 + t166;
t87 = -t98 + t111;
t86 = t97 - t111;
t85 = -t116 - t220;
t80 = t98 - t97;
t79 = t165 - t221;
t76 = t159 * t106 - t210;
t75 = t109 + t84;
t70 = (qJD(4) + t151) * t119 + t179;
t66 = t155 * t90 + t227;
t65 = t97 + t98;
t62 = -t101 * qJD(5) - t180;
t60 = (t101 * t154 - t158 * t99) * t115;
t59 = (-t101 * t158 - t154 * t99) * t115;
t55 = t193 * t99 + t175;
t53 = t63 - t89;
t51 = -t193 * t101 - t180;
t49 = -t101 * t205 + t158 * t63;
t48 = t101 * t204 + t154 * t63;
t47 = -t154 * t62 + t99 * t204;
t46 = t158 * t62 + t99 * t205;
t44 = t155 * t167 - t159 * t75;
t42 = t158 * t86 - t212;
t41 = -t154 * t87 + t228;
t40 = t154 * t86 + t209;
t39 = t158 * t87 + t230;
t37 = t158 * t78 - t212;
t32 = t154 * t69 + t228;
t24 = -t154 * t53 + t158 * t51;
t23 = -t154 * t50 - t158 * t54;
t22 = t154 * t51 + t158 * t53;
t19 = t155 * t38 + t159 * t55;
t17 = t155 * t33 + t159 * t51;
t15 = t155 * t25 + t159 * t65;
t8 = -pkin(4) * t37 + t28;
t7 = -pkin(4) * t32 + t27;
t4 = pkin(4) * t6;
t2 = t155 * t6 - t159 * t29;
t1 = -pkin(4) * t23 - t5;
t3 = [0, 0, 0, 0, 0, qJDD(1), t136, t137, 0, 0, t103, t96, t105, t102, t104, 0, t160 * t136, -t156 * t136, t156 * t112 + t160 * t113, 0, t103, t96, t105, t102, t104, 0, t160 * ((t129 + t130) * pkin(1) + t173) + (-t160 * t161 - t199) * qJ(3), t156 * (-qJ(3) * t178 + t166) + t160 * (-pkin(1) * t127 - qJ(3) * t135), t160 * (qJ(3) * t147 + (-t132 - t200) * pkin(1) + t169) + (t216 + (-t176 + t146) * qJ(3) + t221) * t156, -t156 * qJ(3) * t79 + t160 * (-pkin(1) * t88 + qJ(3) * (-pkin(1) * t200 + t169)), t156 * (-t119 * t203 + t159 * t84) + t160 * (t119 * t202 + t155 * t84), t156 * (t155 * t224 - t159 * t70) + t160 * (-t155 * t70 - t159 * t224), t156 * (-t155 * t108 + t227) + t160 * (t159 * t108 + t229), t156 * (t117 * t202 - t155 * t83) + t160 * (t117 * t203 + t159 * t83), t156 * (t159 * t107 - t210) + t160 * (t155 * t107 + t207), (t156 * (-t117 * t159 + t119 * t155) + t160 * (-t117 * t155 - t119 * t159)) * t151, t156 * (-t215 * t66 - t211) + t160 * (t208 - t218 * t70 + t215 * (t159 * t90 - t229)), t156 * (-t215 * t76 - t208) + t160 * (-t211 + t215 * (-t155 * t106 - t207) + t218 * t224), t156 * (-t215 * t44 - t13) + t160 * (-t218 * t85 + t215 * (t155 * t75 + t159 * t167) + t181), t160 * (t215 * t181 + t218 * t68) - t13 * t185, t156 * (t159 * t49 + t189) + t160 * (t155 * t49 - t188), t156 * (t155 * t80 + t159 * t24) + t160 * (t155 * t24 - t159 * t80), t156 * (t155 * t54 + t159 * t41) + t160 * (t155 * t41 - t159 * t54), t156 * (t159 * t47 - t189) + t160 * (t155 * t47 + t188), t156 * (-t155 * t50 + t159 * t42) + t160 * (t155 * t42 + t159 * t50), t156 * (t155 * t82 + t159 * t60) + t160 * (t155 * t60 - t159 * t82), t156 * (-t155 * t11 + t159 * t7 - t215 * t17) + t160 * (t159 * t11 + t155 * t7 - t218 * t32 + t215 * (-t155 * t51 + t159 * t33)), t156 * (-t155 * t12 + t159 * t8 - t215 * t19) + t160 * (t159 * t12 + t155 * t8 - t218 * t37 + t215 * (-t155 * t55 + t159 * t38)), t156 * (t159 * t1 - t215 * t15) + t160 * (t155 * t1 - t218 * t23 + t215 * (-t155 * t65 + t159 * t25)), t160 * t215 * (t155 * t29 + t159 * t6) - t2 * t185 + (-pkin(4) * t197 + t160 * (-pkin(4) * t155 - t218)) * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, t131, t146, t140, t147, qJDD(2), -t112, -t113, 0, 0, -t140, t131, t146, t140, t147, qJDD(2), 0.2e1 * pkin(1) * t134 + t165 + t174, pkin(1) * t178 + t172 - 0.2e1 * t184, -pkin(1) * t146, pkin(1) * t79, t95, t94, t75, -t95, t167, t150, t218 * t66 - t34, t218 * t76 - t35, t218 * t44, t218 * t13, t48, t22, t39, t46, t40, t59, t218 * t17 + t213, t218 * t19 + t214, t218 * t15 + t190, t218 * t2 + t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, t127, t132, t88, 0, 0, 0, 0, 0, 0, t70, -t224, t85, -t68, 0, 0, 0, 0, 0, 0, t32, t37, t23, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, t94, t75, -t95, t167, t150, -t34, -t35, 0, 0, t48, t22, t39, t46, t40, t59, t213, t214, t190, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t80, t54, -t81, -t50, t82, -t11, -t12, 0, 0;];
tauJ_reg = t3;
