% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% tauc_reg [6x31]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:36:20
% EndTime: 2019-03-09 02:36:27
% DurationCPUTime: 2.65s
% Computational Cost: add. (3397->322), mult. (7874->446), div. (0->0), fcn. (6018->8), ass. (0->162)
t136 = sin(pkin(10));
t137 = cos(pkin(10));
t141 = sin(qJ(4));
t215 = cos(qJ(4));
t153 = -t215 * t136 - t141 * t137;
t220 = t153 * qJD(1);
t233 = qJD(5) - t220;
t177 = t215 * t137;
t163 = qJD(1) * t177;
t192 = qJD(1) * t136;
t173 = t141 * t192;
t101 = t163 - t173;
t140 = sin(qJ(5));
t143 = cos(qJ(5));
t184 = t143 * qJD(4);
t80 = t140 * t101 - t184;
t236 = t233 * t80;
t139 = sin(qJ(6));
t142 = cos(qJ(6));
t108 = t139 * t143 + t142 * t140;
t218 = qJD(5) + qJD(6);
t76 = t218 * t108;
t210 = -t108 * t220 + t76;
t82 = t140 * qJD(4) + t143 * t101;
t159 = t139 * t80 - t142 * t82;
t30 = t139 * t82 + t142 * t80;
t235 = t159 * t30;
t168 = t143 * t233;
t117 = qJD(4) * t173;
t94 = qJD(4) * t163 - t117;
t205 = t140 * t94;
t234 = -t168 * t233 - t205;
t232 = t159 ^ 2 - t30 ^ 2;
t187 = qJD(6) * t142;
t188 = qJD(6) * t139;
t190 = qJD(5) * t140;
t93 = t220 * qJD(4);
t40 = qJD(5) * t184 - t101 * t190 + t143 * t93;
t206 = t140 * t93;
t41 = t82 * qJD(5) + t206;
t8 = -t139 * t41 + t142 * t40 - t80 * t187 - t82 * t188;
t92 = qJD(6) + t233;
t231 = t30 * t92 + t8;
t138 = -pkin(1) - qJ(3);
t221 = t138 * qJD(1);
t118 = qJD(2) + t221;
t169 = -pkin(7) * qJD(1) + t118;
t95 = t169 * t136;
t96 = t169 * t137;
t56 = t141 * t96 + t215 * t95;
t50 = qJD(4) * pkin(8) + t56;
t135 = qJD(1) * qJ(2);
t129 = qJD(3) + t135;
t113 = pkin(3) * t192 + t129;
t51 = -pkin(4) * t220 - t101 * pkin(8) + t113;
t21 = t140 * t51 + t143 * t50;
t13 = -t80 * pkin(9) + t21;
t11 = t13 * t188;
t225 = -t141 * t95 + t215 * t96;
t49 = -qJD(4) * pkin(4) - t225;
t24 = t80 * pkin(5) + t49;
t230 = t24 * t30 + t11;
t145 = t159 * qJD(6) - t139 * t40 - t142 * t41;
t229 = -t159 * t92 + t145;
t148 = t153 * qJD(3);
t26 = qJD(1) * t148 + qJD(4) * t225;
t134 = qJD(1) * qJD(2);
t52 = t94 * pkin(4) - t93 * pkin(8) + t134;
t46 = t143 * t52;
t146 = -t21 * qJD(5) - t140 * t26 + t46;
t2 = t94 * pkin(5) - t40 * pkin(9) + t146;
t189 = qJD(5) * t143;
t155 = t140 * t52 + t143 * t26 + t51 * t189 - t50 * t190;
t5 = -t41 * pkin(9) + t155;
t180 = -t139 * t5 + t142 * t2;
t20 = -t140 * t50 + t143 * t51;
t12 = -t82 * pkin(9) + t20;
t10 = pkin(5) * t233 + t12;
t203 = t142 * t13;
t4 = t139 * t10 + t203;
t228 = -t4 * qJD(6) + t24 * t159 + t180;
t226 = -t143 * t94 + t190 * t233;
t224 = t139 * t190 + t140 * t188;
t214 = -pkin(7) + t138;
t111 = t214 * t136;
t112 = t214 * t137;
t223 = -t141 * t111 + t215 * t112;
t222 = -qJD(6) * t143 - t189;
t193 = t136 ^ 2 + t137 ^ 2;
t219 = t193 * qJD(3);
t207 = t108 * t94;
t107 = t139 * t140 - t142 * t143;
t211 = (-t218 + t220) * t107;
t217 = -t211 * t92 - t207;
t183 = 0.2e1 * t134;
t216 = pkin(8) + pkin(9);
t69 = t101 * pkin(4) - pkin(8) * t220;
t213 = t140 * t69 + t143 * t225;
t123 = t136 * pkin(3) + qJ(2);
t152 = -t141 * t136 + t177;
t68 = -pkin(4) * t153 - pkin(8) * t152 + t123;
t74 = t215 * t111 + t141 * t112;
t70 = t143 * t74;
t212 = t140 * t68 + t70;
t209 = t101 * t30;
t208 = t101 * t80;
t72 = t107 * t94;
t204 = t140 * t220;
t202 = t159 * t101;
t201 = t40 * t140;
t200 = t82 * t101;
t172 = qJD(4) * t215;
t191 = qJD(4) * t141;
t103 = -t136 * t191 + t137 * t172;
t199 = t92 * t103;
t71 = t94 * t153;
t102 = -t136 * t172 - t137 * t191;
t198 = t102 * t140;
t197 = t152 * t140;
t196 = t152 * t143;
t195 = t143 * t102;
t185 = t103 * qJD(4);
t179 = qJD(5) * t216;
t176 = t152 * t189;
t171 = qJD(6) * t10 + t5;
t167 = qJD(1) * t193;
t166 = -qJD(5) * t153 + qJD(1);
t165 = -t210 * t92 - t72;
t27 = qJD(3) * t101 + t95 * t172 + t96 * t191;
t164 = pkin(8) * qJD(5) * t233 + t27;
t162 = -t56 + (t190 - t204) * pkin(5);
t114 = t216 * t140;
t161 = -pkin(9) * t204 + qJD(6) * t114 + t140 * t179 + t213;
t115 = t216 * t143;
t63 = t143 * t69;
t160 = t101 * pkin(5) + qJD(6) * t115 - t140 * t225 + t63 + (-pkin(9) * t220 + t179) * t143;
t158 = t204 * t233 - t226;
t156 = -pkin(8) * t94 + t233 * t49;
t43 = t223 * qJD(4) + t148;
t66 = t103 * pkin(4) - t102 * pkin(8) + qJD(2);
t154 = t140 * t66 + t143 * t43 + t68 * t189 - t74 * t190;
t151 = -t176 - t198;
t150 = t92 * t107;
t44 = t152 * qJD(3) + t74 * qJD(4);
t144 = qJD(1) ^ 2;
t128 = -t143 * pkin(5) - pkin(4);
t97 = t102 * qJD(4);
t65 = t107 * t152;
t64 = t108 * t152;
t61 = t143 * t68;
t59 = t143 * t66;
t47 = pkin(5) * t197 - t223;
t23 = -t151 * pkin(5) + t44;
t22 = -pkin(9) * t197 + t212;
t19 = -pkin(5) * t153 - pkin(9) * t196 - t140 * t74 + t61;
t17 = t41 * pkin(5) + t27;
t15 = t139 * t195 + (t218 * t196 + t198) * t142 - t224 * t152;
t14 = -t107 * t102 - t152 * t76;
t7 = t151 * pkin(9) + t154;
t6 = -pkin(9) * t195 + t103 * pkin(5) - t140 * t43 + t59 + (-t70 + (pkin(9) * t152 - t68) * t140) * qJD(5);
t3 = t142 * t10 - t139 * t13;
t1 = [0, 0, 0, 0, t183, qJ(2) * t183, t136 * t183, t137 * t183, 0.2e1 * qJD(3) * t167 (t129 + t135) * qJD(2) + (-t118 - t221) * t219, t101 * t102 + t152 * t93, -t101 * t103 + t102 * t220 - t152 * t94 + t153 * t93, t97, -t185, 0, -0.2e1 * qJD(2) * t220 - t44 * qJD(4) + t113 * t103 + t123 * t94, -t43 * qJD(4) + t113 * t102 + t123 * t93 + (qJD(1) * t152 + t101) * qJD(2), t82 * t195 - (-t40 * t143 + t82 * t190) * t152 (-t140 * t82 - t143 * t80) * t102 - (t201 + t143 * t41 + (-t140 * t80 + t143 * t82) * qJD(5)) * t152, t82 * t103 - t152 * t226 - t153 * t40 + t195 * t233, -t233 * t198 - t80 * t103 + t41 * t153 - (t189 * t233 + t205) * t152, t103 * t233 - t71 (-t189 * t74 + t59) * t233 + t61 * t94 - (-t189 * t50 + t46) * t153 + t20 * t103 + t44 * t80 - t223 * t41 + t49 * t176 + ((-qJD(5) * t68 - t43) * t233 - t74 * t94 - (-qJD(5) * t51 - t26) * t153 + t27 * t152 + t49 * t102) * t140, -t154 * t233 - t212 * t94 + t155 * t153 - t21 * t103 + t44 * t82 - t223 * t40 + t49 * t195 - (-t27 * t143 + t190 * t49) * t152, -t14 * t159 - t8 * t65, -t14 * t30 - t145 * t65 + t15 * t159 - t8 * t64, -t103 * t159 + t14 * t92 - t153 * t8 - t65 * t94, -t30 * t103 - t145 * t153 - t15 * t92 - t64 * t94, -t71 + t199 (-t139 * t7 + t142 * t6) * t92 + (-t139 * t22 + t142 * t19) * t94 - t180 * t153 + t3 * t103 + t23 * t30 - t47 * t145 + t17 * t64 + t24 * t15 + ((-t139 * t19 - t142 * t22) * t92 + t4 * t153) * qJD(6), -t4 * t103 - t11 * t153 + t24 * t14 - t17 * t65 - t23 * t159 + t47 * t8 + (-(-qJD(6) * t22 + t6) * t92 - t19 * t94 + t2 * t153) * t139 + (-(qJD(6) * t19 + t7) * t92 - t22 * t94 + t171 * t153) * t142; 0, 0, 0, 0, -t144, -t144 * qJ(2), -t144 * t136, -t144 * t137, 0 (-t129 - t219) * qJD(1), 0, 0, 0, 0, 0, qJD(1) * t220 + t97, -qJD(1) * t101 - t185, 0, 0, 0, 0, 0, t140 * t71 - t102 * t80 - t152 * t41 + (-t103 * t140 - t143 * t166) * t233, t143 * t71 - t102 * t82 - t152 * t40 + (-t103 * t143 + t140 * t166) * t233, 0, 0, 0, 0, 0, -t102 * t30 + t152 * t145 - t108 * t199 + qJD(1) * t150 - ((t222 * t142 + t224) * t92 - t207) * t153, t102 * t159 - t152 * t8 + t103 * t150 + t108 * t92 * qJD(1) - (-(t222 * t139 - t140 * t187 - t142 * t190) * t92 + t72) * t153; 0, 0, 0, 0, 0, 0, 0, 0, -t193 * t144, t118 * t167 + t134, 0, 0, 0, 0, 0, -t117 + (t163 + t101) * qJD(4), 0.2e1 * t93, 0, 0, 0, 0, 0, t158 - t208, -t200 + t234, 0, 0, 0, 0, 0, t165 - t209, t202 + t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101 * t220, t101 ^ 2 - t220 ^ 2, 0, t117 + (-t163 + t101) * qJD(4), 0, t56 * qJD(4) - t113 * t101 - t27 (-qJD(3) - t113) * t220, t168 * t82 + t201 (t40 - t236) * t143 + (-t233 * t82 - t41) * t140, -t200 - t234, t158 + t208, -t233 * t101, -pkin(4) * t41 - t20 * t101 - t56 * t80 - t63 * t233 - t164 * t143 + (t225 * t233 + t156) * t140, -pkin(4) * t40 + t21 * t101 + t164 * t140 + t156 * t143 + t213 * t233 - t56 * t82, t8 * t108 - t159 * t211, -t8 * t107 + t108 * t145 + t159 * t210 - t211 * t30, t202 - t217, t165 + t209, -t92 * t101 (-t142 * t114 - t139 * t115) * t94 - t128 * t145 + t17 * t107 - t3 * t101 + (t139 * t161 - t142 * t160) * t92 + t162 * t30 + t210 * t24 -(-t139 * t114 + t142 * t115) * t94 + t128 * t8 + t17 * t108 + t4 * t101 + (t139 * t160 + t142 * t161) * t92 - t162 * t159 + t211 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82 * t80, -t80 ^ 2 + t82 ^ 2, t40 + t236, -t206 + (-qJD(5) + t233) * t82, t94, t21 * t233 - t49 * t82 + t146, t20 * t233 + t49 * t80 - t155, -t235, t232, t231, t229, t94 -(-t139 * t12 - t203) * t92 + (t142 * t94 - t188 * t92 - t82 * t30) * pkin(5) + t228 (-t13 * t92 - t2) * t139 + (t12 * t92 - t171) * t142 + (-t139 * t94 + t159 * t82 - t187 * t92) * pkin(5) + t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t235, t232, t231, t229, t94, t4 * t92 + t228, -t139 * t2 - t142 * t171 + t3 * t92 + t230;];
tauc_reg  = t1;
