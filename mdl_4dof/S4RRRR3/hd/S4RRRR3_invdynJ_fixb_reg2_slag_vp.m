% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:39
% EndTime: 2019-12-31 17:24:43
% DurationCPUTime: 1.97s
% Computational Cost: add. (3218->289), mult. (7860->395), div. (0->0), fcn. (5463->12), ass. (0->150)
t147 = qJ(2) + qJ(3);
t139 = sin(t147);
t140 = cos(t147);
t151 = sin(qJ(1));
t155 = cos(qJ(1));
t177 = g(1) * t155 + g(2) * t151;
t226 = -g(3) * t140 + t177 * t139;
t148 = sin(qJ(4));
t152 = cos(qJ(4));
t153 = cos(qJ(3));
t154 = cos(qJ(2));
t202 = t153 * t154;
t187 = qJD(1) * t202;
t149 = sin(qJ(3));
t150 = sin(qJ(2));
t199 = qJD(1) * t150;
t188 = t149 * t199;
t87 = -t187 + t188;
t99 = t149 * t154 + t150 * t153;
t89 = t99 * qJD(1);
t172 = t148 * t87 - t152 * t89;
t47 = t148 * t89 + t152 * t87;
t225 = t47 * t172;
t194 = qJD(1) * qJD(2);
t184 = t154 * t194;
t193 = t150 * qJDD(1);
t223 = -t184 - t193;
t9 = t172 ^ 2 - t47 ^ 2;
t195 = qJD(4) * t152;
t196 = qJD(4) * t148;
t143 = qJD(2) + qJD(3);
t204 = t149 * t150;
t174 = t143 * t204;
t192 = t154 * qJDD(1);
t179 = -qJD(3) * t187 - t149 * t192 + t223 * t153;
t39 = qJD(1) * t174 + t179;
t173 = t149 * t193 - t153 * t192;
t66 = t143 * t99;
t40 = t66 * qJD(1) + t173;
t11 = t148 * t40 + t152 * t39 + t87 * t195 + t89 * t196;
t138 = qJD(4) + t143;
t5 = t138 * t47 - t11;
t222 = pkin(6) + pkin(5);
t115 = t222 * t154;
t106 = qJD(1) * t115;
t197 = qJD(3) * t153;
t198 = qJD(3) * t149;
t68 = qJDD(2) * pkin(2) + t222 * t223;
t185 = t150 * t194;
t69 = t222 * (-t185 + t192);
t114 = t222 * t150;
t104 = qJD(1) * t114;
t209 = qJD(2) * pkin(2);
t96 = -t104 + t209;
t181 = t106 * t198 - t149 * t68 - t153 * t69 - t96 * t197;
t10 = -pkin(7) * t40 - t181;
t90 = t149 * t106;
t54 = t153 * t96 - t90;
t81 = t89 * pkin(7);
t37 = t54 - t81;
t32 = pkin(3) * t143 + t37;
t221 = pkin(7) * t87;
t94 = t153 * t106;
t55 = t149 * t96 + t94;
t38 = t55 - t221;
t142 = qJDD(2) + qJDD(3);
t22 = -t55 * qJD(3) - t149 * t69 + t153 * t68;
t8 = pkin(3) * t142 + pkin(7) * t39 + t22;
t1 = (qJD(4) * t32 + t10) * t152 + t148 * t8 - t38 * t196;
t141 = qJ(4) + t147;
t133 = sin(t141);
t134 = cos(t141);
t220 = pkin(2) * t154;
t136 = pkin(1) + t220;
t113 = t136 * qJD(1);
t70 = t87 * pkin(3) - t113;
t163 = g(3) * t133 + t177 * t134 + t47 * t70 - t1;
t207 = t152 * t38;
t14 = t148 * t32 + t207;
t2 = -t14 * qJD(4) - t10 * t148 + t152 * t8;
t161 = -g(3) * t134 + t177 * t133 + t172 * t70 + t2;
t165 = t172 * qJD(4) + t148 * t39 - t152 * t40;
t6 = -t138 * t172 + t165;
t72 = -t149 * t114 + t153 * t115;
t216 = g(3) * t154;
t212 = t89 * t87;
t135 = pkin(2) * t153 + pkin(3);
t205 = t148 * t149;
t63 = t104 * t149 - t94;
t42 = t63 + t221;
t64 = -t153 * t104 - t90;
t43 = -t81 + t64;
t211 = -t148 * t42 - t152 * t43 + t135 * t195 + (-t149 * t196 + (t152 * t153 - t205) * qJD(3)) * pkin(2);
t203 = t149 * t152;
t210 = t148 * t43 - t152 * t42 - t135 * t196 + (-t149 * t195 + (-t148 * t153 - t203) * qJD(3)) * pkin(2);
t208 = t148 * t38;
t206 = pkin(5) * qJDD(1);
t145 = t150 ^ 2;
t146 = t154 ^ 2;
t201 = t145 - t146;
t200 = t145 + t146;
t191 = t150 * t209;
t158 = qJD(1) ^ 2;
t190 = t150 * t158 * t154;
t189 = qJD(2) * t222;
t183 = pkin(3) * t140 + t220;
t71 = -t153 * t114 - t115 * t149;
t178 = t150 * t184;
t176 = g(1) * t151 - g(2) * t155;
t13 = t152 * t32 - t208;
t175 = -t13 * t47 - t14 * t172;
t44 = -pkin(7) * t99 + t71;
t98 = -t202 + t204;
t45 = -pkin(7) * t98 + t72;
t25 = -t148 * t45 + t152 * t44;
t26 = t148 * t44 + t152 * t45;
t62 = -t148 * t98 + t152 * t99;
t171 = -0.2e1 * pkin(1) * t194 - pkin(5) * qJDD(2);
t105 = t150 * t189;
t107 = t154 * t189;
t30 = -t153 * t105 - t149 * t107 - t114 * t197 - t115 * t198;
t82 = pkin(2) * t185 - t136 * qJDD(1);
t168 = g(3) * t139 - t113 * t87 + t177 * t140 + t181;
t157 = qJD(2) ^ 2;
t167 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t157 + t176;
t166 = pkin(1) * t158 + t177 - t206;
t31 = -t72 * qJD(3) + t149 * t105 - t153 * t107;
t160 = t113 * t89 + t22 + t226;
t144 = -pkin(7) - t222;
t137 = qJDD(4) + t142;
t103 = pkin(1) + t183;
t84 = pkin(2) * t203 + t135 * t148;
t83 = -pkin(2) * t205 + t135 * t152;
t75 = pkin(3) * t98 - t136;
t73 = pkin(2) * t199 + pkin(3) * t89;
t65 = -qJD(2) * t202 - t154 * t197 + t174;
t61 = t148 * t99 + t152 * t98;
t51 = pkin(3) * t66 + t191;
t41 = -t87 ^ 2 + t89 ^ 2;
t29 = t40 * pkin(3) + t82;
t27 = -t179 + (-t188 + t87) * t143;
t24 = pkin(7) * t65 + t31;
t23 = -pkin(7) * t66 + t30;
t20 = t62 * qJD(4) - t148 * t65 + t152 * t66;
t19 = t148 * t66 + t152 * t65 + t98 * t195 + t99 * t196;
t16 = t152 * t37 - t208;
t15 = -t148 * t37 - t207;
t4 = -t26 * qJD(4) - t148 * t23 + t152 * t24;
t3 = t25 * qJD(4) + t148 * t24 + t152 * t23;
t7 = [0, 0, 0, 0, 0, qJDD(1), t176, t177, 0, 0, qJDD(1) * t145 + 0.2e1 * t178, 0.2e1 * t150 * t192 - 0.2e1 * t201 * t194, qJDD(2) * t150 + t154 * t157, qJDD(1) * t146 - 0.2e1 * t178, qJDD(2) * t154 - t150 * t157, 0, t171 * t150 + t167 * t154, -t167 * t150 + t171 * t154, 0.2e1 * t200 * t206 - t177, -g(1) * (-pkin(1) * t151 + pkin(5) * t155) - g(2) * (pkin(1) * t155 + pkin(5) * t151) + (t200 * pkin(5) ^ 2 + pkin(1) ^ 2) * qJDD(1), -t39 * t99 - t65 * t89, t39 * t98 - t40 * t99 + t65 * t87 - t66 * t89, t142 * t99 - t143 * t65, t40 * t98 + t66 * t87, -t142 * t98 - t143 * t66, 0, -t113 * t66 - t136 * t40 + t140 * t176 + t142 * t71 + t143 * t31 + t87 * t191 + t82 * t98, t113 * t65 + t136 * t39 - t139 * t176 - t142 * t72 - t143 * t30 + t89 * t191 + t82 * t99, t181 * t98 - t22 * t99 - t30 * t87 - t31 * t89 + t39 * t71 - t40 * t72 + t54 * t65 - t55 * t66 - t177, -t181 * t72 + t55 * t30 + t22 * t71 + t54 * t31 - t82 * t136 - t113 * t191 - g(1) * (-t136 * t151 + t155 * t222) - g(2) * (t136 * t155 + t151 * t222), -t11 * t62 + t172 * t19, t11 * t61 + t165 * t62 + t172 * t20 + t19 * t47, t137 * t62 - t138 * t19, -t165 * t61 + t20 * t47, -t137 * t61 - t138 * t20, 0, t134 * t176 + t137 * t25 + t138 * t4 - t165 * t75 + t20 * t70 + t29 * t61 + t47 * t51, -t11 * t75 - t133 * t176 - t137 * t26 - t138 * t3 - t172 * t51 - t19 * t70 + t29 * t62, -t1 * t61 + t11 * t25 + t13 * t19 - t14 * t20 + t165 * t26 + t172 * t4 - t2 * t62 - t3 * t47 - t177, t1 * t26 + t14 * t3 + t2 * t25 + t13 * t4 + t29 * t75 + t70 * t51 - g(1) * (-t103 * t151 - t144 * t155) - g(2) * (t103 * t155 - t144 * t151); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t190, t201 * t158, t193, t190, t192, qJDD(2), t166 * t150 - t216, g(3) * t150 + t166 * t154, 0, 0, t212, t41, t27, -t212, -t173, t142, -t143 * t63 + (t142 * t153 - t143 * t198 - t87 * t199) * pkin(2) + t160, t143 * t64 + (-t142 * t149 - t143 * t197 - t89 * t199) * pkin(2) + t168, (t55 + t63) * t89 + (-t54 + t64) * t87 + (-t149 * t40 + t153 * t39 + (t149 * t89 - t153 * t87) * qJD(3)) * pkin(2), -t54 * t63 - t55 * t64 + (-t216 - t149 * t181 + t153 * t22 + (-t149 * t54 + t153 * t55) * qJD(3) + (qJD(1) * t113 + t177) * t150) * pkin(2), -t225, t9, t5, t225, t6, t137, t137 * t83 + t210 * t138 - t47 * t73 + t161, -t137 * t84 - t211 * t138 + t172 * t73 + t163, t11 * t83 + t165 * t84 + t172 * t210 - t211 * t47 + t175, t1 * t84 + t2 * t83 - t70 * t73 - g(3) * t183 + t211 * t14 + t210 * t13 - t177 * (-pkin(2) * t150 - pkin(3) * t139); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t212, t41, t27, -t212, -t173, t142, t143 * t55 + t160, t143 * t54 + t168, 0, 0, -t225, t9, t5, t225, t6, t137, -t138 * t15 + (t137 * t152 - t138 * t196 - t47 * t89) * pkin(3) + t161, t138 * t16 + (-t137 * t148 - t138 * t195 + t172 * t89) * pkin(3) + t163, -t15 * t172 + t16 * t47 + (t11 * t152 + t165 * t148 + (-t148 * t172 - t152 * t47) * qJD(4)) * pkin(3) + t175, -t13 * t15 - t14 * t16 + (t1 * t148 + t152 * t2 - t70 * t89 + (-t13 * t148 + t14 * t152) * qJD(4) + t226) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t225, t9, t5, t225, t6, t137, t138 * t14 + t161, t13 * t138 + t163, 0, 0;];
tau_reg = t7;
