% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:55
% EndTime: 2022-01-20 11:03:01
% DurationCPUTime: 1.67s
% Computational Cost: add. (2429->249), mult. (3676->313), div. (0->0), fcn. (2827->16), ass. (0->166)
t145 = sin(qJ(5));
t149 = cos(qJ(5));
t141 = qJD(1) + qJD(2);
t144 = cos(pkin(9));
t150 = cos(qJ(4));
t209 = t150 * t144;
t193 = t141 * t209;
t143 = sin(pkin(9));
t146 = sin(qJ(4));
t212 = t146 * t143;
t194 = t141 * t212;
t70 = -t193 + t194;
t91 = t150 * t143 + t146 * t144;
t72 = t91 * t141;
t169 = t145 * t70 - t149 * t72;
t136 = qJDD(1) + qJDD(2);
t192 = qJD(4) * t193 + t91 * t136;
t37 = -qJD(4) * t194 + t192;
t90 = -t209 + t212;
t173 = t90 * t136;
t81 = t91 * qJD(4);
t38 = t141 * t81 + t173;
t156 = t169 * qJD(5) - t145 * t37 - t149 * t38;
t140 = qJD(4) + qJD(5);
t219 = t169 * t140;
t243 = t156 - t219;
t151 = cos(qJ(2));
t204 = qJD(1) * t151;
t195 = pkin(1) * t204;
t175 = qJD(3) - t195;
t30 = t145 * t72 + t149 * t70;
t218 = t30 * t140;
t199 = qJD(5) * t149;
t200 = qJD(5) * t145;
t7 = -t145 * t38 + t149 * t37 - t70 * t199 - t72 * t200;
t242 = t7 + t218;
t241 = t169 * t30;
t205 = t143 ^ 2 + t144 ^ 2;
t147 = sin(qJ(2));
t198 = qJDD(1) * t147;
t202 = qJD(2) * t151;
t59 = t136 * qJ(3) + t141 * qJD(3) + (qJD(1) * t202 + t198) * pkin(1);
t182 = t205 * t59;
t142 = qJ(1) + qJ(2);
t131 = sin(t142);
t124 = g(1) * t131;
t132 = cos(t142);
t227 = g(2) * t132;
t240 = t124 - t227;
t239 = t169 ^ 2 - t30 ^ 2;
t139 = pkin(9) + qJ(4);
t130 = qJ(5) + t139;
t118 = sin(t130);
t225 = t147 * pkin(1);
t196 = qJD(1) * t225;
t99 = t141 * qJ(3) + t196;
t189 = pkin(7) * t141 + t99;
t60 = t189 * t143;
t61 = t189 * t144;
t167 = t146 * t60 - t150 * t61;
t21 = -t70 * pkin(8) - t167;
t119 = cos(t130);
t213 = t119 * t132;
t214 = t119 * t131;
t121 = -t144 * pkin(3) - pkin(2);
t68 = t121 * t141 + t175;
t41 = t70 * pkin(4) + t68;
t190 = pkin(7) * t136 + t59;
t45 = t190 * t143;
t46 = t190 * t144;
t184 = -t146 * t46 - t150 * t45;
t5 = qJDD(4) * pkin(4) - t37 * pkin(8) + t167 * qJD(4) + t184;
t238 = t41 * t30 + g(3) * t118 + t21 * t200 + g(2) * t214 + g(1) * t213 + (-t21 * t140 - t5) * t145;
t236 = -t146 * t61 - t150 * t60;
t20 = -t72 * pkin(8) + t236;
t19 = qJD(4) * pkin(4) + t20;
t220 = t149 * t21;
t172 = -t145 * t19 - t220;
t215 = t118 * t132;
t216 = t118 * t131;
t168 = -t146 * t45 + t150 * t46;
t6 = -t38 * pkin(8) + t236 * qJD(4) + t168;
t237 = g(1) * t215 + g(2) * t216 - g(3) * t119 + t172 * qJD(5) - t145 * t6 + t149 * t5 + t41 * t169;
t104 = (-pkin(7) - qJ(3)) * t143;
t133 = t144 * pkin(7);
t105 = t144 * qJ(3) + t133;
t222 = t146 * t104 + t150 * t105;
t235 = -t222 * qJD(4) - t175 * t91;
t234 = g(1) * t132 + g(2) * t131;
t224 = t151 * pkin(1);
t206 = -qJD(2) * t196 + qJDD(1) * t224;
t233 = t124 + t206;
t201 = qJD(4) * t150;
t232 = (qJD(3) * t143 + qJD(4) * t105) * t146 - t90 * t195 - qJD(3) * t209 - t104 * t201;
t80 = t90 * qJD(4);
t231 = t80 * pkin(8);
t230 = t81 * pkin(4);
t229 = t90 * pkin(4);
t228 = t91 * pkin(8);
t226 = t136 * pkin(2);
t120 = qJ(3) + t225;
t82 = (-pkin(7) - t120) * t143;
t83 = t144 * t120 + t133;
t223 = t146 * t82 + t150 * t83;
t208 = t132 * pkin(2) + t131 * qJ(3);
t203 = qJD(2) * t147;
t197 = pkin(1) * t203;
t191 = t141 * t203;
t187 = qJDD(3) - t206;
t69 = t187 - t226;
t188 = -t69 - t227;
t186 = -t131 * pkin(2) + t132 * qJ(3);
t183 = -t146 * t83 + t150 * t82;
t181 = t150 * t104 - t146 * t105;
t111 = pkin(1) * t202 + qJD(3);
t180 = t111 * t205;
t179 = t205 * t136;
t178 = -t234 + t182;
t42 = t181 - t228;
t78 = t81 * pkin(8);
t177 = -qJD(5) * t42 + t232 + t78;
t85 = t90 * pkin(8);
t43 = -t85 + t222;
t176 = qJD(5) * t43 - t231 - t235;
t27 = t183 - t228;
t28 = -t85 + t223;
t171 = -t145 * t28 + t149 * t27;
t170 = t145 * t27 + t149 * t28;
t51 = t145 * t91 + t149 * t90;
t52 = -t145 * t90 + t149 * t91;
t103 = t121 - t224;
t165 = -t196 + t230;
t22 = -t51 * qJD(5) - t145 * t81 - t149 * t80;
t53 = t121 * t136 + t187;
t24 = t38 * pkin(4) + t53;
t164 = -g(1) * t216 + g(2) * t215 + t41 * t22 + t24 * t52;
t23 = t52 * qJD(5) - t145 * t80 + t149 * t81;
t163 = g(1) * t214 - g(2) * t213 + t41 * t23 + t24 * t51;
t128 = sin(t139);
t162 = -t240 * t128 + t53 * t91 - t68 * t80;
t129 = cos(t139);
t161 = t240 * t129 + t53 * t90 + t68 * t81;
t160 = t141 * t196 - t227;
t158 = t82 * t201 + t111 * t209 + (-qJD(4) * t83 - t111 * t143) * t146;
t157 = t175 * t205;
t154 = -t223 * qJD(4) - t91 * t111;
t152 = cos(qJ(1));
t148 = sin(qJ(1));
t135 = qJDD(4) + qJDD(5);
t126 = -pkin(2) - t224;
t110 = t144 * t124;
t92 = -t141 * pkin(2) + t175;
t67 = t121 + t229;
t62 = t197 + t230;
t57 = t103 + t229;
t50 = -t81 * qJD(4) - t90 * qJDD(4);
t49 = -t80 * qJD(4) + t91 * qJDD(4);
t17 = t154 + t231;
t16 = -t78 + t158;
t15 = t37 * t91 - t72 * t80;
t10 = -t51 * t135 - t23 * t140;
t9 = t52 * t135 + t22 * t140;
t4 = -t37 * t90 - t91 * t38 + t80 * t70 - t72 * t81;
t2 = -t169 * t22 + t7 * t52;
t1 = t156 * t52 + t169 * t23 - t22 * t30 - t7 * t51;
t3 = [qJDD(1), g(1) * t148 - g(2) * t152, g(1) * t152 + g(2) * t148, t136, -t227 + (t136 * t151 - t191) * pkin(1) + t233, ((-qJDD(1) - t136) * t147 + (-qJD(1) - t141) * t202) * pkin(1) + t234, t110 + (-pkin(1) * t191 - t126 * t136 + t188) * t144, t120 * t179 + t141 * t180 + t178, t69 * t126 + t92 * t197 - g(1) * (-t148 * pkin(1) + t186) - g(2) * (t152 * pkin(1) + t208) + t120 * t182 + t99 * t180, t15, t4, t49, t50, 0, t154 * qJD(4) + t183 * qJDD(4) + t103 * t38 + t70 * t197 + t161, -t158 * qJD(4) - t223 * qJDD(4) + t103 * t37 + t72 * t197 + t162, t2, t1, t9, t10, 0, t62 * t30 - t57 * t156 + (-t170 * qJD(5) - t145 * t16 + t149 * t17) * t140 + t171 * t135 + t163, -t62 * t169 + t57 * t7 - (t171 * qJD(5) + t145 * t17 + t149 * t16) * t140 - t170 * t135 + t164; 0, 0, 0, t136, t160 + t233, (-t198 + (-qJD(2) + t141) * t204) * pkin(1) + t234, t110 + (t160 - t69 + t226) * t144, qJ(3) * t179 + t157 * t141 + t178, -t69 * pkin(2) - g(1) * t186 - g(2) * t208 + qJ(3) * t182 + t157 * t99 - t92 * t196, t15, t4, t49, t50, 0, t235 * qJD(4) + t181 * qJDD(4) + t121 * t38 - t70 * t196 + t161, t232 * qJD(4) - t222 * qJDD(4) + t121 * t37 - t72 * t196 + t162, t2, t1, t9, t10, 0, -t67 * t156 + (-t145 * t43 + t149 * t42) * t135 + t165 * t30 + (t177 * t145 - t176 * t149) * t140 + t163, t67 * t7 - (t145 * t42 + t149 * t43) * t135 - t165 * t169 + (t176 * t145 + t177 * t149) * t140 + t164; 0, 0, 0, 0, 0, 0, -t144 * t136, -t205 * t141 ^ 2, -t205 * t99 * t141 - t124 - t188, 0, 0, 0, 0, 0, 0.2e1 * t72 * qJD(4) + t173, (-t70 - t194) * qJD(4) + t192, 0, 0, 0, 0, 0, -t156 - t219, t7 - t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, t72 * t70, -t70 ^ 2 + t72 ^ 2, (t70 - t194) * qJD(4) + t192, -t173, qJDD(4), -g(3) * t129 + t128 * t234 - t68 * t72 + t184, g(3) * t128 + t129 * t234 + t68 * t70 - t168, -t241, t239, t242, t243, t135, -(-t145 * t20 - t220) * t140 + (t149 * t135 - t140 * t200 - t72 * t30) * pkin(4) + t237, (-qJD(5) * t19 + t20 * t140 - t6) * t149 + (-t145 * t135 - t140 * t199 + t169 * t72) * pkin(4) + t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t241, t239, t242, t243, t135, -t172 * t140 + t237, (-t6 + (-qJD(5) + t140) * t19) * t149 + t238;];
tau_reg = t3;
