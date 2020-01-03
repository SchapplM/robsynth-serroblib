% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:51:19
% EndTime: 2019-12-31 19:51:22
% DurationCPUTime: 1.76s
% Computational Cost: add. (3102->284), mult. (4635->315), div. (0->0), fcn. (3104->12), ass. (0->170)
t148 = cos(qJ(2));
t212 = qJD(1) * t148;
t177 = -pkin(1) * t212 + qJD(3);
t143 = cos(pkin(8));
t237 = cos(qJ(4));
t198 = t237 * t143;
t142 = sin(pkin(8));
t145 = sin(qJ(4));
t219 = t145 * t142;
t165 = t198 - t219;
t144 = -pkin(7) - qJ(3);
t104 = t144 * t142;
t133 = t143 * pkin(7);
t105 = qJ(3) * t143 + t133;
t166 = t237 * t104 - t145 * t105;
t231 = qJD(4) * t166 + t177 * t165;
t137 = t142 ^ 2;
t138 = t143 ^ 2;
t213 = t137 + t138;
t136 = qJDD(1) + qJDD(2);
t140 = qJD(1) + qJD(2);
t146 = sin(qJ(2));
t208 = qJDD(1) * t146;
t210 = qJD(2) * t148;
t68 = qJ(3) * t136 + qJD(3) * t140 + (qJD(1) * t210 + t208) * pkin(1);
t188 = t213 * t68;
t141 = qJ(1) + qJ(2);
t132 = cos(t141);
t124 = g(2) * t132;
t236 = pkin(1) * t146;
t202 = qJD(1) * t236;
t234 = pkin(1) * t148;
t214 = -qJD(2) * t202 + qJDD(1) * t234;
t191 = qJDD(3) - t214;
t233 = pkin(2) * t136;
t78 = t191 - t233;
t252 = t78 + t124;
t90 = t237 * t142 + t145 * t143;
t209 = qJD(4) * t145;
t196 = t142 * t209;
t192 = qJD(4) * t237;
t180 = t143 * t192;
t201 = t90 * t136 + t140 * t180;
t41 = t140 * t196 - t201;
t174 = t165 * t136;
t85 = t90 * qJD(4);
t42 = t140 * t85 - t174;
t121 = pkin(3) * t143 + pkin(2);
t60 = -t121 * t136 + t191;
t152 = pkin(4) * t42 + qJ(5) * t41 + t60;
t247 = t90 * t140;
t11 = -qJD(5) * t247 + t152;
t77 = -t121 * t140 + t177;
t200 = t140 * t219;
t79 = -t140 * t198 + t200;
t31 = pkin(4) * t79 - qJ(5) * t247 + t77;
t84 = -t180 + t196;
t251 = -t11 * t90 + t31 * t84;
t250 = -t11 * t165 + t31 * t85;
t249 = -t165 * t60 + t77 * t85;
t248 = t60 * t90 - t77 * t84;
t238 = t247 ^ 2;
t74 = t79 ^ 2;
t246 = -t74 - t238;
t245 = -t74 + t238;
t131 = sin(t141);
t244 = g(1) * t132 + g(2) * t131;
t125 = g(1) * t131;
t243 = t124 - t125;
t195 = pkin(7) * t136 + t68;
t50 = t195 * t142;
t51 = t195 * t143;
t92 = qJ(3) * t140 + t202;
t194 = pkin(7) * t140 + t92;
t69 = t194 * t142;
t206 = -t145 * t50 - t69 * t192 + t237 * t51;
t70 = t194 * t143;
t12 = -t209 * t70 + t206;
t187 = t145 * t51 + t70 * t192 - t69 * t209 + t237 * t50;
t229 = t145 * t70;
t35 = -t237 * t69 - t229;
t36 = -t145 * t69 + t237 * t70;
t242 = t12 * t165 + t187 * t90 + t35 * t84 - t36 * t85;
t225 = qJDD(4) * pkin(4);
t10 = qJDD(5) + t187 - t225;
t217 = qJD(5) - t35;
t32 = -qJD(4) * pkin(4) + t217;
t33 = qJD(4) * qJ(5) + t36;
t207 = qJDD(4) * qJ(5);
t9 = t207 + (qJD(5) - t229) * qJD(4) + t206;
t241 = t10 * t90 + t165 * t9 - t32 * t84 - t33 * t85;
t139 = pkin(8) + qJ(4);
t129 = sin(t139);
t130 = cos(t139);
t240 = -g(3) * t129 - t244 * t130 + t206;
t223 = t129 * t132;
t224 = t129 * t131;
t227 = -g(1) * t224 + g(2) * t223;
t67 = t145 * t104 + t237 * t105;
t239 = -t231 * qJD(4) - qJDD(4) * t67 + t227;
t147 = sin(qJ(1));
t235 = pkin(1) * t147;
t232 = t247 * t79;
t230 = qJD(4) * t67 + t177 * t90;
t228 = t252 * t142;
t226 = qJD(4) * t36;
t222 = t130 * t132;
t221 = t132 * t144;
t220 = t143 * t136;
t216 = t132 * pkin(2) + t131 * qJ(3);
t211 = qJD(2) * t146;
t94 = t132 * t121;
t205 = pkin(4) * t222 + qJ(5) * t223 + t94;
t203 = pkin(1) * t211;
t197 = t140 * t211;
t190 = -pkin(2) * t131 + t132 * qJ(3);
t189 = t35 + t229;
t186 = -t131 * t144 + t94;
t112 = pkin(1) * t210 + qJD(3);
t185 = t112 * t213;
t184 = t213 * t136;
t183 = -t244 + t188;
t182 = t140 * t202;
t181 = -t214 + t243;
t37 = pkin(4) * t85 + qJ(5) * t84 - qJD(5) * t90;
t179 = -t37 + t202;
t178 = -g(2) * t222 + t130 * t125;
t149 = cos(qJ(1));
t175 = g(1) * t147 - g(2) * t149;
t173 = -t165 * t42 + t79 * t85;
t171 = pkin(4) * t130 + qJ(5) * t129;
t168 = -t121 * t131 - t221;
t59 = qJD(4) * t85 - qJDD(4) * t165;
t120 = qJ(3) + t236;
t86 = (-pkin(7) - t120) * t142;
t87 = t120 * t143 + t133;
t167 = -t145 * t87 + t237 * t86;
t55 = t145 * t86 + t237 * t87;
t164 = -t182 - t233;
t127 = -pkin(2) - t234;
t163 = pkin(1) * t197 + t127 * t136;
t26 = qJD(4) * t167 + t112 * t165;
t162 = qJD(4) * t26 + qJDD(4) * t55 - t227;
t56 = -pkin(4) * t165 - qJ(5) * t90 - t121;
t161 = g(1) * t223 + g(2) * t224 - g(3) * t130 - t187;
t27 = qJD(4) * t55 + t112 * t90;
t159 = t167 * t41 + t247 * t27 - t26 * t79 - t55 * t42 - t244;
t2 = -t165 * t41 - t247 * t85 - t42 * t90 + t79 * t84;
t158 = -qJD(4) * t27 + qJDD(4) * t167 + t178;
t157 = t177 * t213;
t156 = -t230 * qJD(4) + qJDD(4) * t166 + t178;
t155 = t247 * t31 + qJDD(5) - t161;
t154 = (-g(1) * (-t121 - t171) + g(2) * t144) * t131;
t153 = t41 * t166 + t230 * t247 - t231 * t79 - t67 * t42 - t244;
t151 = 0.2e1 * t247 * qJD(4) - t174;
t134 = t149 * pkin(1);
t115 = t138 * t136;
t114 = t137 * t136;
t111 = t143 * t125;
t103 = -t121 - t234;
t95 = 0.2e1 * t142 * t220;
t91 = -pkin(2) * t140 + t177;
t57 = -qJD(4) * t84 + qJDD(4) * t90;
t49 = t56 - t234;
t40 = pkin(4) * t247 + qJ(5) * t79;
t34 = t37 + t203;
t30 = (t79 - t200) * qJD(4) + t201;
t29 = (t79 + t200) * qJD(4) - t201;
t14 = -t247 * t84 - t41 * t90;
t1 = [0, 0, 0, 0, 0, qJDD(1), t175, g(1) * t149 + g(2) * t147, 0, 0, 0, 0, 0, 0, 0, t136, (t136 * t148 - t197) * pkin(1) - t181, ((-qJDD(1) - t136) * t146 + (-qJD(1) - t140) * t210) * pkin(1) + t244, 0, (t175 + (t146 ^ 2 + t148 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t114, t95, 0, t115, 0, 0, t111 + (-t163 - t252) * t143, (t163 - t125) * t142 + t228, t120 * t184 + t140 * t185 + t183, t78 * t127 + t91 * t203 - g(1) * (t190 - t235) - g(2) * (t134 + t216) + t120 * t188 + t92 * t185, t14, t2, t57, t173, -t59, 0, t103 * t42 + t203 * t79 + t158 + t249, -t103 * t41 + t203 * t247 - t162 + t248, t159 + t242, t12 * t55 + t36 * t26 - t187 * t167 - t35 * t27 + t60 * t103 + t77 * t203 - g(1) * (t168 - t235) - g(2) * (t134 + t186), t14, t57, -t2, 0, t59, t173, t34 * t79 + t42 * t49 + t158 + t250, t159 + t241, -t247 * t34 + t41 * t49 + t162 + t251, t9 * t55 + t33 * t26 + t11 * t49 + t31 * t34 - t10 * t167 + t32 * t27 - g(1) * (-t221 - t235) - g(2) * (t134 + t205) + t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, -t181 + t182, (-t208 + (-qJD(2) + t140) * t212) * pkin(1) + t244, 0, 0, t114, t95, 0, t115, 0, 0, t111 + (-t164 - t252) * t143, (t164 - t125) * t142 + t228, qJ(3) * t184 + t140 * t157 + t183, -t78 * pkin(2) - g(1) * t190 - g(2) * t216 + qJ(3) * t188 + t157 * t92 - t202 * t91, t14, t2, t57, t173, -t59, 0, -t121 * t42 - t202 * t79 + t156 + t249, t121 * t41 - t202 * t247 + t239 + t248, t153 + t242, -g(1) * t168 - g(2) * t186 + t12 * t67 - t60 * t121 - t166 * t187 - t77 * t202 - t230 * t35 + t231 * t36, t14, t57, -t2, 0, t59, t173, -t179 * t79 + t42 * t56 + t156 + t250, t153 + t241, t179 * t247 + t41 * t56 - t239 + t251, g(1) * t221 - g(2) * t205 - t10 * t166 + t11 * t56 - t179 * t31 + t230 * t32 + t231 * t33 + t9 * t67 + t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t220, t142 * t136, -t213 * t140 ^ 2, -t140 * t213 * t92 + t243 + t78, 0, 0, 0, 0, 0, 0, t151, -t29, t246, t247 * t35 + t36 * t79 + t243 + t60, 0, 0, 0, 0, 0, 0, t151, t246, t29, t33 * t79 + (-qJD(5) - t32) * t247 + t152 + t243; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t232, t245, t30, -t232, t174, qJDD(4), -t247 * t77 + t161 + t226, qJD(4) * t189 + t77 * t79 - t240, 0, 0, t232, t30, -t245, qJDD(4), -t174, -t232, -t40 * t79 - t155 + 0.2e1 * t225 + t226, pkin(4) * t41 - qJ(5) * t42 + (t33 - t36) * t247 + (t32 - t217) * t79, 0.2e1 * t207 - t31 * t79 + t40 * t247 + (0.2e1 * qJD(5) - t189) * qJD(4) + t240, -t10 * pkin(4) - g(3) * t171 + t9 * qJ(5) + t217 * t33 - t31 * t40 - t32 * t36 + t244 * (pkin(4) * t129 - qJ(5) * t130); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) + t232, t30, -qJD(4) ^ 2 - t238, -qJD(4) * t33 + t155 - t225;];
tau_reg = t1;
