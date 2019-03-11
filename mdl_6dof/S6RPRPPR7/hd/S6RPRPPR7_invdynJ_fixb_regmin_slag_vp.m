% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPPR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:57:28
% EndTime: 2019-03-09 02:57:32
% DurationCPUTime: 2.21s
% Computational Cost: add. (2624->346), mult. (5130->421), div. (0->0), fcn. (3457->10), ass. (0->187)
t133 = sin(qJ(1));
t136 = cos(qJ(1));
t205 = g(1) * t133 - g(2) * t136;
t132 = sin(qJ(3));
t202 = qJD(1) * t132;
t84 = pkin(3) * t202 + qJD(1) * qJ(2) + qJD(4);
t253 = -qJD(1) * t84 - t205;
t139 = qJD(1) ^ 2;
t151 = -qJ(2) * t139 - t205;
t135 = cos(qJ(3));
t221 = cos(pkin(9));
t181 = t221 * t135;
t169 = qJD(1) * t181;
t220 = sin(pkin(9));
t180 = t220 * t132;
t94 = qJD(1) * t180;
t77 = t169 - t94;
t67 = qJD(6) + t77;
t252 = t67 ^ 2;
t72 = t77 ^ 2;
t148 = t221 * t132 + t220 * t135;
t73 = t148 * qJD(1);
t251 = -t73 ^ 2 - t72;
t137 = -pkin(1) - pkin(7);
t90 = t137 * qJD(1) + qJD(2);
t65 = -qJ(4) * t202 + t132 * t90;
t56 = t220 * t65;
t201 = qJD(1) * t135;
t66 = -qJ(4) * t201 + t135 * t90;
t40 = t221 * t66 - t56;
t208 = -qJD(5) + t40;
t125 = qJDD(1) * qJ(2);
t168 = g(1) * t136 + g(2) * t133;
t126 = qJD(1) * qJD(2);
t188 = 0.2e1 * t126;
t250 = 0.2e1 * t125 + t188 - t168;
t242 = t73 * pkin(5);
t186 = t221 * t65;
t59 = qJD(3) * pkin(3) + t66;
t32 = t220 * t59 + t186;
t27 = -qJD(3) * qJ(5) - t32;
t16 = -t27 - t242;
t39 = t220 * t66 + t186;
t174 = qJDD(1) * t220;
t175 = qJDD(1) * t221;
t155 = t132 * t174 - t135 * t175;
t195 = qJD(1) * qJD(3);
t48 = t148 * t195 + t155;
t45 = -qJDD(6) + t48;
t106 = -t221 * pkin(3) - pkin(4);
t98 = -pkin(8) + t106;
t249 = -t98 * t45 + (t16 - t39 + t242) * t67;
t131 = sin(qJ(6));
t134 = cos(qJ(6));
t190 = -qJD(3) * t169 - t132 * t175 - t135 * t174;
t47 = qJD(3) * t94 + t190;
t53 = qJD(3) * t134 + t131 * t73;
t15 = qJD(6) * t53 + t131 * qJDD(3) + t134 * t47;
t191 = t135 * qJDD(1);
t194 = qJD(1) * qJD(4);
t200 = qJD(3) * t132;
t89 = t137 * qJDD(1) + qJDD(2);
t82 = t135 * t89;
t30 = -t135 * t194 - t90 * t200 + qJDD(3) * pkin(3) + t82 + (t132 * t195 - t191) * qJ(4);
t199 = qJD(3) * t135;
t36 = (-qJ(4) * qJD(1) + t90) * t199 + (-qJ(4) * qJDD(1) - t194 + t89) * t132;
t10 = -t220 * t36 + t221 * t30;
t11 = t220 * t30 + t221 * t36;
t31 = t221 * t59 - t56;
t177 = qJD(3) * t220;
t178 = qJD(3) * t221;
t75 = t132 * t177 - t135 * t178;
t76 = -t132 * t178 - t135 * t177;
t80 = -t180 + t181;
t248 = -t10 * t80 - t11 * t148 - t31 * t76 + t32 * t75;
t157 = qJD(5) - t31;
t25 = -qJD(3) * pkin(4) + t157;
t189 = qJDD(3) * qJ(5) + t11;
t8 = -qJD(3) * qJD(5) - t189;
t156 = qJDD(5) - t10;
t215 = qJDD(3) * pkin(4);
t9 = t156 - t215;
t247 = t148 * t8 + t25 * t76 - t27 * t75 + t80 * t9;
t244 = pkin(4) + pkin(8);
t243 = t47 * pkin(4);
t241 = t77 * pkin(5);
t123 = qJ(3) + pkin(9);
t113 = sin(t123);
t238 = g(3) * t113;
t237 = g(3) * t132;
t119 = t132 * pkin(3);
t234 = t16 * t148;
t210 = qJ(2) + t119;
t170 = -t80 * qJ(5) + t210;
t26 = t148 * t244 + t170;
t232 = t26 * t45;
t51 = qJD(3) * t131 - t134 * t73;
t228 = t51 * t67;
t227 = t53 * t67;
t226 = t53 * t73;
t225 = t73 * t51;
t224 = t131 * t45;
t42 = t134 * t45;
t197 = qJD(6) * t134;
t198 = qJD(6) * t131;
t14 = -qJD(3) * t198 + t134 * qJDD(3) - t131 * t47 + t73 * t197;
t223 = t14 * t134;
t222 = t48 * qJ(5);
t219 = pkin(1) * qJDD(1);
t214 = t131 * t133;
t213 = t131 * t136;
t212 = t133 * t134;
t211 = t134 * t136;
t209 = qJ(4) - t137;
t207 = t241 - t208;
t206 = t136 * pkin(1) + t133 * qJ(2);
t129 = t135 ^ 2;
t204 = t132 ^ 2 - t129;
t138 = qJD(3) ^ 2;
t203 = -t138 - t139;
t196 = pkin(3) * t199 + qJD(2);
t193 = qJDD(3) * t132;
t192 = t132 * qJDD(1);
t187 = t148 * t197;
t185 = t135 * t195;
t184 = -t133 * pkin(1) + t136 * qJ(2);
t162 = qJDD(4) + t125 + t126 + (t185 + t192) * pkin(3);
t144 = -t77 * qJD(5) + t162 + t222;
t1 = -t244 * t47 + t144;
t13 = -t244 * qJD(3) + t157 + t241;
t183 = qJD(6) * t13 + t1;
t160 = -t77 * qJ(5) + t84;
t17 = t244 * t73 + t160;
t3 = -t48 * pkin(5) - t244 * qJDD(3) + t156;
t182 = -qJD(6) * t17 + t3;
t179 = pkin(3) * t201 + t73 * qJ(5);
t176 = t131 * t67;
t86 = t209 * t135;
t173 = qJD(6) * t80 + qJD(1);
t171 = qJDD(2) - t219;
t166 = t14 * t148 - t75 * t53;
t164 = -t148 * t45 - t67 * t75;
t62 = -t135 * qJD(4) + t209 * t200;
t63 = -qJD(3) * t86 - t132 * qJD(4);
t34 = t220 * t63 - t221 * t62;
t85 = t209 * t132;
t49 = -t220 * t85 + t221 * t86;
t114 = cos(t123);
t163 = pkin(4) * t113 - qJ(5) * t114;
t6 = t13 * t131 + t134 * t17;
t130 = -qJ(4) - pkin(7);
t161 = t136 * t119 + t133 * t130 + t184;
t159 = t133 * t119 - t136 * t130 + t206;
t158 = -t176 * t67 - t42;
t154 = 0.2e1 * qJ(2) * t195 + qJDD(3) * t137;
t37 = t80 * pkin(5) + t49;
t4 = pkin(5) * t47 - t8;
t153 = t148 * t4 - t16 * t75 + t37 * t45;
t150 = -t76 * qJ(5) - t80 * qJD(5) + t196;
t149 = -t252 * t134 + t224;
t35 = t220 * t62 + t221 * t63;
t50 = -t220 * t86 - t221 * t85;
t147 = -g(3) * t114 - t113 * t205;
t146 = t148 * t47 + t48 * t80 + t73 * t75 - t76 * t77;
t145 = t162 - t168;
t143 = t34 * t77 - t35 * t73 + t47 * t50 - t48 * t49 + t205;
t142 = -t137 * t138 + t250;
t33 = pkin(4) * t73 + t160;
t141 = t205 * t114 + t33 * t77 + t156 - t238;
t140 = t4 + (-qJD(6) * t98 + t244 * t77 + t179) * t67 + t147;
t116 = qJDD(3) * t135;
t102 = t220 * pkin(3) + qJ(5);
t71 = -t114 * t214 + t211;
t70 = -t114 * t212 - t213;
t69 = -t114 * t213 - t212;
t68 = -t114 * t211 + t214;
t46 = pkin(4) * t148 + t170;
t41 = pkin(4) * t77 + t179;
t38 = -pkin(5) * t148 + t50;
t23 = -pkin(4) * t75 + t150;
t19 = t75 * pkin(5) + t35;
t18 = t76 * pkin(5) + t34;
t12 = -t244 * t75 + t150;
t7 = t144 - t243;
t5 = t13 * t134 - t131 * t17;
t2 = t134 * t3;
t20 = [qJDD(1), t205, t168, qJDD(2) - t205 - 0.2e1 * t219, t250, -t171 * pkin(1) - g(1) * t184 - g(2) * t206 + (t188 + t125) * qJ(2), qJDD(1) * t129 - 0.2e1 * t132 * t185, -0.2e1 * t132 * t191 + 0.2e1 * t195 * t204, -t132 * t138 + t116, -t135 * t138 - t193, 0, t132 * t142 + t135 * t154, -t132 * t154 + t135 * t142, t143 + t248, -g(1) * t161 - g(2) * t159 - t10 * t49 + t11 * t50 + t162 * t210 + t84 * t196 - t31 * t34 + t32 * t35, t143 + t247, t34 * qJD(3) + t49 * qJDD(3) + t113 * t168 - t148 * t7 - t23 * t73 + t33 * t75 + t46 * t47, t35 * qJD(3) + t50 * qJDD(3) + t114 * t168 - t23 * t77 - t33 * t76 + t46 * t48 - t7 * t80, t7 * t46 + t33 * t23 - t8 * t50 - t27 * t35 + t9 * t49 + t25 * t34 - g(1) * (t136 * t163 + t161) - g(2) * (t133 * t163 + t159) t131 * t166 + t187 * t53 (t131 * t51 - t134 * t53) * t75 - (t131 * t15 - t223 + (t131 * t53 + t134 * t51) * qJD(6)) * t148, t131 * t164 + t14 * t80 + t187 * t67 + t53 * t76, -t148 * t198 * t67 + t134 * t164 - t15 * t80 - t51 * t76, -t45 * t80 + t67 * t76, -g(1) * t69 - g(2) * t71 + t38 * t15 + t19 * t51 + t2 * t80 + t5 * t76 + (-t1 * t80 - t12 * t67 + t232) * t131 + (t18 * t67 - t153) * t134 + ((-t131 * t37 - t134 * t26) * t67 - t6 * t80 + t131 * t234) * qJD(6), -g(1) * t68 - g(2) * t70 + t38 * t14 + t19 * t53 - t6 * t76 + (-(qJD(6) * t37 + t12) * t67 + t232 - t183 * t80 + qJD(6) * t234) * t134 + (-(-qJD(6) * t26 + t18) * t67 - t182 * t80 + t153) * t131; 0, 0, 0, qJDD(1), -t139, t171 + t151, 0, 0, 0, 0, 0, t132 * t203 + t116, t135 * t203 - t193, t146, -t248 + t253, t146, qJD(1) * t73 - qJD(3) * t76 - qJDD(3) * t80, qJD(1) * t77 - qJD(3) * t75 + qJDD(3) * t148, -qJD(1) * t33 - t205 - t247, 0, 0, 0, 0, 0, t80 * t42 + t148 * t15 - t75 * t51 + (t131 * t173 - t134 * t76) * t67, -t80 * t224 + (t131 * t76 + t134 * t173) * t67 + t166; 0, 0, 0, 0, 0, 0, t135 * t139 * t132, -t204 * t139, t191, -t192, qJDD(3), t135 * t151 + t237 + t82, g(3) * t135 + (-t151 - t89) * t132 (t32 - t39) * t77 + (-t31 + t40) * t73 + (t220 * t47 + t221 * t48) * pkin(3), t31 * t39 - t32 * t40 + (t221 * t10 + t220 * t11 + t253 * t135 + t237) * pkin(3), t102 * t47 - t106 * t48 + (-t27 - t39) * t77 + (t25 + t208) * t73, -t39 * qJD(3) + t41 * t73 + (-pkin(4) + t106) * qJDD(3) + t141, t102 * qJDD(3) - t33 * t73 + t41 * t77 + (0.2e1 * qJD(5) - t40) * qJD(3) + t147 + t189, -t8 * t102 + t9 * t106 - t33 * t41 - t25 * t39 - g(3) * (-t119 - t163) + t208 * t27 - t205 * (pkin(3) * t135 + pkin(4) * t114 + qJ(5) * t113) -t176 * t53 + t223 (-t15 - t227) * t134 + (-t14 + t228) * t131, t158 + t226, t149 - t225, t67 * t73, t102 * t15 + t140 * t131 + t249 * t134 + t207 * t51 + t5 * t73, t102 * t14 - t249 * t131 + t140 * t134 + t207 * t53 - t6 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, t31 * t77 + t32 * t73 + t145, t251 (t94 - t77) * qJD(3) + t190, qJD(3) * t73 + t48, -t243 + t222 - t27 * t73 + (-qJD(5) - t25) * t77 + t145, 0, 0, 0, 0, 0, t149 + t225, t131 * t252 + t226 + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, -t73 * t77 + qJDD(3), -t72 - t138, t27 * qJD(3) + t141 - t215, 0, 0, 0, 0, 0, -qJD(3) * t51 + t158, -qJD(3) * t53 + t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t51, -t51 ^ 2 + t53 ^ 2, t14 + t228, -t15 + t227, -t45, -g(1) * t70 + g(2) * t68 - t131 * t1 - t134 * t238 - t16 * t53 + t2 + (-qJD(6) + t67) * t6, g(1) * t71 - g(2) * t69 + t16 * t51 + t5 * t67 - t183 * t134 + (-t182 + t238) * t131;];
tau_reg  = t20;
