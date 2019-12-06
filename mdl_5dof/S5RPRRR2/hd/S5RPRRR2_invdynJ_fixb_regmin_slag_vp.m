% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:12:19
% EndTime: 2019-12-05 18:12:27
% DurationCPUTime: 2.48s
% Computational Cost: add. (3500->267), mult. (8699->349), div. (0->0), fcn. (7152->16), ass. (0->162)
t143 = sin(qJ(5));
t147 = cos(qJ(5));
t144 = sin(qJ(4));
t148 = cos(qJ(4));
t142 = cos(pkin(9));
t149 = cos(qJ(3));
t204 = t149 * t142;
t141 = sin(pkin(9));
t145 = sin(qJ(3));
t208 = t141 * t145;
t99 = -t204 + t208;
t91 = t99 * qJD(1);
t100 = t141 * t149 + t142 * t145;
t92 = t100 * qJD(1);
t170 = -t144 * t91 + t148 * t92;
t68 = t144 * t92 + t148 * t91;
t172 = t143 * t68 - t147 * t170;
t201 = qJD(4) * t148;
t202 = qJD(4) * t144;
t192 = qJD(1) * t208;
t195 = t142 * qJDD(1);
t196 = t141 * qJDD(1);
t198 = t149 * qJD(3);
t193 = qJD(1) * t142 * t198 + t145 * t195 + t149 * t196;
t76 = -qJD(3) * t192 + t193;
t176 = t145 * t196 - t149 * t195;
t94 = t100 * qJD(3);
t77 = qJD(1) * t94 + t176;
t25 = -t144 * t77 + t148 * t76 - t201 * t91 - t202 * t92;
t26 = qJD(4) * t170 + t144 * t76 + t148 * t77;
t159 = qJD(5) * t172 - t143 * t25 - t147 * t26;
t194 = -qJD(4) - qJD(5);
t134 = qJD(3) - t194;
t215 = t134 * t172;
t232 = t159 - t215;
t30 = t143 * t170 + t147 * t68;
t214 = t134 * t30;
t199 = qJD(5) * t147;
t200 = qJD(5) * t143;
t4 = -t143 * t26 + t147 * t25 - t170 * t200 - t199 * t68;
t236 = t4 + t214;
t240 = t172 * t30;
t237 = t172 ^ 2 - t30 ^ 2;
t139 = pkin(9) + qJ(3);
t135 = qJ(4) + t139;
t128 = qJ(5) + t135;
t122 = sin(t128);
t123 = cos(t128);
t219 = pkin(6) + qJ(2);
t113 = t219 * t141;
t101 = qJD(1) * t113;
t114 = t219 * t142;
t102 = qJD(1) * t114;
t168 = t101 * t145 - t102 * t149;
t55 = -pkin(7) * t91 - t168;
t52 = t148 * t55;
t229 = -t101 * t149 - t102 * t145;
t54 = -pkin(7) * t92 + t229;
t53 = qJD(3) * pkin(3) + t54;
t171 = -t144 * t53 - t52;
t242 = pkin(8) * t68;
t14 = -t171 - t242;
t146 = sin(qJ(1));
t150 = cos(qJ(1));
t178 = g(1) * t150 + g(2) * t146;
t127 = -pkin(2) * t142 - pkin(1);
t108 = qJD(1) * t127 + qJD(2);
t80 = pkin(3) * t91 + t108;
t40 = pkin(4) * t68 + t80;
t245 = g(3) * t122 + t123 * t178 + t14 * t200 + t40 * t30;
t136 = qJDD(3) + qJDD(4);
t197 = qJD(1) * qJD(2);
t223 = qJDD(1) * t219 + t197;
t84 = t223 * t141;
t85 = t223 * t142;
t184 = -t145 * t85 - t149 * t84;
t20 = qJDD(3) * pkin(3) - pkin(7) * t76 + qJD(3) * t168 + t184;
t169 = -t145 * t84 + t149 * t85;
t24 = -pkin(7) * t77 + qJD(3) * t229 + t169;
t161 = qJD(4) * t171 - t144 * t24 + t148 * t20;
t2 = pkin(4) * t136 - pkin(8) * t25 + t161;
t224 = (qJD(4) * t53 + t24) * t148 + t144 * t20 - t55 * t202;
t3 = -pkin(8) * t26 + t224;
t244 = -g(3) * t123 + t122 * t178 - t143 * t3 + t147 * t2 + t40 * t172;
t140 = qJD(3) + qJD(4);
t211 = t140 * t68;
t243 = t25 + t211;
t241 = t170 * t68;
t212 = t140 * t170;
t239 = -t26 + t212;
t238 = t170 ^ 2 - t68 ^ 2;
t125 = sin(t135);
t126 = cos(t135);
t235 = g(3) * t125 + t126 * t178 + t80 * t68 - t224;
t234 = (-t134 * t14 - t2) * t143 + t245;
t50 = t144 * t55;
t188 = t148 * t53 - t50;
t230 = pkin(8) * t170;
t13 = t188 - t230;
t11 = pkin(4) * t140 + t13;
t213 = t14 * t147;
t175 = -t11 * t143 - t213;
t233 = t175 * qJD(5) + t244;
t226 = qJ(2) * qJDD(1);
t225 = -g(3) * t126 + t125 * t178 - t170 * t80 + t161;
t222 = pkin(3) * t94;
t218 = t148 * t54 - t50;
t182 = -t113 * t149 - t114 * t145;
t65 = -pkin(7) * t100 + t182;
t216 = -t113 * t145 + t114 * t149;
t66 = -pkin(7) * t99 + t216;
t217 = t144 * t65 + t148 * t66;
t210 = qJDD(1) * pkin(1);
t131 = qJDD(5) + t136;
t207 = t144 * t131;
t206 = t144 * t147;
t205 = t147 * t131;
t203 = t141 ^ 2 + t142 ^ 2;
t191 = -qJD(5) * t11 - t3;
t187 = -t144 * t54 - t52;
t186 = -t144 * t66 + t148 * t65;
t181 = t203 * qJD(1) ^ 2;
t179 = 0.2e1 * t203;
t177 = g(1) * t146 - g(2) * t150;
t79 = t100 * t148 - t144 * t99;
t17 = -pkin(8) * t79 + t186;
t78 = t100 * t144 + t148 * t99;
t18 = -pkin(8) * t78 + t217;
t174 = -t143 * t18 + t147 * t17;
t173 = t143 * t17 + t147 * t18;
t38 = t143 * t79 + t147 * t78;
t39 = -t143 * t78 + t147 * t79;
t82 = pkin(3) * t99 + t127;
t162 = qJD(2) * t204 - t113 * t198 + (-qJD(2) * t141 - qJD(3) * t114) * t145;
t44 = -pkin(7) * t94 + t162;
t153 = -qJD(2) * t100 - qJD(3) * t216;
t93 = t99 * qJD(3);
t45 = pkin(7) * t93 + t153;
t167 = t144 * t45 + t148 * t44 + t201 * t65 - t202 * t66;
t103 = qJDD(1) * t127 + qJDD(2);
t166 = -t177 - t210;
t130 = qJDD(2) - t210;
t164 = -t130 - t166;
t56 = pkin(3) * t77 + t103;
t160 = -qJD(4) * t217 - t144 * t44 + t148 * t45;
t157 = t179 * t197 - t178;
t133 = cos(t139);
t132 = sin(t139);
t129 = pkin(3) * t148 + pkin(4);
t48 = pkin(4) * t78 + t82;
t46 = pkin(3) * t92 + pkin(4) * t170;
t37 = qJD(4) * t79 - t144 * t93 + t148 * t94;
t36 = -qJD(4) * t78 - t144 * t94 - t148 * t93;
t27 = pkin(4) * t37 + t222;
t16 = t218 - t230;
t15 = t187 + t242;
t10 = pkin(4) * t26 + t56;
t9 = qJD(5) * t39 + t143 * t36 + t147 * t37;
t8 = -qJD(5) * t38 - t143 * t37 + t147 * t36;
t7 = -pkin(8) * t36 + t160;
t6 = -pkin(8) * t37 + t167;
t1 = [qJDD(1), t177, t178, t164 * t142, -t164 * t141, t179 * t226 + t157, (-t130 + t177) * pkin(1) + (t203 * t226 + t157) * qJ(2), t100 * t76 - t92 * t93, -t100 * t77 - t76 * t99 + t91 * t93 - t92 * t94, -qJD(3) * t93 + qJDD(3) * t100, -qJD(3) * t94 - qJDD(3) * t99, 0, qJD(3) * t153 + qJDD(3) * t182 + t103 * t99 + t108 * t94 + t127 * t77 + t133 * t177, -qJD(3) * t162 - qJDD(3) * t216 + t103 * t100 - t108 * t93 + t127 * t76 - t132 * t177, t170 * t36 + t25 * t79, -t170 * t37 - t25 * t78 - t26 * t79 - t36 * t68, t136 * t79 + t140 * t36, -t136 * t78 - t140 * t37, 0, t126 * t177 + t136 * t186 + t140 * t160 + t222 * t68 + t82 * t26 + t80 * t37 + t56 * t78, -t125 * t177 - t136 * t217 - t140 * t167 + t170 * t222 + t82 * t25 + t80 * t36 + t56 * t79, -t172 * t8 + t39 * t4, t159 * t39 + t172 * t9 - t30 * t8 - t38 * t4, t131 * t39 + t134 * t8, -t131 * t38 - t134 * t9, 0, t27 * t30 - t48 * t159 + t10 * t38 + t40 * t9 + (-qJD(5) * t173 - t143 * t6 + t147 * t7) * t134 + t174 * t131 + t177 * t123, -t27 * t172 + t48 * t4 + t10 * t39 + t40 * t8 - (qJD(5) * t174 + t143 * t7 + t147 * t6) * t134 - t173 * t131 - t177 * t122; 0, 0, 0, -t195, t196, -t181, -qJ(2) * t181 + qJDD(2) + t166, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t92 + t176, (-t91 - t192) * qJD(3) + t193, 0, 0, 0, 0, 0, t26 + t212, t25 - t211, 0, 0, 0, 0, 0, -t159 - t215, t4 - t214; 0, 0, 0, 0, 0, 0, 0, t92 * t91, -t91 ^ 2 + t92 ^ 2, (t91 - t192) * qJD(3) + t193, -t176, qJDD(3), -g(3) * t133 - t108 * t92 + t132 * t178 + t184, g(3) * t132 + t108 * t91 + t133 * t178 - t169, t241, t238, t243, t239, t136, -t187 * t140 + (t136 * t148 - t140 * t202 - t68 * t92) * pkin(3) + t225, t218 * t140 + (-t136 * t144 - t140 * t201 - t170 * t92) * pkin(3) + t235, -t240, t237, t236, t232, t131, t129 * t205 - t46 * t30 - (-t143 * t16 + t147 * t15) * t134 + (-t143 * t207 + (-t143 * t148 - t206) * t134 * qJD(4)) * pkin(3) + ((-pkin(3) * t206 - t129 * t143) * t134 + t175) * qJD(5) + t244, t46 * t172 + (-t129 * t131 - t2 + (-pkin(3) * t144 * t194 + t15) * t134) * t143 + (-pkin(3) * t207 + (-pkin(3) * t201 - qJD(5) * t129 + t16) * t134 + t191) * t147 + t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241, t238, t243, t239, t136, -t140 * t171 + t225, t140 * t188 + t235, -t240, t237, t236, t232, t131, -(-t13 * t143 - t213) * t134 + (-t134 * t200 - t170 * t30 + t205) * pkin(4) + t233, (t13 * t134 + t191) * t147 + (-t131 * t143 - t134 * t199 + t170 * t172) * pkin(4) + t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t240, t237, t236, t232, t131, -t134 * t175 + t233, (-t3 + (-qJD(5) + t134) * t11) * t147 + t234;];
tau_reg = t1;
