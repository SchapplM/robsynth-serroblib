% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:58:33
% EndTime: 2019-12-05 18:58:37
% DurationCPUTime: 1.53s
% Computational Cost: add. (2163->241), mult. (2960->319), div. (0->0), fcn. (1907->16), ass. (0->172)
t126 = sin(qJ(3));
t131 = cos(qJ(3));
t132 = cos(qJ(2));
t190 = qJD(1) * t132;
t170 = qJD(2) * t190;
t127 = sin(qJ(2));
t184 = qJDD(1) * t127;
t142 = (t170 + t184) * pkin(1);
t218 = pkin(1) * t132;
t105 = qJDD(1) * t218;
t117 = qJDD(1) + qJDD(2);
t208 = pkin(1) * qJD(1);
t181 = t127 * t208;
t50 = pkin(2) * t117 - qJD(2) * t181 + t105;
t119 = qJD(1) + qJD(2);
t180 = pkin(1) * t190;
t69 = pkin(2) * t119 + t180;
t227 = -t126 * t50 - (qJD(3) * t69 + t142) * t131;
t123 = qJ(1) + qJ(2);
t114 = qJ(3) + t123;
t100 = cos(t114);
t99 = sin(t114);
t226 = g(2) * t100 + g(3) * t99;
t109 = qJD(3) + t119;
t124 = sin(qJ(5));
t130 = cos(qJ(4));
t125 = sin(qJ(4));
t129 = cos(qJ(5));
t196 = t125 * t129;
t59 = t124 * t130 + t196;
t46 = t59 * t109;
t188 = qJD(3) * t131;
t88 = t126 * t181;
t54 = t131 * t180 - t88;
t225 = pkin(2) * t188 - t54;
t223 = -g(2) * t99 + g(3) * t100;
t189 = qJD(3) * t126;
t194 = t127 * t131;
t153 = t126 * t132 + t194;
t53 = t153 * t208;
t161 = pkin(2) * t189 - t53;
t118 = qJD(4) + qJD(5);
t101 = pkin(2) * t126 + pkin(8);
t217 = pkin(2) * t131;
t102 = -pkin(3) - t217;
t108 = qJDD(3) + t117;
t134 = qJD(4) ^ 2;
t222 = -t101 * t134 - t102 * t108 - t161 * t109;
t221 = -pkin(8) - pkin(9);
t103 = pkin(2) + t218;
t209 = pkin(1) * t194 + t126 * t103;
t52 = pkin(8) + t209;
t219 = -pkin(9) - t52;
t216 = pkin(3) * t108;
t215 = pkin(3) * t109;
t214 = pkin(4) * t130;
t193 = t129 * t130;
t176 = t109 * t193;
t199 = t124 * t125;
t177 = t109 * t199;
t44 = -t176 + t177;
t213 = t46 * t44;
t212 = -pkin(9) - t101;
t163 = qJD(3) * t181;
t195 = t126 * t127;
t92 = pkin(1) * t195;
t162 = t126 * pkin(1) * t170 + qJDD(1) * t92 + t69 * t189 + (t163 - t50) * t131;
t13 = t162 - t216;
t186 = qJD(4) * t130;
t42 = t131 * t69 - t88;
t36 = -t42 - t215;
t211 = t13 * t125 + t36 * t186;
t187 = qJD(4) * t125;
t106 = pkin(4) * t187;
t210 = t106 + t161;
t43 = t126 * t69 + t131 * t181;
t207 = t109 * t43;
t122 = qJ(4) + qJ(5);
t112 = cos(t122);
t206 = t112 * t99;
t37 = pkin(8) * t109 + t43;
t173 = pkin(9) * t109 + t37;
t26 = t173 * t130;
t204 = t129 * t26;
t30 = t103 * t189 + (t153 * qJD(2) + t127 * t188) * pkin(1);
t203 = t30 * t109;
t201 = t100 * t112;
t200 = t109 * t125;
t197 = t125 * t108;
t192 = t130 * t108;
t120 = t125 ^ 2;
t191 = -t130 ^ 2 + t120;
t185 = qJD(5) * t124;
t183 = t226 * t130 + t36 * t187;
t111 = sin(t123);
t113 = cos(t123);
t182 = g(2) * t113 + g(3) * t111 + t105;
t104 = -pkin(3) - t214;
t175 = qJD(4) * t221;
t174 = t109 * t186;
t172 = -g(2) * t111 + g(3) * t113;
t171 = qJD(4) * t219;
t147 = t109 * t187 - t192;
t10 = pkin(4) * t147 + t13;
t28 = t104 * t109 - t42;
t35 = t118 * t59;
t58 = -t193 + t199;
t169 = g(2) * t201 + g(3) * t206 + t10 * t58 + t28 * t35;
t168 = qJD(4) * t212;
t167 = t103 * t131 - t92;
t166 = qJD(1) * (-qJD(2) + t119);
t165 = qJD(2) * (-qJD(1) - t119);
t82 = t126 * t163;
t164 = t82 + t223;
t51 = -pkin(3) - t167;
t160 = -t43 + t106;
t158 = t124 * t197 - t129 * t192;
t25 = t173 * t125;
t24 = qJD(4) * pkin(4) - t25;
t157 = -t124 * t24 - t204;
t38 = t219 * t125;
t115 = t130 * pkin(9);
t39 = t130 * t52 + t115;
t156 = -t124 * t39 + t129 * t38;
t155 = t124 * t38 + t129 * t39;
t57 = t101 * t130 + t115;
t152 = qJD(5) * t57 + t225 * t125 - t130 * t168;
t85 = pkin(8) * t130 + t115;
t151 = qJD(5) * t85 - t125 * t42 - t130 * t175;
t56 = t212 * t125;
t150 = -qJD(5) * t56 - t125 * t168 - t225 * t130;
t84 = t221 * t125;
t149 = -qJD(5) * t84 - t125 * t175 + t130 * t42;
t148 = -t162 + t226;
t145 = -pkin(8) * t134 + t207 + t216;
t144 = -t108 * t51 - t134 * t52 - t203;
t12 = pkin(8) * t108 - t227 - t82;
t143 = -t109 * t36 - t12 + t223;
t141 = -pkin(8) * qJDD(4) + (t42 - t215) * qJD(4);
t110 = sin(t122);
t34 = t118 * t58;
t140 = t10 * t59 - t110 * t226 - t28 * t34;
t29 = t103 * t188 + (-t127 * t189 + (t131 * t132 - t195) * qJD(2)) * pkin(1);
t139 = -qJDD(4) * t52 + (t109 * t51 - t29) * qJD(4);
t14 = qJD(5) * t176 + t108 * t196 - t118 * t177 + t124 * t192 + t129 * t174;
t138 = -qJDD(4) * t101 + (t102 * t109 - t225) * qJD(4);
t4 = -t37 * t186 + qJDD(4) * pkin(4) - t12 * t125 + (-t174 - t197) * pkin(9);
t137 = -g(2) * t206 + t28 * t44 + t26 * t185 + g(3) * t201 + g(1) * t110 + (-t26 * t118 - t4) * t124;
t5 = -pkin(9) * t147 + t12 * t130 - t37 * t187;
t136 = -g(1) * t112 + t157 * qJD(5) + t223 * t110 - t124 * t5 + t129 * t4 - t28 * t46;
t135 = t164 + t227;
t133 = cos(qJ(1));
t128 = sin(qJ(1));
t116 = qJDD(4) + qJDD(5);
t107 = t109 ^ 2;
t78 = t104 - t217;
t77 = qJDD(4) * t130 - t125 * t134;
t76 = qJDD(4) * t125 + t130 * t134;
t49 = t108 * t120 + 0.2e1 * t125 * t174;
t48 = t51 - t214;
t33 = -0.2e1 * t191 * t109 * qJD(4) + 0.2e1 * t125 * t192;
t27 = t106 + t30;
t22 = -t116 * t58 - t118 * t35;
t21 = t116 * t59 - t118 * t34;
t20 = -t125 * t29 + t130 * t171;
t19 = t125 * t171 + t130 * t29;
t18 = -t44 ^ 2 + t46 ^ 2;
t15 = t35 * t109 + t158;
t8 = t118 * t44 + t14;
t2 = t14 * t59 - t34 * t46;
t1 = -t14 * t58 - t15 * t59 + t34 * t44 - t35 * t46;
t3 = [qJDD(1), g(2) * t133 + g(3) * t128, -g(2) * t128 + g(3) * t133, t117, (t117 * t132 + t127 * t165) * pkin(1) + t182, ((-qJDD(1) - t117) * t127 + t132 * t165) * pkin(1) + t172, t108, t108 * t167 + t148 - t203, -t209 * t108 - t29 * t109 + t135, t49, t33, t76, t77, 0, t139 * t125 + (-t13 + t144) * t130 + t183, t139 * t130 + (-t144 - t226) * t125 + t211, t2, t1, t21, t22, 0, t27 * t44 + t48 * t15 + (-qJD(5) * t155 - t124 * t19 + t129 * t20) * t118 + t156 * t116 + t169, t27 * t46 + t48 * t14 - (qJD(5) * t156 + t124 * t20 + t129 * t19) * t118 - t155 * t116 + t140; 0, 0, 0, t117, pkin(1) * t127 * t166 + t182, (t132 * t166 - t184) * pkin(1) + t172, t108, t109 * t53 + (t108 * t131 - t109 * t189) * pkin(2) + t148, t109 * t54 + (-pkin(2) * t108 - t50) * t126 + ((-pkin(2) * t109 - t69) * qJD(3) - t142) * t131 + t164, t49, t33, t76, t77, 0, t138 * t125 + (-t13 + t222) * t130 + t183, t138 * t130 + (-t226 - t222) * t125 + t211, t2, t1, t21, t22, 0, t78 * t15 + (-t124 * t57 + t129 * t56) * t116 + t210 * t44 + (t124 * t150 - t129 * t152) * t118 + t169, t78 * t14 - (t124 * t56 + t129 * t57) * t116 + t210 * t46 + (t124 * t152 + t129 * t150) * t118 + t140; 0, 0, 0, 0, 0, 0, t108, t148 + t207, t109 * t42 + t135, t49, t33, t76, t77, 0, t141 * t125 + (-t13 + t145) * t130 + t183, t141 * t130 + (-t145 - t226) * t125 + t211, t2, t1, t21, t22, 0, t104 * t15 + (-t124 * t85 + t129 * t84) * t116 + t160 * t44 + (t124 * t149 - t129 * t151) * t118 + t169, t104 * t14 - (t124 * t84 + t129 * t85) * t116 + t160 * t46 + (t124 * t151 + t129 * t149) * t118 + t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125 * t107 * t130, t191 * t107, t197, t192, qJDD(4), -g(1) * t130 + t125 * t143, g(1) * t125 + t130 * t143, t213, t18, t8, -t158, t116, -(t124 * t25 - t204) * t118 + (t129 * t116 - t118 * t185 - t44 * t200) * pkin(4) + t136, (-qJD(5) * t24 - t25 * t118 - t5) * t129 + (-qJD(5) * t129 * t118 - t124 * t116 - t46 * t200) * pkin(4) + t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, t18, t8, -t158, t116, -t118 * t157 + t136, (-t5 + (-qJD(5) + t118) * t24) * t129 + t137;];
tau_reg = t3;
