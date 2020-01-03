% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:04:09
% EndTime: 2019-12-31 19:04:14
% DurationCPUTime: 1.71s
% Computational Cost: add. (1609->253), mult. (3904->371), div. (0->0), fcn. (2644->8), ass. (0->144)
t108 = sin(qJ(5));
t111 = cos(qJ(5));
t109 = sin(qJ(4));
t110 = sin(qJ(3));
t166 = qJD(1) * t110;
t145 = t109 * t166;
t112 = cos(qJ(4));
t158 = t112 * qJD(3);
t76 = t145 - t158;
t159 = t109 * qJD(3);
t78 = t112 * t166 + t159;
t129 = t108 * t76 - t111 * t78;
t27 = t108 * t78 + t111 * t76;
t202 = t129 * t27;
t201 = t129 ^ 2 - t27 ^ 2;
t160 = qJD(5) * t111;
t161 = qJD(5) * t108;
t113 = cos(qJ(3));
t149 = t113 * t158;
t154 = qJD(3) * qJD(4);
t40 = qJD(1) * t149 - qJD(4) * t145 + t112 * t154;
t146 = t113 * t159;
t162 = qJD(4) * t112;
t147 = t110 * t162;
t119 = t146 + t147;
t41 = t119 * qJD(1) + t109 * t154;
t6 = -t108 * t41 + t111 * t40 - t76 * t160 - t78 * t161;
t157 = t113 * qJD(1);
t94 = -qJD(4) + t157;
t93 = -qJD(5) + t94;
t200 = -t27 * t93 + t6;
t96 = -cos(pkin(9)) * pkin(1) - pkin(2);
t71 = -t113 * pkin(3) - t110 * pkin(7) + t96;
t49 = t71 * qJD(1);
t182 = t109 * t49;
t102 = t110 * qJD(2);
t95 = sin(pkin(9)) * pkin(1) + pkin(6);
t88 = t95 * qJD(1);
t54 = t113 * t88 + t102;
t46 = qJD(3) * pkin(7) + t54;
t19 = t112 * t46 + t182;
t14 = -t76 * pkin(8) + t19;
t12 = t14 * t161;
t196 = t113 * qJD(2) - t110 * t88;
t45 = -qJD(3) * pkin(3) - t196;
t25 = t76 * pkin(4) + t45;
t199 = t25 * t27 + t12;
t116 = t129 * qJD(5) - t108 * t40 - t111 * t41;
t198 = t129 * t93 + t116;
t47 = t196 * qJD(3);
t131 = pkin(3) * t110 - pkin(7) * t113;
t85 = t131 * qJD(3);
t70 = qJD(1) * t85;
t139 = t109 * t47 - t112 * t70;
t117 = -t19 * qJD(4) - t139;
t155 = qJD(1) * qJD(3);
t97 = t110 * t155;
t4 = pkin(4) * t97 - t40 * pkin(8) + t117;
t153 = t109 * t70 + t112 * t47 + t49 * t162;
t163 = qJD(4) * t109;
t123 = -t46 * t163 + t153;
t5 = -t41 * pkin(8) + t123;
t152 = -t108 * t5 + t111 * t4;
t18 = -t109 * t46 + t112 * t49;
t13 = -t78 * pkin(8) + t18;
t10 = -t94 * pkin(4) + t13;
t180 = t111 * t14;
t2 = t108 * t10 + t180;
t197 = -t2 * qJD(5) + t25 * t129 + t152;
t80 = t108 * t112 + t111 * t109;
t57 = t80 * t110;
t195 = qJD(4) + qJD(5);
t104 = t110 ^ 2;
t127 = qJD(1) * t104 - t113 * t94;
t148 = t110 * t163;
t194 = -t127 * t158 - t94 * t148;
t193 = pkin(7) + pkin(8);
t192 = t76 * t94;
t191 = t78 * t94;
t171 = t110 * t112;
t173 = t109 * t110;
t16 = -t161 * t173 + (t195 * t171 + t146) * t111 + (t149 - t148) * t108;
t190 = t16 * t93 - t57 * t97;
t79 = t108 * t109 - t111 * t112;
t121 = t79 * t113;
t189 = qJD(1) * t121 - t195 * t79;
t188 = (-t157 + t195) * t80;
t82 = t131 * qJD(1);
t187 = t109 * t82 + t112 * t196;
t186 = t109 * t85 + t71 * t162;
t185 = t110 * t95 * t159 + t112 * t85;
t164 = qJD(3) * t113;
t48 = qJD(3) * t102 + t88 * t164;
t170 = t112 * t113;
t81 = t95 * t170;
t184 = t109 * t71 + t81;
t183 = t109 * t45;
t181 = t110 * t76;
t179 = t112 * t45;
t178 = t112 * t94;
t177 = t113 * t41;
t176 = t40 * t109;
t175 = t48 * t109;
t174 = t48 * t112;
t89 = qJD(1) * t96;
t172 = t109 * t113;
t114 = qJD(3) ^ 2;
t169 = t114 * t110;
t168 = t114 * t113;
t167 = -t113 ^ 2 + t104;
t165 = qJD(3) * t110;
t151 = qJD(4) * t193;
t150 = t109 * t157;
t144 = -t6 * t113 - t129 * t165;
t143 = t94 * t95 + t46;
t141 = qJD(5) * t10 + t5;
t138 = -t109 * t196 + t112 * t82;
t137 = -t40 * t113 + t78 * t165;
t135 = t94 * t147;
t134 = -t54 + (-t150 + t163) * pkin(4);
t125 = pkin(4) * t110 - pkin(8) * t170;
t91 = t193 * t112;
t133 = t125 * qJD(1) + qJD(5) * t91 + t112 * t151 + t138;
t90 = t193 * t109;
t132 = pkin(8) * t150 - qJD(5) * t90 - t109 * t151 - t187;
t60 = t112 * t71;
t22 = -pkin(8) * t171 + t60 + (-t109 * t95 - pkin(4)) * t113;
t24 = -pkin(8) * t173 + t184;
t130 = t108 * t22 + t111 * t24;
t128 = 0.2e1 * qJD(3) * t89;
t124 = -t113 * t116 - t27 * t165;
t122 = t127 * t109;
t15 = -qJD(3) * t121 - t195 * t57;
t58 = t79 * t110;
t120 = t15 * t93 + t58 * t97;
t115 = qJD(1) ^ 2;
t101 = -t112 * pkin(4) - pkin(3);
t63 = (pkin(4) * t109 + t95) * t110;
t33 = t119 * pkin(4) + t95 * t164;
t20 = t41 * pkin(4) + t48;
t9 = (-t110 * t158 - t113 * t163) * t95 - t119 * pkin(8) + t186;
t8 = t125 * qJD(3) + (-t81 + (pkin(8) * t110 - t71) * t109) * qJD(4) + t185;
t1 = t111 * t10 - t108 * t14;
t3 = [0, 0, 0, 0, 0.2e1 * t113 * t97, -0.2e1 * t167 * t155, t168, -t169, 0, t110 * t128 - t95 * t168, t113 * t128 + t95 * t169, t78 * t149 + (t40 * t112 - t78 * t163) * t110, (-t109 * t78 - t112 * t76) * t164 + (-t176 - t112 * t41 + (t109 * t76 - t112 * t78) * qJD(4)) * t110, t137 - t194, t135 + t177 + (-t122 - t181) * qJD(3), (-t94 - t157) * t165, -(-t163 * t71 + t185) * t94 + ((t76 * t95 + t183) * qJD(3) + (t112 * t143 + t182) * qJD(4) + t139) * t113 + (t45 * t162 + t175 + t95 * t41 + ((-t95 * t172 + t60) * qJD(1) + t18) * qJD(3)) * t110, t186 * t94 + (-t143 * t163 + (t78 * t95 + t179) * qJD(3) + t153) * t113 + (-t45 * t163 + t174 + t95 * t40 + (-t184 * qJD(1) - t95 * t178 - t19) * qJD(3)) * t110, -t129 * t15 - t6 * t58, -t116 * t58 + t129 * t16 - t15 * t27 - t6 * t57, -t120 + t144, t124 + t190, (-t93 - t157) * t165, -(-t108 * t9 + t111 * t8) * t93 - t152 * t113 + t33 * t27 - t63 * t116 + t20 * t57 + t25 * t16 + (t113 * t2 + t130 * t93) * qJD(5) + ((-t108 * t24 + t111 * t22) * qJD(1) + t1) * t165, -t12 * t113 + t25 * t15 - t20 * t58 - t33 * t129 + t63 * t6 + ((-qJD(5) * t24 + t8) * t93 + t4 * t113) * t108 + ((qJD(5) * t22 + t9) * t93 + t141 * t113) * t111 + (-t130 * qJD(1) - t2) * t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t169, -t168, 0, 0, 0, 0, 0, t135 - t177 + (-t122 + t181) * qJD(3), t137 + t194, 0, 0, 0, 0, 0, -t124 + t190, t120 + t144; 0, 0, 0, 0, -t110 * t115 * t113, t167 * t115, 0, 0, 0, t54 * qJD(3) - t89 * t166 - t48, -t89 * t157, -t78 * t178 + t176, (t40 + t192) * t112 + (-t41 + t191) * t109, -t94 * t162 + (t94 * t170 + (-t78 + t159) * t110) * qJD(1), t94 * t163 + (-t94 * t172 + (t76 + t158) * t110) * qJD(1), t94 * t166, -pkin(3) * t41 - t174 + t138 * t94 - t54 * t76 + (pkin(7) * t178 + t183) * qJD(4) + (-t18 * t110 + (-pkin(7) * t165 - t113 * t45) * t109) * qJD(1), -pkin(3) * t40 + t175 - t187 * t94 - t54 * t78 + (-t109 * pkin(7) * t94 + t179) * qJD(4) + (-t45 * t170 + (-pkin(7) * t158 + t19) * t110) * qJD(1), -t129 * t189 + t6 * t80, t116 * t80 + t129 * t188 - t189 * t27 - t6 * t79, -t189 * t93 + (qJD(3) * t80 + t129) * t166, t188 * t93 + (-qJD(3) * t79 + t27) * t166, t93 * t166, -t101 * t116 + t20 * t79 + (t108 * t132 + t111 * t133) * t93 + t134 * t27 + t188 * t25 + ((-t108 * t91 - t111 * t90) * qJD(3) - t1) * t166, t101 * t6 + t20 * t80 + (-t108 * t133 + t111 * t132) * t93 - t134 * t129 + t189 * t25 + (-(-t108 * t90 + t111 * t91) * qJD(3) + t2) * t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78 * t76, -t76 ^ 2 + t78 ^ 2, t40 - t192, -t191 - t41, t97, -t19 * t94 - t45 * t78 + t117, -t18 * t94 + t45 * t76 - t123, -t202, t201, t200, t198, t97, (-t108 * t13 - t180) * t93 + (t111 * t97 + t161 * t93 - t78 * t27) * pkin(4) + t197, (t14 * t93 - t4) * t108 + (-t13 * t93 - t141) * t111 + (-t108 * t97 + t129 * t78 + t160 * t93) * pkin(4) + t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202, t201, t200, t198, t97, -t2 * t93 + t197, -t1 * t93 - t108 * t4 - t111 * t141 + t199;];
tauc_reg = t3;
