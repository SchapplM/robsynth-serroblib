% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:30:00
% EndTime: 2019-12-31 19:30:04
% DurationCPUTime: 1.52s
% Computational Cost: add. (1707->275), mult. (3944->350), div. (0->0), fcn. (2831->10), ass. (0->150)
t118 = qJD(2) - qJD(5);
t128 = sin(qJ(5));
t131 = cos(qJ(5));
t125 = sin(pkin(8));
t126 = cos(pkin(8));
t132 = cos(qJ(2));
t186 = t126 * t132;
t172 = qJD(1) * t186;
t129 = sin(qJ(2));
t182 = qJD(1) * t129;
t79 = t125 * t182 - t172;
t92 = t125 * t132 + t126 * t129;
t82 = t92 * qJD(1);
t211 = t128 * t79 + t131 * t82;
t198 = t118 * t211;
t178 = t132 * qJDD(1);
t179 = t129 * qJDD(1);
t159 = t125 * t179 - t126 * t178;
t81 = t92 * qJD(2);
t50 = qJD(1) * t81 + t159;
t180 = qJD(1) * qJD(2);
t171 = t129 * t180;
t143 = t92 * qJDD(1) - t125 * t171;
t170 = t132 * t180;
t51 = t126 * t170 + t143;
t5 = qJD(5) * t211 + t128 * t51 - t131 * t50;
t223 = t5 + t198;
t115 = t132 * pkin(2);
t216 = t115 + pkin(1);
t99 = -qJD(1) * t216 + qJD(3);
t222 = -qJ(4) * t82 + t99;
t145 = pkin(2) * t171 - qJDD(1) * t216 + qJDD(3);
t221 = -qJ(4) * t51 + t145;
t154 = t128 * t82 - t131 * t79;
t166 = t154 * qJD(5) - t128 * t50 - t131 * t51;
t197 = t118 * t154;
t220 = t166 + t197;
t219 = -t154 ^ 2 + t211 ^ 2;
t207 = -pkin(3) - pkin(4);
t14 = t207 * t79 - t222;
t119 = qJ(2) + pkin(8);
t113 = sin(t119);
t114 = cos(t119);
t150 = t113 * t128 + t114 * t131;
t127 = -qJ(3) - pkin(6);
t100 = t127 * t129;
t167 = qJD(2) * t127;
t75 = -qJD(3) * t129 + t132 * t167;
t46 = qJDD(2) * pkin(2) + t75 * qJD(1) + qJDD(1) * t100;
t101 = t127 * t132;
t74 = qJD(3) * t132 + t129 * t167;
t54 = t74 * qJD(1) - qJDD(1) * t101;
t9 = -t125 * t54 + t126 * t46;
t174 = -qJDD(4) + t9;
t2 = -pkin(7) * t51 + t207 * qJDD(2) - t174;
t121 = qJD(2) * qJD(4);
t10 = t125 * t46 + t126 * t54;
t176 = qJDD(2) * qJ(4) + t10;
t7 = t121 + t176;
t3 = pkin(7) * t50 + t7;
t130 = sin(qJ(1));
t212 = -t113 * t131 + t114 * t128;
t59 = t212 * t130;
t133 = cos(qJ(1));
t61 = t212 * t133;
t218 = -g(1) * t61 - g(2) * t59 - g(3) * t150 + t128 * t3 - t131 * t2 + t14 * t211;
t77 = t82 ^ 2;
t215 = -t79 ^ 2 - t77;
t214 = t211 * t154;
t97 = qJD(1) * t101;
t196 = t125 * t97;
t96 = qJD(1) * t100;
t56 = t126 * t96 + t196;
t184 = qJD(4) - t56;
t213 = g(1) * t130 - g(2) * t133;
t162 = g(1) * t133 + g(2) * t130;
t210 = qJD(5) + t118;
t60 = t150 * t130;
t62 = t150 * t133;
t209 = g(1) * t62 + g(2) * t60 - g(3) * t212 - t128 * t2 - t131 * t3 + t14 * t154;
t26 = pkin(3) * t79 + t222;
t208 = -g(3) * t114 + t162 * t113 - t26 * t82 + t174;
t206 = pkin(7) * t79;
t205 = pkin(7) * t82;
t201 = g(3) * t132;
t29 = t125 * t75 + t126 * t74;
t195 = t126 * t97;
t90 = qJD(2) * pkin(2) + t96;
t49 = t125 * t90 - t195;
t58 = t125 * t100 - t126 * t101;
t190 = qJD(4) * t82;
t189 = qJDD(2) * pkin(3);
t185 = -t205 + t184;
t123 = t129 ^ 2;
t183 = -t132 ^ 2 + t123;
t181 = qJD(2) * t129;
t40 = qJD(2) * qJ(4) + t49;
t175 = pkin(2) * t181;
t111 = -pkin(2) * t126 - pkin(3);
t28 = t125 * t74 - t126 * t75;
t55 = t125 * t96 - t195;
t48 = t126 * t90 + t196;
t57 = -t126 * t100 - t101 * t125;
t165 = t118 ^ 2;
t164 = qJ(4) * t92 + t216;
t163 = qJD(4) - t48;
t160 = pkin(3) * t114 + qJ(4) * t113;
t16 = t207 * qJD(2) + t163 - t205;
t22 = t40 + t206;
t158 = t128 * t22 - t131 * t16;
t157 = -t128 * t16 - t131 * t22;
t30 = -pkin(7) * t92 + t57;
t91 = t125 * t129 - t186;
t31 = pkin(7) * t91 + t58;
t156 = -t128 * t31 + t131 * t30;
t155 = t128 * t30 + t131 * t31;
t52 = t128 * t92 - t131 * t91;
t53 = t128 * t91 + t131 * t92;
t107 = -pkin(4) + t111;
t109 = pkin(2) * t125 + qJ(4);
t152 = t107 * t131 - t109 * t128;
t151 = t107 * t128 + t109 * t131;
t149 = -pkin(2) * t182 - qJ(4) * t79;
t148 = -0.2e1 * pkin(1) * t180 - pkin(6) * qJDD(2);
t84 = qJD(2) * t186 - t125 * t181;
t144 = qJ(4) * t84 + qJD(4) * t92 - t175;
t134 = qJD(2) ^ 2;
t142 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t134 + t213;
t135 = qJD(1) ^ 2;
t141 = pkin(1) * t135 - pkin(6) * qJDD(1) + t162;
t139 = t28 * t82 - t29 * t79 - t58 * t50 + t51 * t57 - t162;
t136 = pkin(3) * t50 + t221;
t117 = qJDD(2) - qJDD(5);
t103 = t133 * t216;
t47 = pkin(3) * t91 - t164;
t32 = -qJD(2) * pkin(3) + t163;
t27 = pkin(3) * t82 - t149;
t24 = t55 + t206;
t23 = t207 * t91 + t164;
t21 = pkin(3) * t81 - t144;
t19 = pkin(7) * t81 + t29;
t18 = -pkin(7) * t84 + t28;
t17 = t207 * t82 + t149;
t13 = t207 * t81 + t144;
t12 = qJD(5) * t53 + t128 * t84 - t131 * t81;
t11 = -qJD(5) * t52 + t128 * t81 + t131 * t84;
t8 = -t174 - t189;
t6 = t136 - t190;
t1 = t207 * t50 + t190 - t221;
t4 = [qJDD(1), t213, t162, qJDD(1) * t123 + 0.2e1 * t129 * t170, 0.2e1 * t129 * t178 - 0.2e1 * t183 * t180, qJDD(2) * t129 + t132 * t134, qJDD(2) * t132 - t129 * t134, 0, t129 * t148 + t132 * t142, -t129 * t142 + t132 * t148, -t10 * t91 - t48 * t84 - t49 * t81 - t9 * t92 + t139, t10 * t58 + t49 * t29 - t9 * t57 - t48 * t28 - t145 * t216 + t99 * t175 - g(1) * (-t127 * t133 - t130 * t216) - g(2) * (-t127 * t130 + t103), -qJD(2) * t28 - qJDD(2) * t57 + t114 * t213 + t21 * t79 + t26 * t81 + t47 * t50 + t6 * t91, t32 * t84 - t40 * t81 - t7 * t91 + t8 * t92 + t139, qJD(2) * t29 + qJDD(2) * t58 + t113 * t213 - t21 * t82 - t26 * t84 - t47 * t51 - t6 * t92, -g(2) * t103 + t26 * t21 + t32 * t28 + t40 * t29 + t6 * t47 + t8 * t57 + t7 * t58 + (g(1) * t127 - g(2) * t160) * t133 + (-g(1) * (-t216 - t160) + g(2) * t127) * t130, t11 * t211 - t166 * t53, -t11 * t154 - t12 * t211 + t166 * t52 - t5 * t53, -t11 * t118 - t117 * t53, t117 * t52 + t118 * t12, 0, t13 * t154 + t23 * t5 + t1 * t52 + t14 * t12 - (-qJD(5) * t155 - t128 * t19 + t131 * t18) * t118 - t156 * t117 + g(1) * t60 - g(2) * t62, t13 * t211 - t23 * t166 + t1 * t53 + t14 * t11 + (qJD(5) * t156 + t128 * t18 + t131 * t19) * t118 + t155 * t117 - g(1) * t59 + g(2) * t61; 0, 0, 0, -t129 * t135 * t132, t183 * t135, t179, t178, qJDD(2), t129 * t141 - t201, g(3) * t129 + t132 * t141, (t49 - t55) * t82 + (-t48 + t56) * t79 + (-t125 * t50 - t126 * t51) * pkin(2), t48 * t55 - t49 * t56 + (-t201 + t10 * t125 + t126 * t9 + (-qJD(1) * t99 + t162) * t129) * pkin(2), qJD(2) * t55 - t27 * t79 + (pkin(3) - t111) * qJDD(2) + t208, -t109 * t50 + t111 * t51 + (t40 - t55) * t82 + (t32 - t184) * t79, -g(3) * t113 - qJD(2) * t56 + qJDD(2) * t109 - t162 * t114 - t26 * t79 + t27 * t82 + 0.2e1 * t121 + t176, t7 * t109 + t8 * t111 - t26 * t27 - t32 * t55 - g(3) * (t115 + t160) + t184 * t40 + t162 * (pkin(2) * t129 + pkin(3) * t113 - qJ(4) * t114), -t214, -t219, t220, t223, t117, -t152 * t117 - t17 * t154 + (t128 * t185 + t131 * t24) * t118 + (t118 * t151 - t157) * qJD(5) + t218, t151 * t117 - t17 * t211 + (-t128 * t24 + t131 * t185) * t118 + (t118 * t152 - t158) * qJD(5) - t209; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t215, t48 * t82 + t49 * t79 + t145 - t213, 0.2e1 * qJD(2) * t82 + t159, t215, (t79 - t172) * qJD(2) - t143, t40 * t79 + (-qJD(4) - t32) * t82 + t136 - t213, 0, 0, 0, 0, 0, -t5 + t198, t166 - t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79 * t82 - qJDD(2), (t79 + t172) * qJD(2) + t143, -t77 - t134, -qJD(2) * t40 - t189 - t208, 0, 0, 0, 0, 0, -t131 * t117 - t128 * t165 - t154 * t82, t128 * t117 - t131 * t165 - t211 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t214, t219, -t220, -t223, -t117, t210 * t157 - t218, t210 * t158 + t209;];
tau_reg = t4;
