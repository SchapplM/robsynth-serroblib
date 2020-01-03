% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRR12
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR12_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR12_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR12_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:24
% EndTime: 2019-12-31 18:07:27
% DurationCPUTime: 1.23s
% Computational Cost: add. (1332->220), mult. (2725->288), div. (0->0), fcn. (1992->10), ass. (0->126)
t84 = sin(pkin(8));
t85 = cos(pkin(8));
t88 = sin(qJ(4));
t91 = cos(qJ(4));
t48 = t84 * t91 + t85 * t88;
t40 = t48 * qJD(1);
t165 = qJD(5) + t40;
t90 = cos(qJ(5));
t128 = t90 * qJD(4);
t133 = qJD(1) * t84;
t121 = t88 * t133;
t139 = t91 * t85;
t123 = qJD(1) * t139;
t42 = -t121 + t123;
t87 = sin(qJ(5));
t27 = t42 * t87 - t128;
t167 = t165 * t27;
t86 = -pkin(1) - qJ(3);
t156 = -qJD(1) * qJD(3) + qJDD(1) * t86;
t49 = qJDD(2) + t156;
t112 = -pkin(6) * qJDD(1) + t49;
t30 = t112 * t84;
t31 = t112 * t85;
t105 = t30 * t88 - t31 * t91;
t55 = t86 * qJD(1) + qJD(2);
t118 = -pkin(6) * qJD(1) + t55;
t32 = t118 * t84;
t33 = t118 * t85;
t15 = t32 * t91 + t33 * t88;
t2 = -qJDD(4) * pkin(4) + qJD(4) * t15 + t105;
t89 = sin(qJ(1));
t92 = cos(qJ(1));
t163 = g(1) * t89 - g(2) * t92;
t80 = pkin(8) + qJ(4);
t67 = sin(t80);
t68 = cos(t80);
t96 = g(3) * t67 - t163 * t68;
t166 = -(pkin(4) * t42 + t165 * pkin(7)) * t165 - t2 + t96;
t113 = t90 * t165;
t131 = qJD(4) * t91;
t122 = t85 * t131;
t97 = -qJD(4) * t121 + t48 * qJDD(1);
t24 = qJD(1) * t122 + t97;
t20 = qJDD(5) + t24;
t147 = t20 * t87;
t164 = -t113 * t165 - t147;
t135 = t84 ^ 2 + t85 ^ 2;
t81 = qJDD(1) * qJ(2);
t82 = qJD(1) * qJD(2);
t162 = t81 + t82;
t160 = t135 * t55;
t109 = g(1) * t92 + g(2) * t89;
t53 = qJDD(3) + t162;
t159 = t53 - t109;
t158 = qJD(4) * t40;
t104 = t32 * t88 - t33 * t91;
t11 = -qJD(4) * pkin(4) + t104;
t106 = t30 * t91 + t31 * t88;
t66 = qJD(1) * qJ(2) + qJD(3);
t52 = pkin(3) * t133 + t66;
t13 = pkin(4) * t40 - pkin(7) * t42 + t52;
t116 = qJDD(4) * pkin(7) - qJD(4) * t104 + qJD(5) * t13 + t106;
t47 = t84 * t88 - t139;
t61 = t84 * pkin(3) + qJ(2);
t21 = pkin(4) * t48 + pkin(7) * t47 + t61;
t149 = -pkin(6) + t86;
t50 = t149 * t84;
t51 = t149 * t85;
t26 = t50 * t91 + t51 * t88;
t132 = qJD(4) * t88;
t43 = -t131 * t84 - t132 * t85;
t25 = t50 * t88 - t51 * t91;
t9 = -qJD(3) * t48 - qJD(4) * t25;
t155 = t11 * t43 - (qJD(5) * t21 + t9) * t165 - t116 * t48 - t2 * t47 - t26 * t20;
t154 = 0.2e1 * t82;
t152 = g(3) * t68;
t130 = qJD(5) * t87;
t125 = t85 * qJDD(1);
t126 = t84 * qJDD(1);
t107 = t91 * t125 - t126 * t88;
t23 = t107 - t158;
t7 = qJD(5) * t128 + t87 * qJDD(4) - t130 * t42 + t90 * t23;
t151 = t47 * t7;
t150 = t7 * t87;
t148 = t11 * t47;
t146 = t21 * t20;
t145 = t27 * t42;
t29 = qJD(4) * t87 + t42 * t90;
t144 = t29 * t42;
t143 = t87 * t89;
t142 = t87 * t92;
t141 = t89 * t90;
t16 = t90 * t20;
t140 = t90 * t92;
t138 = t43 * qJD(4) - t47 * qJDD(4);
t137 = t92 * pkin(1) + t89 * qJ(2);
t134 = pkin(1) * qJDD(1);
t129 = qJD(5) * t90;
t120 = t135 * t49;
t119 = -t90 * qJDD(4) + t23 * t87;
t12 = qJD(4) * pkin(7) + t15;
t46 = pkin(3) * t126 + t53;
t6 = pkin(4) * t24 - pkin(7) * t23 + t46;
t115 = qJD(5) * t12 - t6;
t111 = qJD(5) * t48 + qJD(1);
t110 = qJDD(2) - t134;
t103 = t16 + (-t40 * t87 - t130) * t165;
t44 = -t132 * t84 + t122;
t102 = -qJD(4) * t44 - qJDD(4) * t48;
t101 = -t116 + t152;
t100 = t130 * t47 + t43 * t90;
t95 = -pkin(7) * t20 + (-t104 + t11) * t165;
t94 = t159 + t162;
t93 = qJD(1) ^ 2;
t74 = t92 * qJ(2);
t39 = t67 * t140 - t143;
t38 = t67 * t142 + t141;
t37 = t67 * t141 + t142;
t36 = -t67 * t143 + t140;
t18 = pkin(4) * t44 - pkin(7) * t43 + qJD(2);
t10 = -qJD(3) * t47 + qJD(4) * t26;
t8 = qJD(5) * t29 + t119;
t5 = t90 * t6;
t4 = t12 * t90 + t13 * t87;
t3 = -t12 * t87 + t13 * t90;
t1 = [qJDD(1), t163, t109, qJDD(2) - 0.2e1 * t134 - t163, -t109 + 0.2e1 * t81 + t154, -t110 * pkin(1) - g(1) * (-pkin(1) * t89 + t74) - g(2) * t137 + (t81 + t154) * qJ(2), t94 * t84, t94 * t85, t163 + t135 * (-t156 - t49), t53 * qJ(2) + t66 * qJD(2) - g(1) * (t86 * t89 + t74) - g(2) * (qJ(3) * t92 + t137) + t86 * t120 - qJD(3) * t160, -t23 * t47 + t42 * t43, -t23 * t48 + t24 * t47 - t40 * t43 - t42 * t44, t138, t102, 0, qJD(2) * t40 - qJD(4) * t10 - qJDD(4) * t25 - t109 * t67 + t24 * t61 + t44 * t52 + t46 * t48, qJD(2) * t42 - qJD(4) * t9 - qJDD(4) * t26 - t109 * t68 + t23 * t61 + t43 * t52 - t46 * t47, t100 * t29 - t90 * t151, (-t27 * t90 - t29 * t87) * t43 + (t150 + t8 * t90 + (-t27 * t87 + t29 * t90) * qJD(5)) * t47, t100 * t165 - t47 * t16 + t29 * t44 + t48 * t7, t47 * t147 - t27 * t44 - t48 * t8 + (t129 * t47 - t43 * t87) * t165, t165 * t44 + t20 * t48, -g(1) * t39 - g(2) * t37 + t10 * t27 + t25 * t8 + t3 * t44 + t5 * t48 + (t18 * t165 + t146 + (-t12 * t48 - t165 * t26 - t148) * qJD(5)) * t90 + t155 * t87, g(1) * t38 - g(2) * t36 + t10 * t29 + t25 * t7 - t4 * t44 + (-(-qJD(5) * t26 + t18) * t165 - t146 + t115 * t48 + qJD(5) * t148) * t87 + t155 * t90; 0, 0, 0, qJDD(1), -t93, -qJ(2) * t93 + t110 - t163, -t93 * t84, -t93 * t85, -t135 * qJDD(1), -qJD(1) * t66 + t120 - t163, 0, 0, 0, 0, 0, -qJD(1) * t40 + t138, -qJD(1) * t42 + t102, 0, 0, 0, 0, 0, -t48 * t147 - t27 * t43 + t47 * t8 + (-t111 * t90 - t44 * t87) * t165, -t48 * t16 - t29 * t43 + t151 + (t111 * t87 - t44 * t90) * t165; 0, 0, 0, 0, 0, 0, t126, t125, -t135 * t93, qJD(1) * t160 + t159, 0, 0, 0, 0, 0, (t42 + t123) * qJD(4) + t97, t107 - 0.2e1 * t158, 0, 0, 0, 0, 0, t103 - t145, -t144 + t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t40, -t40 ^ 2 + t42 ^ 2, t107, (t42 - t123) * qJD(4) - t97, qJDD(4), -t42 * t52 - t105 + t96, t163 * t67 + t40 * t52 - t106 + t152, t113 * t29 + t150, (t7 - t167) * t90 + (-t165 * t29 - t8) * t87, -t144 - t164, t103 + t145, -t165 * t42, -pkin(4) * t8 - t15 * t27 + t166 * t90 - t3 * t42 + t95 * t87, -pkin(4) * t7 - t15 * t29 - t166 * t87 + t4 * t42 + t95 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 * t27, -t27 ^ 2 + t29 ^ 2, t7 + t167, -t119 + (-qJD(5) + t165) * t29, t20, -g(1) * t36 - g(2) * t38 + t101 * t87 - t11 * t29 - t12 * t129 + t165 * t4 + t5, g(1) * t37 - g(2) * t39 + t101 * t90 + t11 * t27 + t115 * t87 + t165 * t3;];
tau_reg = t1;
