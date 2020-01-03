% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR10_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:18
% EndTime: 2019-12-31 18:04:22
% DurationCPUTime: 1.24s
% Computational Cost: add. (1932->210), mult. (5014->300), div. (0->0), fcn. (3705->6), ass. (0->119)
t110 = sin(qJ(5));
t112 = cos(qJ(5));
t107 = sin(pkin(8));
t108 = cos(pkin(8));
t111 = sin(qJ(4));
t113 = cos(qJ(4));
t78 = t107 * t111 + t108 * t113;
t146 = qJD(1) * t78;
t144 = qJD(1) * t107;
t130 = t113 * t144;
t143 = qJD(1) * t108;
t131 = t111 * t143;
t71 = t130 - t131;
t159 = t110 * t71 + t112 * t146;
t167 = t159 ^ 2;
t64 = -qJD(1) * pkin(1) - pkin(2) * t143 - qJ(3) * t144 + qJD(2);
t49 = pkin(3) * t143 - t64;
t30 = pkin(4) * t146 + t49;
t166 = t30 * t159;
t106 = qJD(4) + qJD(5);
t165 = t106 * t159;
t118 = -t110 * t146 + t112 * t71;
t154 = t159 * t118;
t162 = t118 ^ 2;
t164 = t162 - t167;
t137 = qJD(5) * t112;
t138 = qJD(5) * t110;
t67 = t78 * qJD(4);
t60 = qJD(1) * t67;
t89 = qJD(4) * t131;
t61 = qJD(4) * t130 - t89;
t11 = t110 * t61 + t112 * t60 + t137 * t146 + t71 * t138;
t163 = -t11 + t165;
t161 = t30 * t118;
t160 = t106 * t118;
t12 = t118 * qJD(5) - t110 * t60 + t112 * t61;
t158 = -t12 + t160;
t135 = qJD(1) * qJD(2);
t126 = t113 * t135;
t127 = t111 * t135;
t117 = t107 * t126 - t108 * t127;
t90 = qJ(2) * t144 + qJD(3);
t74 = -pkin(6) * t144 + t90;
t152 = -pkin(6) + qJ(2);
t85 = t152 * t108;
t80 = qJD(1) * t85;
t39 = t111 * t74 + t113 * t80;
t22 = -t39 * qJD(4) + t117;
t15 = t60 * pkin(7) + t22;
t25 = -pkin(7) * t146 + t39;
t123 = t110 * t15 - t25 * t138;
t139 = qJD(4) * t113;
t133 = -t107 * t127 - t108 * t126 - t74 * t139;
t140 = qJD(4) * t111;
t21 = -t80 * t140 - t133;
t14 = -t61 * pkin(7) + t21;
t149 = t111 * t80;
t38 = t113 * t74 - t149;
t24 = -t71 * pkin(7) + t38;
t23 = qJD(4) * pkin(4) + t24;
t1 = (qJD(5) * t23 + t14) * t112 + t123;
t157 = t146 ^ 2;
t156 = t71 ^ 2;
t155 = pkin(4) * t71;
t153 = t71 * t146;
t84 = t152 * t107;
t47 = t111 * t84 + t113 * t85;
t151 = t110 * t25;
t148 = t112 * t25;
t104 = t107 ^ 2;
t105 = t108 ^ 2;
t145 = t104 + t105;
t142 = qJD(2) * t111;
t141 = qJD(2) * t113;
t136 = t107 * qJD(3);
t134 = qJD(1) * qJD(3);
t132 = -t108 * pkin(2) - t107 * qJ(3) - pkin(1);
t129 = -pkin(4) * t106 - t23;
t128 = qJ(2) * t135;
t125 = t107 * t134;
t124 = -t110 * t14 + t112 * t15;
t46 = -t111 * t85 + t113 * t84;
t121 = 0.2e1 * t146;
t115 = qJD(1) ^ 2;
t83 = t145 * t115;
t72 = t108 * pkin(3) - t132;
t119 = t105 * t128;
t6 = t110 * t23 + t148;
t79 = t107 * t113 - t108 * t111;
t28 = -t79 * pkin(7) + t46;
t29 = -t78 * pkin(7) + t47;
t9 = -t110 * t29 + t112 * t28;
t10 = t110 * t28 + t112 * t29;
t41 = -t110 * t78 + t112 * t79;
t82 = t110 * t113 + t112 * t111;
t81 = -t110 * t111 + t112 * t113;
t45 = t61 * pkin(4) + t125;
t26 = t107 * t142 + t108 * t141 + t84 * t139 - t85 * t140;
t2 = -t6 * qJD(5) + t124;
t27 = -t47 * qJD(4) + t107 * t141 - t108 * t142;
t114 = qJD(4) ^ 2;
t91 = t104 * t128;
t68 = t107 * t139 - t108 * t140;
t66 = 0.2e1 * t145 * t135;
t48 = t68 * pkin(4) + t136;
t44 = t106 * t82;
t43 = t106 * t81;
t42 = t78 * pkin(4) + t72;
t40 = t110 * t79 + t112 * t78;
t19 = t67 * pkin(7) + t27;
t18 = -t68 * pkin(7) + t26;
t17 = t41 * qJD(5) - t110 * t67 + t112 * t68;
t16 = t110 * t68 + t112 * t67 + t78 * t137 + t79 * t138;
t8 = t112 * t24 - t151;
t7 = -t110 * t24 - t148;
t5 = t112 * t23 - t151;
t4 = -t10 * qJD(5) - t110 * t18 + t112 * t19;
t3 = t9 * qJD(5) + t110 * t19 + t112 * t18;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0.2e1 * t91 + 0.2e1 * t119, 0, 0, 0, 0, 0, 0, 0.2e1 * t108 * t125, t66, 0.2e1 * t104 * t134, 0.2e1 * t119 + t91 + (t90 * qJD(2) + (-qJD(1) * t132 - t64) * qJD(3)) * t107, -t60 * t79 - t71 * t67, t146 * t67 + t60 * t78 - t79 * t61 - t71 * t68, -t67 * qJD(4), t146 * t68 + t61 * t78, -t68 * qJD(4), 0, t27 * qJD(4) + t121 * t136 + t49 * t68 + t72 * t61, -t26 * qJD(4) - t49 * t67 - t72 * t60 + (qJD(1) * t79 + t71) * t136, -t146 * t26 - t21 * t78 - t22 * t79 - t27 * t71 + t38 * t67 - t39 * t68 + t46 * t60 - t47 * t61, t21 * t47 + t22 * t46 + t39 * t26 + t38 * t27 + (qJD(1) * t72 + t49) * t136, -t11 * t41 - t118 * t16, t11 * t40 - t118 * t17 - t41 * t12 + t159 * t16, -t16 * t106, t12 * t40 + t159 * t17, -t17 * t106, 0, t4 * t106 + t42 * t12 + t159 * t48 + t30 * t17 + t45 * t40, -t3 * t106 - t42 * t11 + t118 * t48 - t30 * t16 + t45 * t41, -t1 * t40 - t10 * t12 + t9 * t11 - t118 * t4 - t159 * t3 + t5 * t16 - t6 * t17 - t2 * t41, t1 * t10 + t2 * t9 + t6 * t3 + t30 * t48 + t5 * t4 + t45 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, -qJ(2) * t83, 0, 0, 0, 0, 0, 0, 0, -t83, 0, -t105 * t115 * qJ(2) + (-qJD(3) - t90) * t144, 0, 0, 0, 0, 0, 0, t89 + (-t71 - t130) * qJD(4), t121 * qJD(4), t156 + t157, -t146 * t39 - t38 * t71 - t125, 0, 0, 0, 0, 0, 0, -t12 - t160, t11 + t165, t162 + t167, -t118 * t5 - t159 * t6 - t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107 * t115 * t108, 0, -t104 * t115, (qJD(2) + t64) * t144, 0, 0, 0, 0, 0, 0, -t114 * t111 - t144 * t146, -t114 * t113 - t71 * t144, -t111 * t61 + t113 * t60 + (t111 * t71 - t113 * t146) * qJD(4), -t49 * t144 + t21 * t111 + t22 * t113 + (-t111 * t38 + t113 * t39) * qJD(4), 0, 0, 0, 0, 0, 0, -t44 * t106 - t144 * t159, -t43 * t106 - t118 * t144, t81 * t11 + t118 * t44 - t82 * t12 - t159 * t43, t1 * t82 - t144 * t30 + t2 * t81 + t6 * t43 - t5 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, t156 - t157, 0, -t153, t89 + (t71 - t130) * qJD(4), 0, -t49 * t71 + t117, t49 * t146 + (t38 + t149) * qJD(4) + t133, 0, 0, t154, t164, t163, -t154, t158, 0, -t159 * t155 - t7 * t106 - t161 + (t110 * t129 - t148) * qJD(5) + t124, -t118 * t155 + t8 * t106 + t166 + (qJD(5) * t129 - t14) * t112 - t123, t6 * t118 + t8 * t159 - t5 * t159 + t7 * t118 + (t11 * t112 - t110 * t12 + (t110 * t118 - t112 * t159) * qJD(5)) * pkin(4), -t5 * t7 - t6 * t8 + (t1 * t110 + t112 * t2 - t30 * t71 + (-t110 * t5 + t112 * t6) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, t164, t163, -t154, t158, 0, t6 * t106 - t161 + t2, t5 * t106 - t1 + t166, 0, 0;];
tauc_reg = t13;
