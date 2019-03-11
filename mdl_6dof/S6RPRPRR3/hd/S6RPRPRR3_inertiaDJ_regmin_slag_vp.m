% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:42:26
% EndTime: 2019-03-09 03:42:30
% DurationCPUTime: 1.34s
% Computational Cost: add. (1893->192), mult. (4619->347), div. (0->0), fcn. (4378->10), ass. (0->120)
t102 = sin(qJ(5));
t104 = cos(qJ(5));
t139 = qJD(5) * t104;
t99 = cos(pkin(11));
t146 = t104 * t99;
t156 = pkin(8) + qJ(4);
t98 = sin(pkin(11));
t80 = t156 * t98;
t81 = t156 * t99;
t28 = (qJD(4) * t98 + qJD(5) * t81) * t102 - qJD(4) * t146 + t80 * t139;
t92 = -t99 * pkin(4) - pkin(3);
t162 = 0.2e1 * t92;
t147 = t104 * t98;
t73 = t102 * t99 + t147;
t65 = t73 * qJD(5);
t161 = pkin(5) * t65;
t160 = cos(qJ(6));
t103 = sin(qJ(3));
t159 = pkin(3) * t103;
t105 = cos(qJ(3));
t158 = t105 * pkin(5);
t90 = sin(pkin(10)) * pkin(1) + pkin(7);
t157 = t90 * t98;
t148 = t103 * t99;
t143 = t103 * qJ(4);
t91 = -cos(pkin(10)) * pkin(1) - pkin(2);
t66 = -t105 * pkin(3) - t143 + t91;
t60 = t99 * t66;
t36 = -pkin(8) * t148 + t60 + (-pkin(4) - t157) * t105;
t149 = t103 * t98;
t144 = t105 * t99;
t74 = t90 * t144;
t49 = t98 * t66 + t74;
t41 = -pkin(8) * t149 + t49;
t155 = t102 * t36 + t104 * t41;
t93 = qJD(3) * t103;
t130 = t90 * t93;
t62 = -t103 * qJD(4) + (-qJ(4) * t105 + t159) * qJD(3);
t44 = t98 * t130 + t99 * t62;
t153 = -t102 * t80 + t104 * t81;
t142 = qJD(3) * t105;
t129 = t98 * t142;
t75 = t90 * t142;
t56 = pkin(4) * t129 + t75;
t61 = pkin(4) * t149 + t103 * t90;
t152 = t102 * t148 + t103 * t147;
t151 = t98 ^ 2 + t99 ^ 2;
t150 = t103 ^ 2 - t105 ^ 2;
t145 = t105 * t98;
t141 = qJD(4) * t105;
t140 = qJD(5) * t102;
t101 = sin(qJ(6));
t138 = qJD(6) * t101;
t137 = 0.2e1 * qJD(3) * t91;
t136 = t90 * t145;
t135 = t160 * pkin(5);
t126 = t104 * t142;
t127 = t102 * t142;
t106 = -t98 * t126 - t99 * t127 - t139 * t148 + t140 * t149;
t22 = (pkin(4) * t103 - pkin(8) * t144) * qJD(3) + t44;
t54 = t98 * t62;
t35 = t54 + (-pkin(8) * t145 - t90 * t148) * qJD(3);
t8 = -t102 * t22 - t104 * t35 - t36 * t139 + t41 * t140;
t5 = t106 * pkin(9) - t8;
t134 = t160 * t5;
t133 = pkin(5) * t93;
t132 = pkin(5) * t138;
t14 = -t152 * pkin(9) + t155;
t131 = t160 * t14;
t40 = t103 * t65 - t99 * t126 + t98 * t127;
t9 = -t155 * qJD(5) - t102 * t35 + t104 * t22;
t4 = t40 * pkin(9) + t133 + t9;
t128 = -t101 * t5 + t160 * t4;
t125 = t103 * t142;
t124 = t151 * t105;
t123 = -t102 * t41 + t104 * t36;
t122 = -t102 * t81 - t104 * t80;
t121 = t151 * qJD(4);
t120 = 0.2e1 * t125;
t119 = qJD(6) * t135;
t118 = 0.2e1 * t121;
t117 = t160 * t152;
t45 = -t99 * t130 + t54;
t116 = -t44 * t98 + t45 * t99;
t48 = t60 - t136;
t115 = -t48 * t98 + t49 * t99;
t72 = t102 * t98 - t146;
t58 = t72 * t103;
t13 = t58 * pkin(9) + t123 - t158;
t113 = t101 * t14 - t160 * t13;
t112 = t101 * t13 + t131;
t37 = -t73 * pkin(9) + t122;
t38 = -t72 * pkin(9) + t153;
t111 = t101 * t38 - t160 * t37;
t110 = t101 * t37 + t160 * t38;
t109 = -t101 * t73 - t160 * t72;
t47 = -t101 * t72 + t160 * t73;
t64 = t72 * qJD(5);
t16 = t47 * qJD(6) - t101 * t64 + t160 * t65;
t108 = -t105 * t16 - t109 * t93;
t107 = -t105 * t65 + t72 * t93;
t34 = -t101 * t152 - t160 * t58;
t29 = -t73 * qJD(4) - t153 * qJD(5);
t86 = -0.2e1 * t125;
t53 = t72 * pkin(5) + t92;
t43 = t152 * pkin(5) + t61;
t42 = t64 * t105 + t73 * t93;
t33 = -t101 * t58 + t117;
t19 = t64 * pkin(9) + t29;
t18 = -t65 * pkin(9) - t28;
t17 = -t106 * pkin(5) + t56;
t15 = t109 * qJD(6) - t101 * t65 - t160 * t64;
t12 = -t15 * t105 + t47 * t93;
t11 = t34 * qJD(6) - t101 * t40 - t160 * t106;
t10 = qJD(6) * t117 - t101 * t106 - t58 * t138 + t160 * t40;
t7 = -t110 * qJD(6) - t101 * t18 + t160 * t19;
t6 = t111 * qJD(6) - t101 * t19 - t160 * t18;
t2 = -t112 * qJD(6) + t128;
t1 = t113 * qJD(6) - t101 * t4 - t134;
t3 = [0, 0, 0, 0, t120, -0.2e1 * t150 * qJD(3), 0, 0, 0, t103 * t137, t105 * t137, -0.2e1 * t44 * t105 + 0.2e1 * (t48 + 0.2e1 * t136) * t93, 0.2e1 * t45 * t105 + 0.2e1 * (-t49 + 0.2e1 * t74) * t93, 0.2e1 * (-t44 * t99 - t45 * t98) * t103 + 0.2e1 * (-t48 * t99 - t49 * t98) * t142, 0.2e1 * t125 * t90 ^ 2 + 0.2e1 * t48 * t44 + 0.2e1 * t49 * t45, 0.2e1 * t58 * t40, -0.2e1 * t58 * t106 + 0.2e1 * t40 * t152, 0.2e1 * t40 * t105 - 0.2e1 * t58 * t93, -0.2e1 * t105 * t106 - 0.2e1 * t152 * t93, t86, -0.2e1 * t9 * t105 - 0.2e1 * t61 * t106 + 0.2e1 * t123 * t93 + 0.2e1 * t56 * t152, -0.2e1 * t8 * t105 - 0.2e1 * t155 * t93 - 0.2e1 * t61 * t40 - 0.2e1 * t56 * t58, -0.2e1 * t34 * t10, 0.2e1 * t10 * t33 - 0.2e1 * t34 * t11, 0.2e1 * t10 * t105 + 0.2e1 * t34 * t93, 0.2e1 * t105 * t11 - 0.2e1 * t33 * t93, t86, -0.2e1 * t2 * t105 + 0.2e1 * t43 * t11 - 0.2e1 * t113 * t93 + 0.2e1 * t17 * t33, -0.2e1 * t1 * t105 - 0.2e1 * t43 * t10 - 0.2e1 * t112 * t93 + 0.2e1 * t17 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116 * t103 + (t115 * t105 + t150 * t90) * qJD(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-0.1e1 + t151) * t120, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t142, -t93, 0, -t75, t130, t98 * t141 + (-t98 * t143 + (-pkin(3) * t98 - t90 * t99) * t105) * qJD(3), t99 * t141 + (-t99 * t143 + (-pkin(3) * t99 + t157) * t105) * qJD(3), t116, -pkin(3) * t75 + qJ(4) * t116 + qJD(4) * t115, -t40 * t73 + t58 * t64, t73 * t106 + t64 * t152 + t40 * t72 + t58 * t65, t42, -t107, 0, -t29 * t105 - t106 * t92 + t122 * t93 + t56 * t72 + t61 * t65, -t28 * t105 - t153 * t93 - t92 * t40 + t56 * t73 - t61 * t64, -t10 * t47 + t34 * t15, -t10 * t109 - t47 * t11 - t15 * t33 - t34 * t16, t12, -t108, 0, -t7 * t105 - t109 * t17 + t53 * t11 - t111 * t93 + t43 * t16 + t33 * t161, -t53 * t10 - t6 * t105 - t110 * t93 + t43 * t15 + t34 * t161 + t17 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, -t142, -t99 * t93, t98 * t93, qJD(3) * t124, t103 * t121 + (qJ(4) * t124 - t159) * qJD(3), 0, 0, 0, 0, 0, t107, t42, 0, 0, 0, 0, 0, t108, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, qJ(4) * t118, -0.2e1 * t73 * t64, 0.2e1 * t64 * t72 - 0.2e1 * t73 * t65, 0, 0, 0, t65 * t162, -t64 * t162, 0.2e1 * t47 * t15, 0.2e1 * t109 * t15 - 0.2e1 * t47 * t16, 0, 0, 0, -0.2e1 * t109 * t161 + 0.2e1 * t53 * t16, 0.2e1 * t53 * t15 + 0.2e1 * t47 * t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, t99 * t142, 0, t75, 0, 0, 0, 0, 0, -t106, -t40, 0, 0, 0, 0, 0, t11, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t64, 0, 0, 0, 0, 0, t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, t106, t93, t9, t8, 0, 0, -t10, -t11, t93, t160 * t133 + (-t131 + (-t13 + t158) * t101) * qJD(6) + t128, -t134 + (-t4 - t133) * t101 + (t105 * t135 + t113) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t40, 0, 0, 0, 0, 0, -t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t65, 0, t29, t28, 0, 0, t15, -t16, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t132, -0.2e1 * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, t93, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t16, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132, -t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
