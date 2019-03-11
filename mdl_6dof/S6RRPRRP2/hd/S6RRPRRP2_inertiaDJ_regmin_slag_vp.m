% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:45:12
% EndTime: 2019-03-09 11:45:17
% DurationCPUTime: 1.60s
% Computational Cost: add. (4895->200), mult. (10488->348), div. (0->0), fcn. (10676->8), ass. (0->135)
t142 = cos(pkin(10));
t85 = sin(pkin(10));
t88 = sin(qJ(2));
t90 = cos(qJ(2));
t68 = t142 * t88 + t85 * t90;
t87 = sin(qJ(4));
t120 = t142 * t90;
t104 = t85 * t88 - t120;
t160 = cos(qJ(4));
t96 = t160 * t104;
t47 = t87 * t68 + t96;
t48 = -t87 * t104 + t160 * t68;
t125 = -t90 * pkin(2) - pkin(1);
t57 = t104 * pkin(3) + t125;
t28 = t47 * pkin(4) - t48 * pkin(9) + t57;
t151 = -qJ(3) - pkin(7);
t75 = t151 * t88;
t76 = t151 * t90;
t51 = t142 * t75 + t85 * t76;
t43 = -t68 * pkin(8) + t51;
t52 = -t142 * t76 + t85 * t75;
t44 = -t104 * pkin(8) + t52;
t30 = t160 * t44 + t87 * t43;
t86 = sin(qJ(5));
t89 = cos(qJ(5));
t168 = t86 * t28 + t89 * t30;
t114 = t142 * pkin(2) + pkin(3);
t164 = pkin(2) * t85;
t167 = -t160 * t114 + t87 * t164;
t139 = qJD(4) * t87;
t141 = qJD(2) * t88;
t66 = qJD(2) * t120 - t85 * t141;
t97 = qJD(2) * t68;
t31 = -qJD(4) * t96 - t68 * t139 + t160 * t66 - t87 * t97;
t82 = qJD(5) * t89;
t107 = t86 * t31 + t48 * t82;
t83 = t86 ^ 2;
t84 = t89 ^ 2;
t145 = t83 - t84;
t118 = t145 * qJD(5);
t121 = qJD(4) * t160;
t119 = qJD(2) * t151;
t100 = -t88 * qJD(3) + t90 * t119;
t99 = t90 * qJD(3) + t88 * t119;
t42 = t85 * t100 + t142 * t99;
t91 = -pkin(8) * t97 + t42;
t41 = t142 * t100 - t85 * t99;
t92 = -t66 * pkin(8) + t41;
t11 = -t43 * t121 + t44 * t139 - t160 * t91 - t87 * t92;
t32 = t48 * qJD(4) + t160 * t97 + t87 * t66;
t81 = pkin(2) * t141;
t53 = pkin(3) * t97 + t81;
t16 = t32 * pkin(4) - t31 * pkin(9) + t53;
t5 = -qJD(5) * t168 + t86 * t11 + t89 * t16;
t7 = t47 * qJ(6) + t168;
t111 = t89 * t28 - t86 * t30;
t8 = -t47 * pkin(5) - t111;
t116 = t7 * t89 + t8 * t86;
t136 = t47 * qJD(6);
t143 = t32 * qJ(6);
t137 = qJD(5) * t86;
t4 = t89 * t11 + t30 * t137 - t86 * t16 - t28 * t82;
t2 = t136 - t4 + t143;
t161 = t32 * pkin(5);
t3 = -t161 - t5;
t166 = qJD(5) * t116 + t2 * t86 - t3 * t89;
t165 = 0.2e1 * qJD(6);
t163 = pkin(9) * t32;
t162 = pkin(9) * t47;
t93 = t87 * t114 + t160 * t164;
t62 = pkin(9) + t93;
t159 = t47 * t62;
t158 = t48 * t31;
t157 = t48 * t86;
t156 = t48 * t89;
t154 = t86 * t32;
t153 = t89 * t31;
t152 = t89 * t32;
t12 = t44 * t121 + t43 * t139 - t160 * t92 + t87 * t91;
t29 = -t160 * t43 + t87 * t44;
t150 = t12 * t86 + t29 * t82;
t59 = t93 * qJD(4);
t135 = t86 * qJD(6);
t64 = -pkin(5) * t137 + qJ(6) * t82 + t135;
t45 = t59 - t64;
t148 = -t45 + t64;
t61 = -pkin(4) + t167;
t147 = t59 * t86 + t61 * t82;
t146 = -t83 - t84;
t144 = qJ(6) * t89;
t140 = qJD(2) * t90;
t112 = pkin(5) * t86 - t144;
t18 = t112 * t48 + t29;
t138 = qJD(5) * t18;
t134 = t89 * qJD(6);
t133 = -0.2e1 * pkin(1) * qJD(2);
t132 = pkin(4) * t137;
t131 = pkin(4) * t82;
t130 = pkin(9) * t137;
t129 = pkin(9) * t82;
t128 = t48 * t137;
t126 = t86 * t82;
t58 = t167 * qJD(4);
t36 = t146 * t58;
t124 = -0.4e1 * t86 * t156;
t122 = t61 * t137 - t59 * t89;
t115 = -t7 * t86 + t8 * t89;
t113 = -t89 * pkin(5) - t86 * qJ(6);
t109 = -t32 * t62 + t47 * t58;
t108 = -t48 * t61 + t159;
t106 = t128 - t153;
t22 = t47 * t82 + t154;
t105 = t47 * t137 - t152;
t74 = -pkin(4) + t113;
t103 = t31 * t74 - t48 * t64 - t163;
t6 = t107 * pkin(5) + qJ(6) * t128 - t48 * t134 - t31 * t144 + t12;
t101 = -t6 + (t48 * t74 - t162) * qJD(5);
t50 = t113 + t61;
t98 = -t6 + (t48 * t50 - t159) * qJD(5);
t95 = t31 * t50 + t45 * t48 + t109;
t94 = t31 * t61 + t48 * t59 + t109;
t63 = t113 * qJD(5) + t134;
t1 = t115 * qJD(5) + t2 * t89 + t3 * t86;
t77 = 0.2e1 * t126;
t72 = -0.2e1 * t118;
t67 = t74 * t137;
t49 = t50 * t137;
t46 = t48 ^ 2;
t38 = -t86 * t58 + t62 * t82;
t37 = t62 * t137 + t89 * t58;
t24 = t29 * t137;
t19 = -t48 * t118 + t86 * t153;
t17 = t18 * t137;
t13 = qJD(5) * t124 - t145 * t31;
t9 = [0, 0, 0, 0.2e1 * t88 * t140, 0.2e1 * (-t88 ^ 2 + t90 ^ 2) * qJD(2), 0, 0, 0, t88 * t133, t90 * t133, -0.2e1 * t104 * t42 - 0.2e1 * t41 * t68 - 0.2e1 * t51 * t66 - 0.2e1 * t52 * t97, 0.2e1 * t125 * t81 + 0.2e1 * t51 * t41 + 0.2e1 * t52 * t42, 0.2e1 * t158, -0.2e1 * t31 * t47 - 0.2e1 * t48 * t32, 0, 0, 0, 0.2e1 * t57 * t32 + 0.2e1 * t53 * t47, 0.2e1 * t57 * t31 + 0.2e1 * t53 * t48, -0.2e1 * t126 * t46 + 0.2e1 * t84 * t158, 0.2e1 * t46 * t118 + t31 * t124, -0.2e1 * t106 * t47 + 0.2e1 * t152 * t48, -0.2e1 * t107 * t47 - 0.2e1 * t154 * t48, 0.2e1 * t47 * t32, 0.2e1 * t107 * t29 + 0.2e1 * t111 * t32 + 0.2e1 * t12 * t157 + 0.2e1 * t5 * t47, -0.2e1 * t106 * t29 + 0.2e1 * t12 * t156 - 0.2e1 * t168 * t32 + 0.2e1 * t4 * t47, 0.2e1 * t107 * t18 + 0.2e1 * t6 * t157 - 0.2e1 * t3 * t47 - 0.2e1 * t8 * t32, 0.2e1 * t115 * t31 - 0.2e1 * t166 * t48, 0.2e1 * t106 * t18 - 0.2e1 * t6 * t156 + 0.2e1 * t2 * t47 + 0.2e1 * t7 * t32, 0.2e1 * t18 * t6 + 0.2e1 * t7 * t2 + 0.2e1 * t8 * t3; 0, 0, 0, 0, 0, t140, -t141, 0, -pkin(7) * t140, pkin(7) * t141 (-t142 * t66 - t85 * t97) * pkin(2) (t142 * t41 + t42 * t85) * pkin(2), 0, 0, t31, -t32, 0, -t12, t11, t19, t13, t22, -t105, 0, t24 + (-qJD(5) * t108 - t12) * t89 + t94 * t86, t108 * t137 + t89 * t94 + t150, t86 * t95 + t89 * t98 + t17, t1, t98 * t86 + (-t95 - t138) * t89, t1 * t62 - t116 * t58 + t18 * t45 + t6 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t59, 0.2e1 * t58, t77, t72, 0, 0, 0, 0.2e1 * t122, 0.2e1 * t147, -0.2e1 * t45 * t89 + 0.2e1 * t49, 0.2e1 * t36, -0.2e1 * t45 * t86 - 0.2e1 * t50 * t82, 0.2e1 * t36 * t62 + 0.2e1 * t50 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, 0, 0, 0, 0, 0, t32, t31, 0, 0, 0, 0, 0, -t105, -t22, -t105, t146 * t31, t22, t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t32, 0, -t12, t11, t19, t13, t22, -t105, 0, t24 + (-pkin(4) * t31 - t163) * t86 + (-t12 + (-pkin(4) * t48 - t162) * qJD(5)) * t89, pkin(4) * t106 + pkin(9) * t105 + t150, t101 * t89 + t103 * t86 + t17, t1, t101 * t86 + (-t103 - t138) * t89, pkin(9) * t1 - t18 * t64 + t6 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t58, t77, t72, 0, 0, 0, t122 - t132, -t131 + t147, t148 * t89 + t49 + t67, t36, t148 * t86 + (-t50 - t74) * t82, pkin(9) * t36 + t45 * t74 - t50 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t72, 0, 0, 0, -0.2e1 * t132, -0.2e1 * t131, 0.2e1 * t64 * t89 + 0.2e1 * t67, 0, 0.2e1 * t64 * t86 - 0.2e1 * t74 * t82, -0.2e1 * t74 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t107, t32, t5, t4, t5 + 0.2e1 * t161, t113 * t31 + (qJD(5) * t112 - t135) * t48, 0.2e1 * t136 - t4 + 0.2e1 * t143, -t3 * pkin(5) + t2 * qJ(6) + t7 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, -t137, 0, -t38, t37, -t38, t63, -t37, t112 * t58 + t62 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137, -t82, -t137, 0, t82, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, -t137, 0, -t129, t130, -t129, t63, -t130, t63 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, qJ(6) * t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t106, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
