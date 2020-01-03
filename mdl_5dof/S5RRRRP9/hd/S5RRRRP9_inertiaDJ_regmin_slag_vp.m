% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRP9_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:06:04
% EndTime: 2019-12-31 22:06:09
% DurationCPUTime: 1.53s
% Computational Cost: add. (1838->219), mult. (4544->378), div. (0->0), fcn. (3708->6), ass. (0->111)
t81 = cos(qJ(2));
t123 = t81 * qJD(2);
t78 = sin(qJ(3));
t110 = t78 * t123;
t80 = cos(qJ(3));
t126 = qJD(3) * t80;
t79 = sin(qJ(2));
t147 = t79 * t126 + t110;
t146 = -0.4e1 * t79;
t139 = cos(qJ(4));
t135 = t79 * t80;
t140 = pkin(6) * t78;
t94 = -t81 * pkin(2) - t79 * pkin(7);
t47 = -pkin(1) + t94;
t42 = t80 * t47;
t24 = -pkin(8) * t135 + t42 + (-pkin(3) - t140) * t81;
t134 = t80 * t81;
t62 = pkin(6) * t134;
t130 = t78 * t47 + t62;
t136 = t78 * t79;
t28 = -pkin(8) * t136 + t130;
t77 = sin(qJ(4));
t145 = t139 * t28 + t77 * t24;
t101 = t139 * qJD(4);
t144 = t139 * qJD(3) + t101;
t74 = t79 ^ 2;
t99 = (-t81 ^ 2 + t74) * qJD(2);
t75 = t80 ^ 2;
t129 = t78 ^ 2 - t75;
t100 = t129 * qJD(3);
t143 = qJD(3) + qJD(4);
t93 = pkin(2) * t79 - pkin(7) * t81;
t45 = t93 * qJD(2);
t125 = qJD(3) * t81;
t114 = t78 * t125;
t71 = t79 * qJD(2);
t85 = t80 * t71 + t114;
t19 = t85 * pkin(6) - t47 * t126 - t78 * t45;
t111 = t78 * t71;
t131 = pkin(6) * t111 + t80 * t45;
t141 = pkin(3) * t79;
t11 = (-pkin(8) * t134 + t141) * qJD(2) + (-t62 + (pkin(8) * t79 - t47) * t78) * qJD(3) + t131;
t13 = -pkin(8) * t147 - t19;
t4 = -qJD(4) * t145 + t139 * t11 - t77 * t13;
t82 = 2 * qJD(5);
t142 = -pkin(8) - pkin(7);
t138 = t77 * t78;
t137 = t77 * t80;
t46 = pkin(3) * t136 + t79 * pkin(6);
t127 = qJD(3) * t78;
t124 = qJD(4) * t77;
t122 = t81 * qJD(5);
t121 = -0.2e1 * pkin(1) * qJD(2);
t120 = -0.2e1 * pkin(2) * qJD(3);
t119 = t77 * t136;
t69 = pkin(6) * t123;
t31 = pkin(3) * t147 + t69;
t118 = pkin(4) * t71;
t117 = pkin(3) * t127;
t116 = pkin(3) * t124;
t115 = t79 * t127;
t112 = t80 * t125;
t109 = t78 * t126;
t108 = t79 * t123;
t107 = t80 * t123;
t68 = -t80 * pkin(3) - pkin(2);
t106 = t139 * t80;
t104 = t142 * qJD(3);
t103 = qJD(2) * t139;
t98 = 0.2e1 * t108;
t97 = t78 * t107;
t96 = t142 * t139;
t95 = t81 * t103;
t92 = t78 * t96;
t91 = qJD(3) * t96;
t51 = t142 * t80;
t16 = -qJD(4) * t92 - t104 * t137 - t51 * t124 - t78 * t91;
t30 = t142 * t138 - t139 * t51;
t90 = t16 * t81 + t30 * t71;
t17 = -t51 * t101 - t80 * t91 + (qJD(4) * t142 + t104) * t138;
t29 = -t77 * t51 - t92;
t89 = t17 * t81 - t29 * t71;
t88 = t139 * t24 - t77 * t28;
t44 = t139 * t78 + t137;
t3 = -t24 * t101 - t77 * t11 + t28 * t124 - t139 * t13;
t65 = qJ(5) * t71;
t86 = -t3 + t65;
t83 = t81 * t116 + t4;
t27 = t143 * t44;
t70 = pkin(3) * t101;
t67 = -t139 * pkin(3) - pkin(4);
t64 = t77 * pkin(3) + qJ(5);
t63 = -0.2e1 * t116;
t58 = t70 + qJD(5);
t57 = -0.2e1 * t108;
t43 = -t106 + t138;
t33 = t79 * t106 - t119;
t32 = t44 * t79;
t26 = t143 * t138 - t144 * t80;
t23 = t43 * pkin(4) - t44 * qJ(5) + t68;
t20 = -t130 * qJD(3) + t131;
t18 = t32 * pkin(4) - t33 * qJ(5) + t46;
t15 = t78 * t95 - t77 * t115 - qJD(4) * t119 + (t77 * t123 + t144 * t79) * t80;
t14 = t77 * t110 + t27 * t79 - t80 * t95;
t10 = t81 * pkin(4) - t88;
t9 = -t81 * qJ(5) + t145;
t6 = t27 * pkin(4) + t26 * qJ(5) - t44 * qJD(5) + t117;
t5 = t15 * pkin(4) + t14 * qJ(5) - t33 * qJD(5) + t31;
t2 = -t118 - t4;
t1 = t86 - t122;
t7 = [0, 0, 0, t98, -0.2e1 * t99, 0, 0, 0, t79 * t121, t81 * t121, 0.2e1 * t75 * t108 - 0.2e1 * t74 * t109, 0.2e1 * t74 * t100 + t97 * t146, 0.2e1 * t79 * t114 + 0.2e1 * t80 * t99, 0.2e1 * t79 * t112 - 0.2e1 * t78 * t99, t57, 0.2e1 * t42 * t71 - 0.2e1 * t20 * t81 + 0.2e1 * (t78 * t108 + t74 * t126) * pkin(6), -0.2e1 * t19 * t81 - 0.2e1 * t130 * t71 + 0.2e1 * (-t74 * t127 + t80 * t98) * pkin(6), -0.2e1 * t33 * t14, 0.2e1 * t14 * t32 - 0.2e1 * t33 * t15, 0.2e1 * t14 * t81 + 0.2e1 * t33 * t71, 0.2e1 * t15 * t81 - 0.2e1 * t32 * t71, t57, 0.2e1 * t46 * t15 + 0.2e1 * t31 * t32 - 0.2e1 * t4 * t81 + 0.2e1 * t88 * t71, -0.2e1 * t46 * t14 - 0.2e1 * t145 * t71 - 0.2e1 * t3 * t81 + 0.2e1 * t31 * t33, -0.2e1 * t10 * t71 + 0.2e1 * t18 * t15 + 0.2e1 * t2 * t81 + 0.2e1 * t5 * t32, -0.2e1 * t1 * t32 - 0.2e1 * t10 * t14 - 0.2e1 * t9 * t15 + 0.2e1 * t2 * t33, -0.2e1 * t1 * t81 + 0.2e1 * t18 * t14 - 0.2e1 * t5 * t33 + 0.2e1 * t71 * t9, 0.2e1 * t9 * t1 + 0.2e1 * t10 * t2 + 0.2e1 * t18 * t5; 0, 0, 0, 0, 0, t123, -t71, 0, -t69, pkin(6) * t71, -t79 * t100 + t97, t109 * t146 - t129 * t123, t111 - t112, t85, 0, (pkin(7) * t134 + (-pkin(2) * t80 + t140) * t79) * qJD(3) + (t94 * t78 - t62) * qJD(2), (pkin(6) * t135 + t93 * t78) * qJD(3) + (t81 * t140 + t94 * t80) * qJD(2), -t14 * t44 - t33 * t26, t14 * t43 - t44 * t15 + t26 * t32 - t33 * t27, t26 * t81 + t44 * t71, t27 * t81 - t43 * t71, 0, t32 * t117 + t68 * t15 + t46 * t27 + t31 * t43 + t89, t33 * t117 - t68 * t14 - t46 * t26 + t31 * t44 - t90, t23 * t15 + t18 * t27 + t6 * t32 + t5 * t43 + t89, -t1 * t43 - t10 * t26 - t29 * t14 - t30 * t15 + t16 * t32 + t17 * t33 + t2 * t44 - t9 * t27, t23 * t14 + t18 * t26 - t6 * t33 - t5 * t44 + t90, t1 * t30 + t10 * t17 - t9 * t16 + t18 * t6 + t2 * t29 + t5 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t109, -0.2e1 * t100, 0, 0, 0, t78 * t120, t80 * t120, -0.2e1 * t44 * t26, 0.2e1 * t26 * t43 - 0.2e1 * t44 * t27, 0, 0, 0, 0.2e1 * t43 * t117 + 0.2e1 * t68 * t27, 0.2e1 * t44 * t117 - 0.2e1 * t68 * t26, 0.2e1 * t23 * t27 + 0.2e1 * t6 * t43, 0.2e1 * t16 * t43 + 0.2e1 * t17 * t44 - 0.2e1 * t29 * t26 - 0.2e1 * t30 * t27, 0.2e1 * t23 * t26 - 0.2e1 * t6 * t44, -0.2e1 * t30 * t16 + 0.2e1 * t29 * t17 + 0.2e1 * t23 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107 - t115, -t147, t71, t20, t19, 0, 0, -t14, -t15, t71, t103 * t141 + t83, (t81 * t101 - t77 * t71) * pkin(3) + t3, (pkin(4) - t67) * t71 + t83, t33 * t116 - t67 * t14 - t64 * t15 - t58 * t32, t64 * t71 + (-qJD(5) - t58) * t81 + t86, t1 * t64 + t10 * t116 + t2 * t67 + t9 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, -t127, 0, -pkin(7) * t126, pkin(7) * t127, 0, 0, -t26, -t27, 0, -t17, t16, -t17, t44 * t116 - t67 * t26 - t64 * t27 - t58 * t43, -t16, t116 * t29 - t16 * t64 + t17 * t67 + t30 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -0.2e1 * t70, t63, 0, 0.2e1 * t58, 0.2e1 * t116 * t67 + 0.2e1 * t64 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t15, t71, t4, t3, t4 + 0.2e1 * t118, pkin(4) * t14 - t15 * qJ(5) - t32 * qJD(5), -t3 + 0.2e1 * t65 - 0.2e1 * t122, -t2 * pkin(4) + t1 * qJ(5) + t9 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t27, 0, -t17, t16, -t17, pkin(4) * t26 - t27 * qJ(5) - t43 * qJD(5), -t16, -t17 * pkin(4) - t16 * qJ(5) + t30 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, -t70, -t116, 0, t82 + t70, -pkin(4) * t116 + t58 * qJ(5) + t64 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, qJ(5) * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t14, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
