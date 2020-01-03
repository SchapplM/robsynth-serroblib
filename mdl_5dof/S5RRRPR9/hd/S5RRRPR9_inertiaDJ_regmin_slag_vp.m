% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR9_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:24:38
% EndTime: 2019-12-31 21:24:43
% DurationCPUTime: 1.30s
% Computational Cost: add. (1687->178), mult. (4335->352), div. (0->0), fcn. (3827->8), ass. (0->116)
t102 = sin(qJ(2));
t101 = sin(qJ(3));
t105 = cos(qJ(2));
t131 = t105 * qJD(2);
t124 = t101 * t131;
t104 = cos(qJ(3));
t134 = qJD(3) * t104;
t149 = t102 * t134 + t124;
t148 = -0.4e1 * t102;
t95 = t102 ^ 2;
t118 = (-t105 ^ 2 + t95) * qJD(2);
t96 = t104 ^ 2;
t141 = t101 ^ 2 - t96;
t119 = t141 * qJD(3);
t98 = sin(pkin(9));
t147 = pkin(3) * t98;
t146 = pkin(6) * t101;
t145 = -qJ(4) - pkin(7);
t132 = t104 * qJD(4);
t137 = t104 * t105;
t92 = t102 * qJD(2);
t125 = t101 * t92;
t114 = pkin(2) * t102 - pkin(7) * t105;
t74 = t114 * qJD(2);
t143 = pkin(6) * t125 + t104 * t74;
t115 = -t105 * pkin(2) - t102 * pkin(7);
t77 = -pkin(1) + t115;
t87 = pkin(6) * t137;
t16 = -t102 * t132 + (pkin(3) * t102 - qJ(4) * t137) * qJD(2) + (-t87 + (qJ(4) * t102 - t77) * t101) * qJD(3) + t143;
t138 = t102 * t104;
t144 = -t101 * t74 - t77 * t134;
t20 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t138 + (-qJD(4) * t102 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t105) * t101 - t144;
t99 = cos(pkin(9));
t9 = t98 * t16 + t99 * t20;
t70 = t104 * t77;
t41 = -qJ(4) * t138 + t70 + (-pkin(3) - t146) * t105;
t139 = t101 * t102;
t142 = t101 * t77 + t87;
t44 = -qJ(4) * t139 + t142;
t22 = t98 * t41 + t99 * t44;
t120 = qJD(3) * t145;
t59 = t101 * t120 + t132;
t60 = -t101 * qJD(4) + t104 * t120;
t34 = t99 * t59 + t98 * t60;
t78 = t145 * t101;
t79 = t145 * t104;
t46 = t98 * t78 - t99 * t79;
t75 = pkin(3) * t139 + t102 * pkin(6);
t136 = qJD(3) * t101;
t135 = qJD(3) * t102;
t133 = qJD(3) * t105;
t130 = -0.2e1 * pkin(1) * qJD(2);
t129 = -0.2e1 * pkin(2) * qJD(3);
t90 = pkin(6) * t131;
t48 = t149 * pkin(3) + t90;
t91 = pkin(3) * t136;
t89 = -t104 * pkin(3) - pkin(2);
t128 = t101 * t133;
t126 = t104 * t133;
t123 = t101 * t134;
t122 = t102 * t131;
t121 = t104 * t131;
t8 = t99 * t16 - t98 * t20;
t21 = t99 * t41 - t98 * t44;
t33 = -t98 * t59 + t99 * t60;
t45 = t99 * t78 + t98 * t79;
t117 = 0.2e1 * t122;
t116 = t101 * t121;
t100 = sin(qJ(5));
t103 = cos(qJ(5));
t58 = t99 * t138 - t98 * t139;
t12 = -t105 * pkin(4) - t58 * pkin(8) + t21;
t68 = t99 * t101 + t98 * t104;
t57 = t68 * t102;
t13 = -t57 * pkin(8) + t22;
t113 = t100 * t13 - t103 * t12;
t112 = t100 * t12 + t103 * t13;
t35 = -t68 * pkin(8) + t45;
t109 = t98 * t101 - t99 * t104;
t36 = -pkin(8) * t109 + t46;
t111 = t100 * t36 - t103 * t35;
t110 = t100 * t35 + t103 * t36;
t29 = t100 * t58 + t103 * t57;
t30 = -t100 * t57 + t103 * t58;
t39 = t100 * t68 + t103 * t109;
t40 = -t100 * t109 + t103 * t68;
t88 = t99 * pkin(3) + pkin(4);
t108 = t100 * t88 + t103 * t147;
t107 = t100 * t147 - t103 * t88;
t106 = t104 * t92 + t128;
t32 = -t99 * t121 + t98 * t124 + t68 * t135;
t6 = pkin(4) * t92 + t32 * pkin(8) + t8;
t31 = t109 * t135 - t68 * t131;
t7 = t31 * pkin(8) + t9;
t2 = -t112 * qJD(5) - t100 * t7 + t103 * t6;
t1 = t113 * qJD(5) - t100 * t6 - t103 * t7;
t84 = -0.2e1 * t122;
t63 = t99 * t134 - t98 * t136;
t62 = t68 * qJD(3);
t55 = t108 * qJD(5);
t54 = t107 * qJD(5);
t50 = pkin(4) * t109 + t89;
t47 = t62 * pkin(4) + t91;
t42 = t57 * pkin(4) + t75;
t27 = -t142 * qJD(3) + t143;
t26 = t106 * pkin(6) + t144;
t25 = -t62 * pkin(8) + t34;
t24 = -t63 * pkin(8) + t33;
t23 = -t31 * pkin(4) + t48;
t18 = t40 * qJD(5) + t100 * t63 + t103 * t62;
t17 = -t39 * qJD(5) - t100 * t62 + t103 * t63;
t11 = t30 * qJD(5) - t100 * t32 - t103 * t31;
t10 = -t29 * qJD(5) + t100 * t31 - t103 * t32;
t4 = -t110 * qJD(5) - t100 * t25 + t103 * t24;
t3 = t111 * qJD(5) - t100 * t24 - t103 * t25;
t5 = [0, 0, 0, t117, -0.2e1 * t118, 0, 0, 0, t102 * t130, t105 * t130, 0.2e1 * t96 * t122 - 0.2e1 * t95 * t123, t116 * t148 + 0.2e1 * t95 * t119, 0.2e1 * t102 * t128 + 0.2e1 * t104 * t118, -0.2e1 * t101 * t118 + 0.2e1 * t102 * t126, t84, 0.2e1 * t70 * t92 - 0.2e1 * t27 * t105 + 0.2e1 * (t101 * t122 + t95 * t134) * pkin(6), -0.2e1 * t26 * t105 - 0.2e1 * t142 * t92 + 0.2e1 * (t104 * t117 - t136 * t95) * pkin(6), 0.2e1 * t21 * t32 + 0.2e1 * t22 * t31 - 0.2e1 * t9 * t57 - 0.2e1 * t8 * t58, 0.2e1 * t21 * t8 + 0.2e1 * t22 * t9 + 0.2e1 * t75 * t48, 0.2e1 * t30 * t10, -0.2e1 * t10 * t29 - 0.2e1 * t30 * t11, -0.2e1 * t10 * t105 + 0.2e1 * t30 * t92, 0.2e1 * t11 * t105 - 0.2e1 * t29 * t92, t84, -0.2e1 * t2 * t105 + 0.2e1 * t42 * t11 - 0.2e1 * t113 * t92 + 0.2e1 * t23 * t29, -0.2e1 * t1 * t105 + 0.2e1 * t42 * t10 - 0.2e1 * t112 * t92 + 0.2e1 * t23 * t30; 0, 0, 0, 0, 0, t131, -t92, 0, -t90, pkin(6) * t92, -t102 * t119 + t116, t123 * t148 - t141 * t131, t125 - t126, t106, 0, (pkin(7) * t137 + (-pkin(2) * t104 + t146) * t102) * qJD(3) + (t115 * t101 - t87) * qJD(2), (pkin(6) * t138 + t114 * t101) * qJD(3) + (t115 * t104 + t105 * t146) * qJD(2), -t109 * t9 - t21 * t63 - t22 * t62 + t46 * t31 + t45 * t32 - t33 * t58 - t34 * t57 - t8 * t68, t21 * t33 + t22 * t34 + t8 * t45 + t9 * t46 + t48 * t89 + t75 * t91, t10 * t40 + t30 * t17, -t10 * t39 - t40 * t11 - t17 * t29 - t30 * t18, -t17 * t105 + t40 * t92, t18 * t105 - t39 * t92, 0, -t4 * t105 + t50 * t11 - t111 * t92 + t42 * t18 + t23 * t39 + t47 * t29, t50 * t10 - t3 * t105 - t110 * t92 + t42 * t17 + t23 * t40 + t47 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t123, -0.2e1 * t119, 0, 0, 0, t101 * t129, t104 * t129, -0.2e1 * t109 * t34 - 0.2e1 * t33 * t68 - 0.2e1 * t45 * t63 - 0.2e1 * t46 * t62, 0.2e1 * t45 * t33 + 0.2e1 * t46 * t34 + 0.2e1 * t89 * t91, 0.2e1 * t40 * t17, -0.2e1 * t17 * t39 - 0.2e1 * t40 * t18, 0, 0, 0, 0.2e1 * t50 * t18 + 0.2e1 * t47 * t39, 0.2e1 * t50 * t17 + 0.2e1 * t47 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101 * t135 + t121, -t149, t92, t27, t26, (t31 * t98 + t32 * t99) * pkin(3), (t8 * t99 + t9 * t98) * pkin(3), 0, 0, t10, -t11, t92, t55 * t105 - t107 * t92 + t2, -t54 * t105 - t108 * t92 + t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, -t136, 0, -pkin(7) * t134, pkin(7) * t136, (-t62 * t98 - t63 * t99) * pkin(3), (t33 * t99 + t34 * t98) * pkin(3), 0, 0, t17, -t18, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t55, 0.2e1 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, 0, 0, 0, 0, t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, 0, 0, 0, 0, 0, t18, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t11, t92, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t18, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
