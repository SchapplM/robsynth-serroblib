% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:01:16
% EndTime: 2019-03-09 06:01:20
% DurationCPUTime: 1.54s
% Computational Cost: add. (2012->203), mult. (4741->356), div. (0->0), fcn. (4121->8), ass. (0->111)
t98 = cos(qJ(3));
t126 = t98 * qJD(3);
t94 = sin(qJ(4));
t114 = t94 * t126;
t97 = cos(qJ(4));
t130 = qJD(4) * t97;
t95 = sin(qJ(3));
t153 = t95 * t130 + t114;
t152 = -0.4e1 * t95;
t142 = t95 * t97;
t82 = sin(pkin(10)) * pkin(1) + pkin(7);
t145 = t82 * t94;
t102 = -t98 * pkin(3) - t95 * pkin(8);
t83 = -cos(pkin(10)) * pkin(1) - pkin(2);
t56 = t102 + t83;
t45 = t97 * t56;
t29 = -pkin(9) * t142 + t45 + (-pkin(4) - t145) * t98;
t143 = t94 * t95;
t44 = t94 * t56;
t140 = t97 * t98;
t64 = t82 * t140;
t151 = -t44 - t64;
t34 = -pkin(9) * t143 - t151;
t96 = cos(qJ(5));
t30 = t96 * t34;
t93 = sin(qJ(5));
t137 = t93 * t29 + t30;
t149 = pkin(8) + pkin(9);
t73 = t149 * t94;
t74 = t149 * t97;
t134 = -t93 * t73 + t96 * t74;
t110 = t97 * t126;
t131 = qJD(4) * t94;
t150 = -t95 * t131 + t110;
t89 = t95 ^ 2;
t104 = (-t98 ^ 2 + t89) * qJD(3);
t90 = t97 ^ 2;
t133 = t94 ^ 2 - t90;
t105 = t133 * qJD(4);
t125 = qJD(4) + qJD(5);
t129 = qJD(4) * t98;
t117 = t94 * t129;
t87 = t95 * qJD(3);
t48 = t97 * t87 + t117;
t101 = pkin(3) * t95 - pkin(8) * t98;
t67 = t101 * qJD(3);
t18 = -t56 * t130 + t48 * t82 - t94 * t67;
t123 = t93 * t143;
t21 = -qJD(5) * t123 + (t125 * t142 + t114) * t96 + t150 * t93;
t148 = t21 * pkin(5);
t147 = t96 * pkin(4);
t60 = t93 * t97 + t96 * t94;
t37 = t125 * t60;
t20 = -t96 * t110 + t93 * t114 + t37 * t95;
t141 = t96 * t97;
t43 = t95 * t141 - t123;
t146 = t43 * t20;
t144 = t93 * t94;
t59 = -t141 + t144;
t138 = t20 * t59 - t43 * t37;
t111 = t82 * t87;
t135 = t94 * t111 + t97 * t67;
t46 = pkin(4) * t143 + t95 * t82;
t128 = qJD(5) * t93;
t127 = qJD(5) * t96;
t124 = -0.2e1 * pkin(3) * qJD(4);
t122 = 0.2e1 * qJD(3) * t83;
t68 = t82 * t126;
t38 = t153 * pkin(4) + t68;
t121 = pkin(4) * t131;
t120 = pkin(4) * t128;
t119 = pkin(4) * t127;
t115 = t97 * t129;
t113 = t94 * t130;
t112 = t95 * t126;
t86 = -t97 * pkin(4) - pkin(3);
t109 = qJD(4) * t149;
t12 = (pkin(4) * t95 - pkin(9) * t140) * qJD(3) + (-t64 + (pkin(9) * t95 - t56) * t94) * qJD(4) + t135;
t16 = -pkin(9) * t153 - t18;
t108 = t96 * t12 - t93 * t16;
t107 = t96 * t29 - t93 * t34;
t106 = -t96 * t73 - t93 * t74;
t103 = t94 * t110;
t36 = t125 * t144 - t97 * t127 - t96 * t130;
t42 = t60 * t95;
t100 = -t60 * t21 + t42 * t36;
t99 = -t98 * t37 + t59 * t87;
t3 = -t93 * t12 - t29 * t127 + t34 * t128 - t96 * t16;
t65 = t94 * t109;
t66 = t97 * t109;
t22 = t73 * t127 + t74 * t128 + t96 * t65 + t93 * t66;
t4 = -t137 * qJD(5) + t108;
t23 = -t134 * qJD(5) + t93 * t65 - t96 * t66;
t85 = pkin(5) + t147;
t78 = -0.2e1 * t112;
t50 = t94 * t87 - t115;
t40 = t59 * pkin(5) + t86;
t35 = t42 * pkin(5) + t46;
t33 = t37 * pkin(5) + t121;
t32 = -t59 * qJ(6) + t134;
t31 = -t60 * qJ(6) + t106;
t24 = t36 * t98 + t60 * t87;
t19 = t151 * qJD(4) + t135;
t9 = t38 + t148;
t8 = t36 * qJ(6) - t60 * qJD(6) + t23;
t7 = -t37 * qJ(6) - t59 * qJD(6) - t22;
t6 = -t42 * qJ(6) + t137;
t5 = -t98 * pkin(5) - t43 * qJ(6) + t107;
t2 = -t21 * qJ(6) - t42 * qJD(6) - t3;
t1 = pkin(5) * t87 + t20 * qJ(6) - t43 * qJD(6) + t4;
t10 = [0, 0, 0, 0, 0.2e1 * t112, -0.2e1 * t104, 0, 0, 0, t95 * t122, t98 * t122, 0.2e1 * t90 * t112 - 0.2e1 * t89 * t113, t103 * t152 + 0.2e1 * t89 * t105, 0.2e1 * t97 * t104 + 0.2e1 * t95 * t117, -0.2e1 * t94 * t104 + 0.2e1 * t95 * t115, t78, 0.2e1 * t45 * t87 - 0.2e1 * t19 * t98 + 0.2e1 * (t94 * t112 + t89 * t130) * t82, -0.2e1 * t89 * t82 * t131 - 0.2e1 * t18 * t98 + 0.2e1 * (-t44 + t64) * t87, -0.2e1 * t146, 0.2e1 * t42 * t20 - 0.2e1 * t43 * t21, 0.2e1 * t20 * t98 + 0.2e1 * t43 * t87, 0.2e1 * t98 * t21 - 0.2e1 * t42 * t87, t78, 0.2e1 * t107 * t87 + 0.2e1 * t46 * t21 + 0.2e1 * t38 * t42 - 0.2e1 * t4 * t98, -0.2e1 * t137 * t87 - 0.2e1 * t46 * t20 - 0.2e1 * t3 * t98 + 0.2e1 * t38 * t43, -0.2e1 * t1 * t43 - 0.2e1 * t2 * t42 + 0.2e1 * t5 * t20 - 0.2e1 * t6 * t21, 0.2e1 * t5 * t1 + 0.2e1 * t6 * t2 + 0.2e1 * t35 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t42 + t2 * t43 - t6 * t20 - t5 * t21 + t35 * t87 - t9 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t42 * t21 - 0.2e1 * t112 - 0.2e1 * t146; 0, 0, 0, 0, 0, 0, t126, -t87, 0, -t68, t111, -t95 * t105 + t103, t113 * t152 - t133 * t126, t50, t48, 0 (pkin(8) * t140 + (-pkin(3) * t97 + t145) * t95) * qJD(4) + (t102 * t94 - t64) * qJD(3) (t101 * t94 + t82 * t142) * qJD(4) + (t102 * t97 + t98 * t145) * qJD(3), -t20 * t60 - t43 * t36, t100 + t138, t24, -t99, 0, t106 * t87 + t42 * t121 + t86 * t21 - t23 * t98 + t46 * t37 + t38 * t59, t43 * t121 - t134 * t87 - t86 * t20 - t22 * t98 - t46 * t36 + t38 * t60, -t1 * t60 - t2 * t59 + t31 * t20 - t32 * t21 + t5 * t36 - t6 * t37 - t7 * t42 - t8 * t43, t1 * t31 + t2 * t32 + t35 * t33 + t9 * t40 + t5 * t8 + t6 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t126, 0, 0, 0, 0, 0, -t48, t50, 0, 0, 0, 0, 0, t99, t24, -t100 + t138, -t20 * t32 - t21 * t31 - t98 * t33 + t40 * t87 - t42 * t8 + t43 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t113, -0.2e1 * t105, 0, 0, 0, t94 * t124, t97 * t124, -0.2e1 * t60 * t36, 0.2e1 * t36 * t59 - 0.2e1 * t60 * t37, 0, 0, 0, 0.2e1 * t59 * t121 + 0.2e1 * t86 * t37, 0.2e1 * t60 * t121 - 0.2e1 * t86 * t36, 0.2e1 * t31 * t36 - 0.2e1 * t32 * t37 - 0.2e1 * t7 * t59 - 0.2e1 * t8 * t60, 0.2e1 * t31 * t8 + 0.2e1 * t32 * t7 + 0.2e1 * t40 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, -t153, t87, t19, t18, 0, 0, -t20, -t21, t87, t87 * t147 + (-t30 + (t98 * pkin(4) - t29) * t93) * qJD(5) + t108 (t98 * t127 - t93 * t87) * pkin(4) + t3, t85 * t20 + (-t21 * t93 + (-t42 * t96 + t43 * t93) * qJD(5)) * pkin(4), t1 * t85 + (t2 * t93 + (-t5 * t93 + t6 * t96) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t153, -t150, 0, 0, 0, 0, 0, -t21, t20, 0, -t21 * t85 + (-t20 * t93 + (t42 * t93 + t43 * t96) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, -t131, 0, -pkin(8) * t130, pkin(8) * t131, 0, 0, -t36, -t37, 0, t23, t22, t85 * t36 + (-t37 * t93 + (-t59 * t96 + t60 * t93) * qJD(5)) * pkin(4), t8 * t85 + (t7 * t93 + (-t31 * t93 + t32 * t96) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t120, -0.2e1 * t119, 0, 0.2e1 * (-t85 + t147) * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t21, t87, t4, t3, pkin(5) * t20, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t20, 0, -t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t37, 0, t23, t22, pkin(5) * t36, t8 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, -t119, 0, -pkin(5) * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
