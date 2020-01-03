% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR11_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:35:13
% EndTime: 2019-12-31 21:35:18
% DurationCPUTime: 1.41s
% Computational Cost: add. (907->192), mult. (2331->352), div. (0->0), fcn. (1798->6), ass. (0->119)
t79 = cos(qJ(2));
t132 = qJD(3) * t79;
t78 = cos(qJ(3));
t116 = t78 * t132;
t76 = sin(qJ(2));
t66 = t76 * qJD(2);
t75 = sin(qJ(3));
t156 = t75 * t66 - t116;
t155 = -0.4e1 * t76;
t72 = t78 ^ 2;
t139 = t75 ^ 2 - t72;
t107 = t139 * qJD(3);
t134 = qJD(3) * t75;
t74 = sin(qJ(5));
t77 = cos(qJ(5));
t43 = t74 * t75 + t77 * t78;
t68 = qJD(3) * t78;
t20 = t43 * qJD(5) - t74 * t134 - t77 * t68;
t147 = pkin(7) * t79;
t99 = pkin(2) * t76 - t147;
t45 = t99 * qJD(2);
t145 = t76 * pkin(7);
t100 = -t79 * pkin(2) - t145;
t47 = -pkin(1) + t100;
t141 = t75 * t45 + t47 * t68;
t117 = t75 * t132;
t84 = t78 * t66 + t117;
t15 = t84 * pkin(6) - t141;
t137 = qJ(4) * t78;
t95 = pkin(3) * t75 - t137;
t88 = pkin(6) + t95;
t26 = t88 * t76;
t128 = t75 * qJD(4);
t31 = t95 * qJD(3) - t128;
t136 = t75 * qJ(4);
t96 = t78 * pkin(3) + t136;
t46 = -pkin(2) - t96;
t154 = (-t46 * t79 + t145) * qJD(2) - qJD(3) * t26 - t31 * t76;
t104 = -t156 * pkin(6) + t47 * t134 - t78 * t45;
t133 = qJD(3) * t76;
t118 = t75 * t133;
t142 = t78 * t79;
t150 = pkin(3) + pkin(4);
t5 = pkin(8) * t118 + (-pkin(8) * t142 - t150 * t76) * qJD(2) + t104;
t144 = t76 * t78;
t124 = qJ(4) * qJD(2);
t65 = t76 * t124;
t6 = t65 + (-pkin(6) * qJD(2) + pkin(8) * qJD(3)) * t144 + (-qJD(4) + (-pkin(6) * qJD(3) + pkin(8) * qJD(2)) * t75) * t79 + t141;
t146 = pkin(8) * t76;
t148 = pkin(6) * t75;
t63 = t79 * t148;
t69 = t79 * pkin(3);
t18 = t79 * pkin(4) + t63 + t69 + (-t47 - t146) * t78;
t64 = pkin(6) * t142;
t140 = t75 * t47 + t64;
t23 = -t79 * qJ(4) + t140;
t19 = t75 * t146 + t23;
t94 = t74 * t18 + t77 * t19;
t2 = -t94 * qJD(5) + t77 * t5 - t74 * t6;
t153 = -t150 * t78 - t136;
t127 = t78 * qJD(4);
t152 = t96 * qJD(3) - t127;
t151 = 0.2e1 * qJD(4);
t149 = pkin(7) - pkin(8);
t143 = t77 * t75;
t71 = t76 ^ 2;
t138 = -t79 ^ 2 + t71;
t135 = qJD(2) * t78;
t131 = qJD(5) * t74;
t130 = qJD(5) * t77;
t129 = qJD(5) * t79;
t126 = t79 * qJD(2);
t125 = t79 * qJD(4);
t123 = -0.2e1 * pkin(1) * qJD(2);
t122 = -0.2e1 * pkin(2) * qJD(3);
t121 = pkin(3) * t66;
t120 = pkin(7) * t134;
t119 = pkin(7) * t68;
t52 = t149 * t78;
t113 = t75 * t126;
t112 = t75 * t68;
t111 = t76 * t126;
t109 = t78 * t126;
t108 = t78 * t47 - t63;
t106 = t138 * qJD(2);
t105 = 0.2e1 * t111;
t103 = t75 * t109;
t102 = qJD(3) * t52;
t101 = t149 * t134;
t24 = -t108 + t69;
t93 = -t23 * t75 + t24 * t78;
t51 = t149 * t75;
t92 = t74 * t51 + t77 * t52;
t91 = t74 * t78 - t143;
t90 = t77 * qJ(4) - t150 * t74;
t1 = -t18 * t130 + t19 * t131 - t74 * t5 - t77 * t6;
t87 = -t150 * t75 + t137;
t85 = -pkin(6) + t87;
t12 = t65 - t15 - t125;
t13 = t104 - t121;
t81 = t93 * qJD(3) + t12 * t78 + t13 * t75;
t57 = -0.2e1 * t111;
t56 = pkin(7) * t116;
t40 = pkin(2) - t153;
t32 = t109 - t118;
t30 = t43 * t76;
t29 = -t76 * t143 + t74 * t144;
t28 = t74 * qJD(4) + t90 * qJD(5);
t27 = qJ(4) * t131 - t77 * qJD(4) + t130 * t150;
t25 = t87 * qJD(3) + t128;
t22 = t85 * t76;
t21 = t75 * t130 - t78 * t131 - t77 * t134 + t74 * t68;
t14 = t88 * t126 + t152 * t76;
t11 = t92 * qJD(5) - t74 * t101 - t77 * t102;
t10 = t77 * t101 - t74 * t102 - t51 * t130 + t52 * t131;
t9 = t43 * t126 + (qJD(3) - qJD(5)) * t76 * t91;
t8 = t74 * t109 - t77 * t113 + t20 * t76;
t7 = (t153 * qJD(3) + t127) * t76 + t85 * t126;
t3 = [0, 0, 0, t105, -0.2e1 * t106, 0, 0, 0, t76 * t123, t79 * t123, 0.2e1 * t72 * t111 - 0.2e1 * t71 * t112, t103 * t155 + 0.2e1 * t71 * t107, 0.2e1 * t76 * t117 + 0.2e1 * t138 * t135, -0.2e1 * t75 * t106 + 0.2e1 * t76 * t116, t57, 0.2e1 * t104 * t79 + 0.2e1 * t108 * t66 + 0.2e1 * (t75 * t105 + t71 * t68) * pkin(6), -0.2e1 * t15 * t79 - 0.2e1 * t140 * t66 + 0.2e1 * (t78 * t105 - t71 * t134) * pkin(6), 0.2e1 * (qJD(2) * t26 * t75 + t13) * t79 + 0.2e1 * (-qJD(2) * t24 + t14 * t75 + t26 * t68) * t76, 0.2e1 * t93 * t126 + 0.2e1 * (-t12 * t75 + t13 * t78 + (-t23 * t78 - t24 * t75) * qJD(3)) * t76, 0.2e1 * (-t26 * t135 - t12) * t79 + 0.2e1 * (qJD(2) * t23 + t26 * t134 - t14 * t78) * t76, 0.2e1 * t23 * t12 + 0.2e1 * t24 * t13 + 0.2e1 * t26 * t14, 0.2e1 * t30 * t9, -0.2e1 * t9 * t29 - 0.2e1 * t30 * t8, -0.2e1 * t30 * t66 + 0.2e1 * t9 * t79, 0.2e1 * t29 * t66 - 0.2e1 * t8 * t79, t57, 0.2e1 * t2 * t79 - 0.2e1 * (t77 * t18 - t74 * t19) * t66 + 0.2e1 * t7 * t29 + 0.2e1 * t22 * t8, 0.2e1 * t1 * t79 + 0.2e1 * t22 * t9 + 0.2e1 * t7 * t30 + 0.2e1 * t66 * t94; 0, 0, 0, 0, 0, t126, -t66, 0, -pkin(6) * t126, pkin(6) * t66, -t76 * t107 + t103, t112 * t155 - t139 * t126, t156, t84, 0, t56 + (-pkin(2) * t78 + t148) * t133 + (t100 * t75 - t64) * qJD(2), (pkin(6) * t144 + t99 * t75) * qJD(3) + (t100 * t78 + t63) * qJD(2), t56 + (t46 * t133 - t14) * t78 - t154 * t75, t81, (-t14 + (t46 * t76 + t147) * qJD(3)) * t75 + t154 * t78, pkin(7) * t81 + t14 * t46 + t26 * t31, -t30 * t20 - t9 * t91, t20 * t29 - t30 * t21 - t9 * t43 + t8 * t91, -t20 * t79 + t66 * t91, -t21 * t79 + t43 * t66, 0, -t11 * t79 - (t77 * t51 - t74 * t52) * t66 + t25 * t29 + t40 * t8 + t7 * t43 + t22 * t21, t10 * t79 - t22 * t20 + t25 * t30 + t40 * t9 + t66 * t92 - t7 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t112, -0.2e1 * t107, 0, 0, 0, t75 * t122, t78 * t122, 0.2e1 * t46 * t134 - 0.2e1 * t31 * t78, 0, -0.2e1 * t31 * t75 - 0.2e1 * t46 * t68, 0.2e1 * t46 * t31, 0.2e1 * t91 * t20, 0.2e1 * t20 * t43 + 0.2e1 * t21 * t91, 0, 0, 0, 0.2e1 * t40 * t21 + 0.2e1 * t25 * t43, -0.2e1 * t40 * t20 - 0.2e1 * t25 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t76 * t68 - t113, t66, -t104, t15, -t104 + 0.2e1 * t121, (-pkin(3) * t126 - qJ(4) * t133) * t78 + (-t79 * t124 + (pkin(3) * qJD(3) - qJD(4)) * t76) * t75, 0.2e1 * t65 - t15 - 0.2e1 * t125, -t13 * pkin(3) + t12 * qJ(4) + t23 * qJD(4), 0, 0, -t9, t8, t66, -t28 * t79 - (-t74 * qJ(4) - t150 * t77) * t66 - t2, t27 * t79 + t66 * t90 - t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, -t134, 0, -t119, t120, -t119, -t152, -t120, -t152 * pkin(7), 0, 0, t20, t21, 0, t11, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, qJ(4) * t151, 0, 0, 0, 0, 0, 0.2e1 * t28, -0.2e1 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, t32, 0, t13, 0, 0, 0, 0, 0, -t129 * t74 - t66 * t77, -t129 * t77 + t66 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, t119, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, -t66, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t21, 0, -t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, -t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
