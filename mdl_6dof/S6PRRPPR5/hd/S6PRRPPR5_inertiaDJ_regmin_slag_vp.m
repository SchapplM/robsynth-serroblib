% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPPR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:21:15
% EndTime: 2019-03-08 21:21:18
% DurationCPUTime: 1.09s
% Computational Cost: add. (798->171), mult. (2175->329), div. (0->0), fcn. (1956->10), ass. (0->103)
t120 = pkin(4) + pkin(8);
t74 = sin(qJ(3));
t107 = t74 * qJ(4);
t72 = -pkin(3) - qJ(5);
t77 = cos(qJ(3));
t84 = -t72 * t77 + t107;
t68 = sin(pkin(11));
t70 = cos(pkin(11));
t73 = sin(qJ(6));
t76 = cos(qJ(6));
t123 = -t73 * t68 + t76 * t70;
t51 = (t68 ^ 2 + t70 ^ 2) * qJD(5);
t104 = qJD(4) * t77;
t122 = t84 * qJD(3) + qJD(5) * t74 - t104;
t121 = 0.2e1 * qJD(4);
t119 = -pkin(9) + t72;
t71 = cos(pkin(6));
t78 = cos(qJ(2));
t105 = qJD(2) * t78;
t69 = sin(pkin(6));
t95 = t69 * t105;
t75 = sin(qJ(2));
t115 = t69 * t75;
t99 = t74 * t115;
t25 = -qJD(3) * t99 + (qJD(3) * t71 + t95) * t77;
t39 = t77 * t115 + t71 * t74;
t118 = t39 * t25;
t117 = t68 * t74;
t116 = t68 * t77;
t114 = t69 * t78;
t113 = t70 * t74;
t112 = t70 * t77;
t102 = t74 * qJD(3);
t93 = pkin(3) * t102 - t74 * qJD(4);
t27 = -t77 * qJD(5) + (-qJ(4) * t77 + qJ(5) * t74) * qJD(3) + t93;
t63 = t77 * qJD(3);
t61 = pkin(8) * t63;
t48 = pkin(4) * t63 + t61;
t12 = t70 * t27 + t68 * t48;
t42 = -pkin(2) - t84;
t55 = t120 * t74;
t21 = t70 * t42 + t68 * t55;
t56 = t120 * t77;
t106 = qJD(2) * t75;
t103 = qJD(6) * t77;
t101 = qJ(4) * qJD(4);
t100 = -0.2e1 * pkin(2) * qJD(3);
t98 = pkin(8) * t102;
t97 = t68 * t102;
t96 = t69 * t106;
t94 = t70 * t102;
t11 = -t68 * t27 + t70 * t48;
t92 = -t77 * pkin(3) - t107;
t5 = t11 * t70 + t12 * t68;
t26 = t39 * qJD(3) + t74 * t95;
t13 = t26 * t70 - t68 * t96;
t14 = t26 * t68 + t70 * t96;
t6 = t13 * t70 + t14 * t68;
t46 = t70 * t55;
t15 = t74 * pkin(5) + t46 + (pkin(9) * t77 - t42) * t68;
t18 = -pkin(9) * t112 + t21;
t91 = t76 * t15 - t73 * t18;
t90 = t73 * t15 + t76 * t18;
t38 = -t71 * t77 + t99;
t23 = t68 * t114 + t38 * t70;
t24 = -t70 * t114 + t38 * t68;
t89 = t76 * t23 - t73 * t24;
t88 = t73 * t23 + t76 * t24;
t49 = t119 * t68;
t50 = t119 * t70;
t87 = t76 * t49 + t73 * t50;
t86 = t73 * t49 - t76 * t50;
t43 = t76 * t68 + t73 * t70;
t83 = t25 * qJ(4) + t39 * qJD(4);
t37 = t123 * qJD(6);
t81 = -t37 * t74 - t43 * t63;
t80 = t92 * qJD(3) + t104;
t79 = t25 * t77 + t26 * t74 + (t38 * t77 - t39 * t74) * qJD(3);
t60 = t68 * pkin(5) + qJ(4);
t58 = 0.2e1 * t74 * t63;
t52 = -pkin(2) + t92;
t47 = t120 * t102;
t36 = t43 * qJD(6);
t35 = pkin(5) * t112 + t56;
t34 = -qJ(4) * t63 + t93;
t32 = t43 * t77;
t31 = t123 * t77;
t30 = (-pkin(5) * t70 - t120) * t102;
t29 = (t78 * t102 + t77 * t106) * t69;
t28 = (t74 * t106 - t78 * t63) * t69;
t20 = -t68 * t42 + t46;
t19 = t123 * t63 - t36 * t74;
t17 = t43 * t103 - t73 * t97 + t76 * t94;
t16 = t43 * t102 - t103 * t123;
t10 = -qJD(5) * t123 - t87 * qJD(6);
t9 = t43 * qJD(5) + t86 * qJD(6);
t8 = pkin(9) * t94 + t12;
t7 = (pkin(5) * t77 - pkin(9) * t117) * qJD(3) + t11;
t4 = -t88 * qJD(6) + t76 * t13 - t73 * t14;
t3 = -t89 * qJD(6) - t73 * t13 - t76 * t14;
t2 = -t90 * qJD(6) + t76 * t7 - t73 * t8;
t1 = -t91 * qJD(6) - t73 * t7 - t76 * t8;
t22 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t69 ^ 2 * t75 * t105 + 0.2e1 * t38 * t26 + 0.2e1 * t118, 0, 0, 0, 0.2e1 * t23 * t13 + 0.2e1 * t24 * t14 + 0.2e1 * t118, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t96, -t95, 0, 0, 0, 0, 0, -t29, t28, t79, t29, -t28 (t52 * t106 - t34 * t78) * t69 + t79 * pkin(8), t25 * t112 + t13 * t74 + (-t39 * t113 + t23 * t77) * qJD(3), -t25 * t116 - t14 * t74 + (t39 * t117 - t24 * t77) * qJD(3) (t13 * t68 - t14 * t70) * t77 + (-t23 * t68 + t24 * t70) * t102, t23 * t11 + t24 * t12 + t13 * t20 + t14 * t21 + t25 * t56 - t39 * t47, 0, 0, 0, 0, 0, -t39 * t17 + t25 * t31 + t4 * t74 + t89 * t63, t39 * t16 - t25 * t32 + t3 * t74 - t88 * t63; 0, 0, 0, 0, t58, 0.2e1 * (-t74 ^ 2 + t77 ^ 2) * qJD(3), 0, 0, 0, t74 * t100, t77 * t100, 0, -0.2e1 * t52 * t102 + 0.2e1 * t34 * t77, -0.2e1 * t34 * t74 - 0.2e1 * t52 * t63, 0.2e1 * t52 * t34, -0.2e1 * t47 * t112 + 0.2e1 * t11 * t74 + 0.2e1 * (-t56 * t113 + t20 * t77) * qJD(3), 0.2e1 * t47 * t116 - 0.2e1 * t12 * t74 + 0.2e1 * (t56 * t117 - t21 * t77) * qJD(3), 0.2e1 * (t11 * t68 - t12 * t70) * t77 + 0.2e1 * (-t20 * t68 + t21 * t70) * t102, 0.2e1 * t20 * t11 + 0.2e1 * t21 * t12 - 0.2e1 * t56 * t47, -0.2e1 * t32 * t16, -0.2e1 * t16 * t31 - 0.2e1 * t32 * t17, 0.2e1 * t16 * t74 - 0.2e1 * t32 * t63, 0.2e1 * t17 * t74 - 0.2e1 * t31 * t63, t58, -0.2e1 * t35 * t17 + 0.2e1 * t2 * t74 + 0.2e1 * t30 * t31 + 0.2e1 * t91 * t63, 0.2e1 * t1 * t74 + 0.2e1 * t35 * t16 - 0.2e1 * t30 * t32 - 0.2e1 * t90 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t25, 0, t26, t25, -t26 * pkin(3) + t83, t25 * t68, t25 * t70, -t6, t6 * t72 + (-t23 * t70 - t24 * t68) * qJD(5) + t83, 0, 0, 0, 0, 0, t25 * t43 + t39 * t37, t123 * t25 - t39 * t36; 0, 0, 0, 0, 0, 0, t63, -t102, 0, -t61, t98, t80, t61, -t98, t80 * pkin(8), -t122 * t70 - t47 * t68, t122 * t68 - t47 * t70, -t5, -t47 * qJ(4) + t56 * qJD(4) + t5 * t72 + (-t20 * t70 - t21 * t68) * qJD(5), t123 * t16 + t32 * t36, t123 * t17 - t16 * t43 + t36 * t31 + t32 * t37, t19, t81, 0, qJD(4) * t31 + t10 * t74 - t60 * t17 + t30 * t43 + t35 * t37 - t86 * t63, -qJD(4) * t32 + t123 * t30 + t60 * t16 - t35 * t36 - t87 * t63 + t9 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, 0.2e1 * t101, t68 * t121, t70 * t121, 0.2e1 * t51, -0.2e1 * t51 * t72 + 0.2e1 * t101, -0.2e1 * t123 * t36, -0.2e1 * t123 * t37 + 0.2e1 * t36 * t43, 0, 0, 0, 0.2e1 * qJD(4) * t43 + 0.2e1 * t60 * t37, 0.2e1 * qJD(4) * t123 - 0.2e1 * t60 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, 0, t61, t70 * t63, -t68 * t63, 0, t5, 0, 0, 0, 0, 0, t19, t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, t97, 0, -t47, 0, 0, 0, 0, 0, -t17, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), 0, 0, 0, 0, 0, t37, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t17, t63, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t37, 0, t10, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t22;
