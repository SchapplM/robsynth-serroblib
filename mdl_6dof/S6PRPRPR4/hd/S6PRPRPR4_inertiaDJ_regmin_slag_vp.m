% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRPR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:41:17
% EndTime: 2019-03-08 19:41:20
% DurationCPUTime: 1.13s
% Computational Cost: add. (1363->156), mult. (3648->306), div. (0->0), fcn. (3781->12), ass. (0->99)
t71 = sin(pkin(12));
t74 = cos(pkin(12));
t77 = sin(qJ(6));
t80 = cos(qJ(6));
t123 = -t77 * t71 + t80 * t74;
t54 = t80 * t71 + t77 * t74;
t47 = t54 * qJD(6);
t122 = t123 * qJD(6);
t121 = 0.2e1 * t122;
t75 = cos(pkin(11));
t66 = -t75 * pkin(3) - pkin(2);
t120 = 0.2e1 * t66;
t72 = sin(pkin(11));
t78 = sin(qJ(4));
t81 = cos(qJ(4));
t55 = t81 * t72 + t78 * t75;
t119 = t55 * t71;
t112 = t81 * t75;
t53 = t78 * t72 - t112;
t48 = t53 * qJD(4);
t118 = t71 * t48;
t73 = sin(pkin(6));
t79 = sin(qJ(2));
t117 = t73 * t79;
t82 = cos(qJ(2));
t116 = t73 * t82;
t115 = t74 * t48;
t111 = pkin(8) + qJ(3);
t110 = pkin(9) + qJ(5);
t49 = t55 * qJD(4);
t25 = t49 * pkin(4) + t48 * qJ(5) - t55 * qJD(5);
t105 = qJD(4) * t81;
t58 = t111 * t72;
t60 = t111 * t75;
t29 = t58 * t105 - qJD(3) * t112 + (qJD(3) * t72 + qJD(4) * t60) * t78;
t8 = t71 * t25 - t74 * t29;
t38 = t53 * pkin(4) - t55 * qJ(5) + t66;
t41 = -t78 * t58 + t81 * t60;
t16 = t71 * t38 + t74 * t41;
t109 = t71 ^ 2 + t74 ^ 2;
t108 = t72 ^ 2 + t75 ^ 2;
t107 = qJD(2) * t73;
t106 = qJD(2) * t79;
t102 = t73 * t106;
t101 = t82 * t107;
t100 = t108 * t82;
t7 = t74 * t25 + t71 * t29;
t15 = t74 * t38 - t71 * t41;
t40 = t81 * t58 + t78 * t60;
t99 = 0.2e1 * t109 * qJD(5);
t98 = 0.2e1 * t108 * qJD(3);
t97 = t7 * t74 + t8 * t71;
t96 = -t7 * t71 + t8 * t74;
t10 = -pkin(9) * t119 + t16;
t9 = -t74 * t55 * pkin(9) + t53 * pkin(5) + t15;
t95 = t80 * t10 + t77 * t9;
t94 = t77 * t10 - t80 * t9;
t76 = cos(pkin(6));
t44 = -t72 * t117 + t76 * t75;
t45 = t75 * t117 + t76 * t72;
t17 = -t44 * t105 - t101 * t112 + (qJD(4) * t45 + t72 * t101) * t78;
t13 = t74 * t102 + t71 * t17;
t14 = t71 * t102 - t74 * t17;
t93 = t13 * t74 + t14 * t71;
t92 = -t13 * t71 + t14 * t74;
t33 = t78 * t44 + t81 * t45;
t18 = t33 * qJD(4) + t55 * t101;
t32 = -t81 * t44 + t78 * t45;
t91 = t18 * t55 - t32 * t48;
t23 = -t74 * t116 - t71 * t33;
t24 = -t71 * t116 + t74 * t33;
t90 = t80 * t23 - t77 * t24;
t89 = t77 * t23 + t80 * t24;
t30 = t55 * qJD(3) + t41 * qJD(4);
t88 = t30 * t55 - t40 * t48;
t87 = -t44 * t72 + t45 * t75;
t86 = -t122 * t53 - t54 * t49;
t57 = t110 * t71;
t59 = t110 * t74;
t85 = -t80 * t57 - t77 * t59;
t84 = -t77 * t57 + t80 * t59;
t83 = pkin(4) * t48 - qJ(5) * t49 - qJD(5) * t53;
t65 = -t74 * pkin(5) - pkin(4);
t35 = t123 * t55;
t34 = t54 * t55;
t31 = pkin(5) * t119 + t40;
t28 = -t54 * qJD(5) - t84 * qJD(6);
t27 = -qJD(5) * t123 - t85 * qJD(6);
t20 = t123 * t49 - t47 * t53;
t19 = -pkin(5) * t118 + t30;
t12 = t122 * t55 - t54 * t48;
t11 = -t123 * t48 - t47 * t55;
t6 = pkin(9) * t118 + t8;
t5 = t49 * pkin(5) + pkin(9) * t115 + t7;
t4 = -t89 * qJD(6) + t80 * t13 - t77 * t14;
t3 = -t90 * qJD(6) - t77 * t13 - t80 * t14;
t2 = -t95 * qJD(6) + t80 * t5 - t77 * t6;
t1 = t94 * qJD(6) - t77 * t5 - t80 * t6;
t21 = [0, 0, 0, 0, 0, 0, 0, 0.2e1 * (t87 - t117) * t101, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t23 * t13 + 0.2e1 * t24 * t14 + 0.2e1 * t32 * t18, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t102, -t101, -t75 * t102, t72 * t102, t100 * t107, t87 * qJD(3) + (-pkin(2) * t79 + qJ(3) * t100) * t107, 0, 0, 0, 0, 0 (t53 * t106 - t49 * t82) * t73 (t55 * t106 + t48 * t82) * t73, t13 * t53 + t23 * t49 + t91 * t71, -t14 * t53 - t24 * t49 + t91 * t74, -t93 * t55 - (-t23 * t74 - t24 * t71) * t48, t13 * t15 + t14 * t16 + t18 * t40 + t23 * t7 + t24 * t8 + t32 * t30, 0, 0, 0, 0, 0, t32 * t12 + t18 * t34 + t4 * t53 + t90 * t49, t32 * t11 + t18 * t35 + t3 * t53 - t89 * t49; 0, 0, 0, 0, 0, 0, t98, qJ(3) * t98, -0.2e1 * t55 * t48, 0.2e1 * t48 * t53 - 0.2e1 * t55 * t49, 0, 0, 0, t49 * t120, -t48 * t120, 0.2e1 * t15 * t49 + 0.2e1 * t7 * t53 + 0.2e1 * t88 * t71, -0.2e1 * t16 * t49 - 0.2e1 * t8 * t53 + 0.2e1 * t88 * t74, -0.2e1 * t97 * t55 - 0.2e1 * (-t15 * t74 - t16 * t71) * t48, 0.2e1 * t15 * t7 + 0.2e1 * t16 * t8 + 0.2e1 * t40 * t30, 0.2e1 * t35 * t11, -0.2e1 * t11 * t34 - 0.2e1 * t35 * t12, 0.2e1 * t11 * t53 + 0.2e1 * t35 * t49, -0.2e1 * t12 * t53 - 0.2e1 * t34 * t49, 0.2e1 * t53 * t49, 0.2e1 * t31 * t12 + 0.2e1 * t19 * t34 + 0.2e1 * t2 * t53 - 0.2e1 * t94 * t49, 0.2e1 * t1 * t53 + 0.2e1 * t31 * t11 + 0.2e1 * t19 * t35 - 0.2e1 * t95 * t49; 0, 0, 0, 0, 0, 0, 0, t102, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t48, t74 * t49, -t71 * t49, t109 * t48, t97, 0, 0, 0, 0, 0, t20, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t17, -t18 * t74, t18 * t71, t92, -t18 * pkin(4) + (-t23 * t71 + t24 * t74) * qJD(5) + t92 * qJ(5), 0, 0, 0, 0, 0, -t123 * t18 + t32 * t47, t122 * t32 + t18 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t49, 0, -t30, t29, -t30 * t74 + t83 * t71, t30 * t71 + t83 * t74, t96, -t30 * pkin(4) + (-t15 * t71 + t16 * t74) * qJD(5) + t96 * qJ(5), t11 * t54 + t122 * t35, t11 * t123 - t54 * t12 - t122 * t34 - t35 * t47, -t86, t20, 0, t65 * t12 - t123 * t19 + t28 * t53 + t31 * t47 + t85 * t49, t65 * t11 + t122 * t31 + t19 * t54 + t27 * t53 - t84 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, qJ(5) * t99, t54 * t121, 0.2e1 * t122 * t123 - 0.2e1 * t54 * t47, 0, 0, 0, 0.2e1 * t65 * t47, t65 * t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t118, -t115, 0, t30, 0, 0, 0, 0, 0, t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t12, t49, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, -t47, 0, t28, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t21;
