% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRRP8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:16:24
% EndTime: 2019-03-09 02:16:26
% DurationCPUTime: 0.91s
% Computational Cost: add. (1536->138), mult. (3084->252), div. (0->0), fcn. (2921->6), ass. (0->88)
t56 = sin(pkin(9));
t57 = cos(pkin(9));
t109 = sin(qJ(4));
t79 = qJD(4) * t109;
t110 = cos(qJ(4));
t80 = qJD(4) * t110;
t33 = -t56 * t79 + t57 * t80;
t38 = t109 * t57 + t110 * t56;
t105 = t38 * t33;
t32 = -t56 * t80 - t57 * t79;
t66 = -t109 * t56 + t110 * t57;
t107 = t66 * t32;
t35 = t38 ^ 2;
t36 = t66 ^ 2;
t60 = cos(qJ(5));
t59 = sin(qJ(5));
t96 = qJD(5) * t59;
t120 = (t35 + t36) * t96 + 0.2e1 * (-t105 - t107) * t60;
t47 = t56 * pkin(3) + qJ(2);
t22 = t38 * pkin(4) - pkin(8) * t66 + t47;
t58 = -pkin(1) - qJ(3);
t111 = -pkin(7) + t58;
t40 = t111 * t56;
t41 = t111 * t57;
t28 = t109 * t41 + t110 * t40;
t118 = t59 * t22 + t60 * t28;
t54 = t59 ^ 2;
t55 = t60 ^ 2;
t99 = t54 - t55;
t78 = qJD(5) * t99;
t42 = (t56 ^ 2 + t57 ^ 2) * qJD(3);
t10 = t38 * qJD(3) + t40 * t79 - t41 * t80;
t21 = t33 * pkin(4) - t32 * pkin(8) + qJD(2);
t4 = -qJD(5) * t118 + t59 * t10 + t60 * t21;
t116 = 0.2e1 * qJD(2);
t115 = 2 * qJD(6);
t114 = pkin(8) * t33;
t113 = pkin(8) * t38;
t112 = t33 * pkin(5);
t11 = t66 * qJD(3) + t28 * qJD(4);
t108 = t11 * t59;
t106 = t66 * t60;
t104 = t59 * t32;
t103 = t59 * t33;
t102 = t60 * t32;
t98 = t54 + t55;
t97 = t33 * qJ(6);
t50 = qJD(5) * t60;
t95 = t38 * qJD(6);
t94 = t59 * qJD(6);
t93 = qJ(2) * qJD(2);
t92 = -0.2e1 * t103;
t90 = -0.2e1 * pkin(4) * qJD(5);
t88 = -0.2e1 * t104 * t66 - t36 * t50;
t87 = pkin(8) * t96;
t86 = pkin(8) * t50;
t85 = t38 * t50;
t84 = t59 * t50;
t83 = t98 * t33;
t82 = -0.4e1 * t59 * t106;
t77 = -pkin(4) * t32 - t114;
t76 = -pkin(4) * t66 - t113;
t6 = t38 * qJ(6) + t118;
t71 = t60 * t22 - t59 * t28;
t7 = -t38 * pkin(5) - t71;
t75 = t59 * t7 + t6 * t60;
t74 = -t59 * t6 + t60 * t7;
t27 = t109 * t40 - t110 * t41;
t73 = -t60 * pkin(5) - t59 * qJ(6);
t72 = pkin(5) * t59 - qJ(6) * t60;
t30 = -pkin(5) * t96 + qJ(6) * t50 + t94;
t43 = -pkin(4) + t73;
t69 = t30 * t66 - t32 * t43;
t19 = -t50 * t66 - t104;
t67 = -t66 * t96 + t102;
t18 = t85 + t103;
t16 = -t60 * t33 + t38 * t96;
t3 = t60 * t10 - t59 * t21 - t22 * t50 + t28 * t96;
t63 = t73 * qJD(5) + t60 * qJD(6);
t5 = t72 * t32 - t63 * t66 + t11;
t65 = -t5 + (t43 * t66 - t113) * qJD(5);
t8 = t66 * t72 + t27;
t64 = -qJD(5) * t8 + t114 + t69;
t1 = -t3 + t95 + t97;
t2 = -t112 - t4;
t62 = t75 * qJD(5) + t1 * t59 - t2 * t60;
t61 = t74 * qJD(5) + t1 * t60 + t2 * t59;
t9 = [0, 0, 0, 0, t116, 0.2e1 * t93, t56 * t116, t57 * t116, 0.2e1 * t42, -0.2e1 * t42 * t58 + 0.2e1 * t93, 0.2e1 * t107, -0.2e1 * t32 * t38 - 0.2e1 * t33 * t66, 0, 0, 0, 0.2e1 * qJD(2) * t38 + 0.2e1 * t47 * t33, 0.2e1 * qJD(2) * t66 + 0.2e1 * t47 * t32, 0.2e1 * t55 * t107 - 0.2e1 * t36 * t84, t32 * t82 + 0.2e1 * t36 * t78, 0.2e1 * t38 * t102 - 0.2e1 * t16 * t66, -0.2e1 * t38 * t104 - 0.2e1 * t18 * t66, 0.2e1 * t105, 0.2e1 * t108 * t66 - 0.2e1 * t19 * t27 + 0.2e1 * t71 * t33 + 0.2e1 * t4 * t38, 0.2e1 * t11 * t106 - 0.2e1 * t118 * t33 + 0.2e1 * t67 * t27 + 0.2e1 * t3 * t38, 0.2e1 * t8 * t104 - 0.2e1 * t2 * t38 - 0.2e1 * t7 * t33 - 0.2e1 * (-t5 * t59 - t8 * t50) * t66, 0.2e1 * t74 * t32 - 0.2e1 * t62 * t66, -0.2e1 * t8 * t102 + 0.2e1 * t1 * t38 + 0.2e1 * t6 * t33 - 0.2e1 * (t5 * t60 - t8 * t96) * t66, 0.2e1 * t6 * t1 + 0.2e1 * t7 * t2 + 0.2e1 * t8 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35 * t50 + t38 * t92 + t88, t120 (-t85 + t92) * t38 + t88, 0, -t120, -t8 * t32 + t75 * t33 + t61 * t38 - t5 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t38 * t83 + 0.2e1 * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, t33, t32, 0, 0, 0, 0, 0, -t16, -t18, -t16, -t98 * t32, t18, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t33, 0, -t11, t10, t59 * t102 - t66 * t78, qJD(5) * t82 - t99 * t32, t18, -t16, 0, -t11 * t60 + t77 * t59 + (t27 * t59 + t76 * t60) * qJD(5), t108 + t77 * t60 + (t27 * t60 - t76 * t59) * qJD(5), -t64 * t59 + t65 * t60, t61, t65 * t59 + t64 * t60, t61 * pkin(8) - t8 * t30 + t5 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t33, 0, 0, 0, 0, 0, t67, t19, t67, t83, -t19, pkin(8) * t83 + t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t84, -0.2e1 * t78, 0, 0, 0, t59 * t90, t60 * t90, 0.2e1 * t30 * t60 + 0.2e1 * t43 * t96, 0, 0.2e1 * t30 * t59 - 0.2e1 * t43 * t50, -0.2e1 * t43 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t19, t33, t4, t3, t4 + 0.2e1 * t112, t73 * t32 - (-t72 * qJD(5) + t94) * t66, -t3 + 0.2e1 * t95 + 0.2e1 * t97, -t2 * pkin(5) + t1 * qJ(6) + t6 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t16, -t18, 0, -t16, -t72 * t33 + t63 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t50, -t96, 0, t50, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t96, 0, -t86, t87, -t86, t63, -t87, t63 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, qJ(6) * t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t67, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
