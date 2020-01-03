% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRP4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:20
% EndTime: 2019-12-31 17:52:21
% DurationCPUTime: 0.49s
% Computational Cost: add. (681->116), mult. (1339->167), div. (0->0), fcn. (664->4), ass. (0->90)
t51 = sin(pkin(7));
t57 = qJD(4) ^ 2;
t58 = qJD(1) ^ 2;
t110 = t51 * (t57 + t58);
t52 = cos(pkin(7));
t56 = -pkin(1) - pkin(2);
t34 = t56 * qJD(1) + qJD(2);
t87 = qJD(1) * qJ(2);
t21 = t51 * t34 + t52 * t87;
t17 = -qJD(1) * pkin(6) + t21;
t54 = sin(qJ(4));
t55 = cos(qJ(4));
t10 = t55 * qJD(3) - t54 * t17;
t103 = t55 * t17;
t11 = t54 * qJD(3) + t103;
t63 = t10 * t54 - t11 * t55;
t109 = t63 * t52;
t86 = qJD(1) * qJD(2);
t79 = t52 * t86;
t4 = -t11 * qJD(4) - t54 * t79;
t84 = qJD(3) * qJD(4);
t89 = qJD(4) * t54;
t83 = t17 * t89 + (-t79 - t84) * t55;
t59 = -(t10 * t55 + t11 * t54) * qJD(4) - t83 * t55 - t4 * t54;
t80 = qJ(5) * t89;
t1 = (-t55 * qJD(5) + t80) * qJD(1) - t83;
t85 = qJD(1) * qJD(4);
t77 = t55 * t85;
t92 = qJD(1) * t54;
t100 = qJ(5) * t77 + qJD(5) * t92;
t2 = t4 + t100;
t88 = qJ(5) * qJD(1);
t8 = t54 * t88 + t10;
t95 = qJD(4) * pkin(4);
t5 = t8 + t95;
t9 = -t55 * t88 + t11;
t108 = -t1 * t55 + t2 * t54 + (t5 * t55 + t54 * t9) * qJD(4);
t37 = t51 * t86;
t107 = 0.2e1 * t37;
t106 = t5 - t8;
t105 = t55 * t9;
t104 = t52 * t58;
t102 = t55 * t58;
t101 = t57 * t55;
t99 = t52 * qJ(2) + t51 * t56;
t49 = t54 ^ 2;
t50 = t55 ^ 2;
t98 = t49 - t50;
t97 = t49 + t50;
t27 = -pkin(6) + t99;
t94 = qJ(5) - t27;
t76 = -t51 * qJ(2) + t52 * t56;
t26 = pkin(3) - t76;
t93 = qJD(1) * t26;
t91 = qJD(1) * t55;
t90 = qJD(2) * t52;
t82 = 0.2e1 * t86;
t81 = 0.2e1 * t85;
t78 = t54 * t85;
t20 = t52 * t34 - t51 * t87;
t16 = qJD(1) * pkin(3) - t20;
t12 = pkin(4) * t91 + qJD(5) + t16;
t22 = t55 * pkin(4) + t26;
t75 = -qJD(1) * t22 - t12;
t25 = -pkin(4) * t78 + t37;
t29 = -pkin(4) * t89 + t51 * qJD(2);
t74 = qJD(1) * t29 + t25;
t73 = -t16 - t90;
t72 = qJD(4) * t94;
t71 = t52 * t81;
t70 = -qJD(5) + t90;
t69 = t54 * t77;
t65 = t5 * t54 - t105;
t62 = t20 * t51 - t21 * t52;
t61 = -t27 * t57 + t107;
t60 = qJD(4) * (t73 - t93);
t48 = t57 * t54;
t36 = t54 * t102;
t32 = -0.2e1 * t69;
t31 = 0.2e1 * t69;
t30 = t98 * t58;
t24 = t98 * t81;
t23 = t97 * t104;
t19 = t94 * t55;
t18 = t94 * t54;
t15 = t54 * t110 + t55 * t71;
t14 = -t55 * t110 + t54 * t71;
t7 = -t70 * t54 + t55 * t72;
t6 = t54 * t72 + t70 * t55;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, qJ(2) * t82, 0, 0, 0, 0, 0, 0, t107, 0.2e1 * t79, 0, ((-t51 * t76 + t52 * t99) * qJD(1) - t62) * qJD(2), t31, -t24, -t101, t32, t48, 0, t54 * t60 + t61 * t55, -t61 * t54 + t55 * t60, -t97 * t79 - t59, t59 * t27 + (-t109 + (t16 + t93) * t51) * qJD(2), t31, -t24, -t101, t32, t48, 0, t74 * t55 + (t75 * t54 + t7) * qJD(4), -t74 * t54 + (t75 * t55 - t6) * qJD(4), (t54 * t7 - t55 * t6 + (t18 * t55 - t19 * t54) * qJD(4)) * qJD(1) + t108, -t1 * t19 + t12 * t29 + t2 * t18 + t25 * t22 + t5 * t7 + t9 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t58 * qJ(2), 0, 0, 0, 0, 0, 0, -t51 * t58, -t104, 0, t62 * qJD(1), 0, 0, 0, 0, 0, 0, t14, t15, t23, qJD(1) * t109 + (qJD(1) * t73 + t59) * t51, 0, 0, 0, 0, 0, 0, t14, t15, t23, (qJD(1) * t65 - t25) * t52 + (-qJD(1) * t12 - t108) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t101, 0, -t63 * qJD(4) + t4 * t55 - t54 * t83, 0, 0, 0, 0, 0, 0, -t48, -t101, 0, -qJD(4) * t65 + t1 * t54 + t2 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t30, 0, t36, 0, 0, (t16 - t90) * t92, t10 * qJD(4) + t16 * t91 + t83, 0, 0, -t36, t30, 0, t36, 0, 0, (t9 - t103) * qJD(4) + (pkin(4) * t102 - t84 + (t12 - t90) * qJD(1)) * t54 + t100, -t49 * t58 * pkin(4) + t8 * qJD(4) + (-t80 + (qJD(5) + t12) * t55) * qJD(1) + t83, (t95 - t106) * t91, t106 * t9 + (t12 * t92 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t78, -0.2e1 * t77, -t97 * t58, t37 + (t105 + (-t5 - t95) * t54) * qJD(1);];
tauc_reg = t3;
