% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% tauc_reg [5x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:14
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:14:04
% EndTime: 2021-01-15 11:14:06
% DurationCPUTime: 0.52s
% Computational Cost: add. (1013->146), mult. (2433->198), div. (0->0), fcn. (1522->6), ass. (0->88)
t60 = sin(pkin(7)) * pkin(1) + pkin(6);
t103 = qJ(4) + t60;
t66 = sin(pkin(8));
t68 = cos(pkin(8));
t70 = sin(qJ(3));
t71 = cos(qJ(3));
t50 = t66 * t71 + t68 * t70;
t102 = qJD(1) * t50;
t40 = t102 ^ 2;
t100 = qJD(1) * t70;
t107 = t68 * t71;
t89 = qJD(1) * t107;
t42 = t66 * t100 - t89;
t113 = -t42 ^ 2 - t40;
t83 = t103 * qJD(1);
t33 = t71 * qJD(2) - t83 * t70;
t34 = t70 * qJD(2) + t83 * t71;
t112 = pkin(3) * t70;
t92 = qJD(1) * qJD(4);
t25 = t33 * qJD(3) + t71 * t92;
t74 = -t34 * qJD(3) - t70 * t92;
t2 = t66 * t25 - t68 * t74;
t48 = t103 * t71;
t86 = t103 * t70;
t22 = t66 * t48 + t68 * t86;
t111 = t2 * t22;
t49 = t66 * t70 - t107;
t110 = t2 * t49;
t62 = -cos(pkin(7)) * pkin(1) - pkin(2);
t51 = -t71 * pkin(3) + t62;
t101 = qJD(1) * t51;
t41 = qJD(4) + t101;
t16 = t42 * pkin(4) - qJ(5) * t102 + t41;
t109 = t16 * t102;
t108 = t66 * t34;
t29 = t68 * t34;
t72 = qJD(3) ^ 2;
t106 = t72 * t70;
t105 = t72 * t71;
t3 = t68 * t25 + t66 * t74;
t31 = qJD(3) * pkin(3) + t33;
t9 = t66 * t31 + t29;
t104 = t70 ^ 2 - t71 ^ 2;
t53 = qJD(1) * t62;
t44 = t50 * qJD(3);
t99 = t44 * qJD(3);
t95 = t70 * qJD(3);
t47 = qJD(3) * t107 - t66 * t95;
t98 = t47 * qJD(3);
t97 = t53 * qJD(1);
t13 = t68 * t33 - t108;
t94 = qJD(5) - t13;
t93 = qJD(1) * qJD(3);
t91 = pkin(3) * t95;
t90 = pkin(3) * t100;
t88 = t70 * t93;
t87 = t71 * t93;
t85 = 0.2e1 * t102;
t84 = qJD(3) * t103;
t8 = t68 * t31 - t108;
t81 = 0.2e1 * qJD(3) * t53;
t36 = qJD(1) * t44;
t54 = t66 * t88;
t37 = t68 * t87 - t54;
t57 = pkin(3) * t88;
t80 = t36 * pkin(4) - t37 * qJ(5) + t57;
t12 = t66 * t33 + t29;
t79 = t12 * qJD(3) - t2;
t78 = t102 * t44 - t50 * t36 + t49 * t37 - t47 * t42;
t77 = -t70 * qJD(4) - t71 * t84;
t35 = t71 * qJD(4) - t70 * t84;
t14 = t66 * t35 - t68 * t77;
t15 = t68 * t35 + t66 * t77;
t23 = t68 * t48 - t66 * t86;
t76 = t102 * t14 - t15 * t42 + t2 * t50 + t22 * t37 - t23 * t36;
t75 = t85 * qJD(3);
t73 = qJD(1) ^ 2;
t61 = -t68 * pkin(3) - pkin(4);
t58 = t66 * pkin(3) + qJ(5);
t24 = -t54 + (-t42 + t89) * qJD(3);
t18 = t49 * pkin(4) - t50 * qJ(5) + t51;
t17 = pkin(4) * t102 + t42 * qJ(5) + t90;
t10 = t44 * pkin(4) - t47 * qJ(5) - t50 * qJD(5) + t91;
t7 = qJD(3) * qJ(5) + t9;
t6 = -qJD(3) * pkin(4) + qJD(5) - t8;
t5 = -qJD(5) * t102 + t80;
t1 = qJD(3) * qJD(5) + t3;
t4 = [0, 0, 0, 0, 0.2e1 * t70 * t87, -0.2e1 * t104 * t93, t105, -t106, 0, -t60 * t105 + t70 * t81, t60 * t106 + t71 * t81, t51 * t36 + t41 * t44 + (-t14 + (qJD(1) * t49 + t42) * t112) * qJD(3), t51 * t37 + t41 * t47 + (t85 * t112 - t15) * qJD(3), -t3 * t49 - t9 * t44 - t8 * t47 + t76, -t8 * t14 + t9 * t15 + t111 + t3 * t23 + (t41 + t101) * t91, -t14 * qJD(3) + t10 * t42 + t16 * t44 + t18 * t36 + t5 * t49, -t1 * t49 - t7 * t44 + t6 * t47 + t76, t15 * qJD(3) - t10 * t102 - t16 * t47 - t18 * t37 - t5 * t50, t1 * t23 + t16 * t10 + t6 * t14 + t7 * t15 + t5 * t18 + t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t105, -t99, -t98, t78, t3 * t50 - t8 * t44 + t9 * t47 + t110, -t99, t78, t98, t1 * t50 + t6 * t44 + t7 * t47 + t110; 0, 0, 0, 0, -t70 * t73 * t71, t104 * t73, 0, 0, 0, -t70 * t97, -t71 * t97, -t102 * t41 - t42 * t90 + t79, t13 * qJD(3) - t102 * t90 + t41 * t42 - t3, (-t12 + t9) * t102 + (t13 - t8) * t42 + (-t36 * t66 - t37 * t68) * pkin(3), t8 * t12 - t9 * t13 + (-t41 * t100 - t2 * t68 + t3 * t66) * pkin(3), -t17 * t42 - t109 + t79, -t58 * t36 + t61 * t37 + (-t12 + t7) * t102 + (t6 - t94) * t42, -t16 * t42 + t17 * t102 + (0.2e1 * qJD(5) - t13) * qJD(3) + t3, t1 * t58 - t6 * t12 - t16 * t17 + t2 * t61 + t94 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t24, t113, t102 * t8 + t9 * t42 + t57, t75, t113, -t24, t7 * t42 + (-qJD(5) - t6) * t102 + t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102 * t42, -t54 + (t42 + t89) * qJD(3), -t40 - t72, -t7 * qJD(3) + t109 + t2;];
tauc_reg = t4;
