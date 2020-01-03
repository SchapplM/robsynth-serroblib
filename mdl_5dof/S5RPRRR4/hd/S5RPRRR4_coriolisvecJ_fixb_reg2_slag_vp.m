% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:52:21
% EndTime: 2020-01-03 11:52:26
% DurationCPUTime: 0.82s
% Computational Cost: add. (1812->113), mult. (4046->156), div. (0->0), fcn. (2322->8), ass. (0->81)
t43 = cos(pkin(9)) * pkin(1) + pkin(2);
t41 = t43 * qJD(1);
t56 = sin(qJ(3));
t59 = cos(qJ(3));
t102 = sin(pkin(9)) * pkin(1);
t83 = qJD(1) * t102;
t29 = t59 * t41 - t56 * t83;
t49 = qJD(1) + qJD(3);
t22 = t49 * pkin(3) + t29;
t55 = sin(qJ(4));
t30 = t56 * t41 + t59 * t83;
t58 = cos(qJ(4));
t93 = t58 * t30;
t16 = t55 * t22 + t93;
t47 = qJD(4) + t49;
t14 = t47 * pkin(8) + t16;
t54 = sin(qJ(5));
t57 = cos(qJ(5));
t10 = t54 * qJD(2) + t57 * t14;
t25 = t29 * qJD(3);
t26 = t30 * qJD(3);
t87 = qJD(4) * t55;
t78 = -t55 * t26 - t30 * t87;
t5 = (qJD(4) * t22 + t25) * t58 + t78;
t9 = t57 * qJD(2) - t54 * t14;
t2 = t9 * qJD(5) + t57 * t5;
t3 = -t10 * qJD(5) - t54 * t5;
t105 = t2 * t57 - t3 * t54 + (-t10 * t54 - t57 * t9) * qJD(5);
t94 = t55 * t30;
t18 = t58 * t29 - t94;
t88 = pkin(3) * qJD(4);
t104 = t58 * t88 - t18;
t74 = -t56 * t102 + t59 * t43;
t35 = pkin(3) + t74;
t36 = t59 * t102 + t56 * t43;
t91 = t55 * t35 + t58 * t36;
t17 = t55 * t29 + t93;
t44 = t55 * pkin(3) + pkin(8);
t60 = qJD(5) ^ 2;
t103 = (pkin(3) * t87 - t17) * t47 + t44 * t60;
t79 = t55 * t25 + t58 * t26;
t6 = t16 * qJD(4) + t79;
t101 = t47 * pkin(4);
t33 = t74 * qJD(3);
t34 = t36 * qJD(3);
t69 = t58 * t35 - t55 * t36;
t7 = qJD(4) * t69 + t58 * t33 - t55 * t34;
t100 = t7 * t47;
t8 = t91 * qJD(4) + t55 * t33 + t58 * t34;
t99 = t8 * t47;
t15 = t58 * t22 - t94;
t13 = -t15 - t101;
t86 = qJD(5) * t57;
t98 = t13 * t86 + t6 * t54;
t97 = t15 * t47;
t96 = t16 * t47;
t92 = t60 * t54;
t50 = t54 ^ 2;
t51 = t57 ^ 2;
t90 = t50 - t51;
t89 = t50 + t51;
t46 = t47 ^ 2;
t84 = t54 * t46 * t57;
t81 = -pkin(3) * t47 - t22;
t80 = -t13 * t47 - t5;
t75 = t54 * t47 * t86;
t72 = pkin(8) * t60 - t96;
t71 = t10 * t57 - t54 * t9;
t20 = pkin(8) + t91;
t70 = t20 * t60 + t99;
t68 = qJD(5) * (t15 - t101);
t19 = -pkin(4) - t69;
t67 = qJD(5) * (t19 * t47 - t7);
t45 = -t58 * pkin(3) - pkin(4);
t66 = qJD(5) * (t45 * t47 - t104);
t48 = t60 * t57;
t38 = -0.2e1 * t75;
t37 = 0.2e1 * t75;
t28 = -0.2e1 * t90 * t47 * qJD(5);
t11 = t13 * qJD(5) * t54;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34 * t49 - t26, -t33 * t49 - t25, 0, t25 * t36 - t26 * t74 - t29 * t34 + t30 * t33, 0, 0, 0, 0, 0, 0, -t6 - t99, -t5 - t100, 0, -t15 * t8 + t16 * t7 + t5 * t91 - t6 * t69, t37, t28, t48, t38, -t92, 0, t11 + t54 * t67 + (-t6 - t70) * t57, t70 * t54 + t57 * t67 + t98, t89 * t100 + t105, t105 * t20 + t13 * t8 + t6 * t19 + t71 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, -t48, 0, t71 * qJD(5) + t2 * t54 + t3 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 * t49 - t26, t29 * t49 - t25, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t47 + (t81 * t55 - t93) * qJD(4) - t79, t18 * t47 + (t81 * qJD(4) - t25) * t58 - t78, 0, t15 * t17 - t16 * t18 + (t5 * t55 - t58 * t6 + (-t15 * t55 + t16 * t58) * qJD(4)) * pkin(3), t37, t28, t48, t38, -t92, 0, t11 + t54 * t66 + (-t103 - t6) * t57, t103 * t54 + t57 * t66 + t98, t104 * t47 * t89 + t105, -t13 * t17 + t6 * t45 - t71 * t18 + (t13 * t55 + t71 * t58) * t88 + t105 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6 + t96, -t5 + t97, 0, 0, t37, t28, t48, t38, -t92, 0, t11 + t54 * t68 + (-t6 - t72) * t57, t72 * t54 + t57 * t68 + t98, -t89 * t97 + t105, -t6 * pkin(4) + pkin(8) * t105 - t13 * t16 - t71 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, t90 * t46, 0, t84, 0, 0, t80 * t54, t80 * t57, 0, 0;];
tauc_reg = t1;
