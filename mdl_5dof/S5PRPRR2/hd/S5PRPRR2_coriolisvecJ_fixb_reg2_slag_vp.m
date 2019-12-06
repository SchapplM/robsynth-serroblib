% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:45:13
% EndTime: 2019-12-05 15:45:16
% DurationCPUTime: 0.51s
% Computational Cost: add. (1289->100), mult. (2893->152), div. (0->0), fcn. (2192->8), ass. (0->80)
t53 = sin(pkin(9));
t103 = pkin(2) * t53;
t54 = cos(pkin(9));
t47 = t54 * pkin(2) + pkin(3);
t56 = sin(qJ(4));
t59 = cos(qJ(4));
t88 = t59 * t103 + t56 * t47;
t34 = pkin(7) + t88;
t61 = qJD(5) ^ 2;
t50 = qJD(2) + qJD(4);
t57 = sin(qJ(2));
t60 = cos(qJ(2));
t40 = t53 * t60 + t54 * t57;
t36 = t40 * qJD(1);
t72 = t53 * t57 - t54 * t60;
t38 = t72 * qJD(1);
t89 = -t88 * qJD(4) + t59 * t36 - t56 * t38;
t81 = t89 * t50;
t106 = t34 * t61 - t81;
t45 = qJD(2) * pkin(2) + t60 * qJD(1);
t85 = qJD(1) * t57;
t24 = t54 * t45 - t53 * t85;
t22 = qJD(2) * pkin(3) + t24;
t25 = t53 * t45 + t54 * t85;
t16 = t56 * t22 + t59 * t25;
t14 = t50 * pkin(7) + t16;
t55 = sin(qJ(5));
t58 = cos(qJ(5));
t10 = t55 * qJD(3) + t58 * t14;
t35 = t40 * qJD(2);
t31 = qJD(1) * t35;
t37 = t72 * qJD(2);
t32 = qJD(1) * t37;
t92 = t56 * t25;
t5 = (qJD(4) * t22 - t32) * t59 - qJD(4) * t92 - t56 * t31;
t9 = t58 * qJD(3) - t55 * t14;
t2 = t9 * qJD(5) + t58 * t5;
t3 = -t10 * qJD(5) - t55 * t5;
t105 = t2 * t58 - t3 * t55 + (-t10 * t55 - t58 * t9) * qJD(5);
t69 = -t56 * t103 + t59 * t47;
t90 = t69 * qJD(4) + t56 * t36 + t59 * t38;
t104 = t90 * t50;
t75 = t10 * t58 - t55 * t9;
t6 = t16 * qJD(4) + t59 * t31 - t56 * t32;
t102 = t50 * pkin(4);
t73 = -t56 * t40 - t59 * t72;
t100 = t6 * t73;
t7 = t73 * qJD(4) - t56 * t35 - t59 * t37;
t99 = t7 * t50;
t20 = t59 * t40 - t56 * t72;
t8 = t20 * qJD(4) + t59 * t35 - t56 * t37;
t98 = t8 * t50;
t15 = t59 * t22 - t92;
t13 = -t15 - t102;
t84 = qJD(5) * t58;
t97 = t13 * t84 + t6 * t55;
t95 = t15 * t50;
t94 = t16 * t50;
t91 = t61 * t55;
t51 = t55 ^ 2;
t52 = t58 ^ 2;
t87 = t51 - t52;
t86 = t51 + t52;
t49 = t50 ^ 2;
t83 = t55 * t49 * t58;
t82 = -t13 * t50 - t5;
t77 = t55 * t50 * t84;
t76 = pkin(7) * t61 - t94;
t74 = t20 * t61 + t98;
t71 = qJD(5) * (t15 - t102);
t70 = qJD(5) * (-t50 * t73 - t7);
t33 = -pkin(4) - t69;
t68 = qJD(5) * (t33 * t50 - t90);
t62 = qJD(2) ^ 2;
t48 = t61 * t58;
t43 = -0.2e1 * t77;
t42 = 0.2e1 * t77;
t28 = -0.2e1 * t87 * t50 * qJD(5);
t11 = t13 * qJD(5) * t55;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62 * t57, -t62 * t60, 0, 0, 0, 0, 0, 0, 0, 0, -t35 * qJD(2), t37 * qJD(2), 0, -t24 * t35 - t25 * t37 + t31 * t72 - t32 * t40, 0, 0, 0, 0, 0, 0, -t98, -t99, 0, -t15 * t8 + t16 * t7 + t5 * t20 - t100, 0, 0, 0, 0, 0, 0, t55 * t70 - t74 * t58, t74 * t55 + t58 * t70, t86 * t99, t105 * t20 + t13 * t8 + t75 * t7 - t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 * t36 + t25 * t38 + (-t31 * t54 - t32 * t53) * pkin(2), 0, 0, 0, 0, 0, 0, t81 - t6, -t5 - t104, 0, t89 * t15 + t90 * t16 + t5 * t88 - t6 * t69, t42, t28, t48, t43, -t91, 0, t11 + t55 * t68 + (-t6 - t106) * t58, t106 * t55 + t58 * t68 + t97, t86 * t104 + t105, t105 * t34 - t89 * t13 + t6 * t33 + t75 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, -t48, 0, t75 * qJD(5) + t2 * t55 + t3 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6 + t94, -t5 + t95, 0, 0, t42, t28, t48, t43, -t91, 0, t11 + t55 * t71 + (-t6 - t76) * t58, t76 * t55 + t58 * t71 + t97, -t86 * t95 + t105, -t6 * pkin(4) + pkin(7) * t105 - t13 * t16 - t75 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, t87 * t49, 0, t83, 0, 0, t82 * t55, t82 * t58, 0, 0;];
tauc_reg = t1;
