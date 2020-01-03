% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPPR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:15
% EndTime: 2019-12-31 17:45:17
% DurationCPUTime: 0.42s
% Computational Cost: add. (647->82), mult. (1309->120), div. (0->0), fcn. (854->6), ass. (0->60)
t83 = 2 * qJD(3);
t42 = sin(pkin(7)) * pkin(1) + qJ(3);
t82 = qJD(1) * t42;
t54 = cos(qJ(5));
t51 = cos(pkin(8));
t68 = qJD(1) * t51;
t65 = t54 * t68;
t53 = sin(qJ(5));
t49 = sin(pkin(8));
t69 = qJD(1) * t49;
t66 = t53 * t69;
t27 = t65 - t66;
t81 = -t53 * t49 + t54 * t51;
t41 = -cos(pkin(7)) * pkin(1) - pkin(2) - qJ(4);
t80 = t41 * qJD(1);
t79 = t27 ^ 2;
t47 = qJD(3) * qJD(1);
t64 = 0.2e1 * t47;
t77 = -pkin(6) + t41;
t38 = qJD(5) * t66;
t16 = qJD(5) * t65 - t38;
t32 = t54 * t49 + t53 * t51;
t29 = t32 * qJD(5);
t71 = qJD(1) * t32;
t76 = -t16 * t81 + t29 * t71;
t75 = t27 * t71;
t31 = qJD(3) + t80;
t14 = t51 * qJD(2) + t49 * t31;
t72 = t49 ^ 2 + t51 ^ 2;
t30 = t81 * qJD(5);
t18 = t30 * qJD(5);
t67 = qJD(1) * qJD(4);
t34 = qJD(4) + t82;
t13 = -t49 * qJD(2) + t51 * t31;
t11 = -pkin(6) * t68 + t13;
t12 = -pkin(6) * t69 + t14;
t3 = t54 * t11 - t53 * t12;
t4 = t53 * t11 + t54 * t12;
t63 = t13 * t51 + t14 * t49;
t15 = qJD(1) * t29;
t62 = -t15 * t81 - t27 * t29;
t61 = -t32 * t15 + t27 * t30;
t60 = t32 * t16 + t30 * t71;
t23 = t77 * t49;
t24 = t77 * t51;
t8 = t54 * t23 + t53 * t24;
t7 = -t53 * t23 + t54 * t24;
t58 = t27 * qJD(4);
t57 = t32 * qJD(4);
t1 = -qJD(1) * t57 + t3 * qJD(5);
t2 = -t4 * qJD(5) - t58;
t56 = t1 * t32 + t2 * t81 - t3 * t29 + t4 * t30;
t55 = qJD(1) ^ 2;
t35 = t49 * pkin(4) + t42;
t22 = t71 ^ 2;
t21 = pkin(4) * t69 + t34;
t17 = t29 * qJD(5);
t6 = -qJD(4) * t81 - t8 * qJD(5);
t5 = t7 * qJD(5) - t57;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t82 * t83, 0, 0, 0, 0, 0, 0, t49 * t64, t51 * t64, 0.2e1 * t72 * t67, (t34 + t82) * qJD(3) + (-t72 * t80 - t63) * qJD(4), t62, -t61 + t76, -t17, t60, -t18, 0, t6 * qJD(5) + t35 * t16 + t21 * t30 + t71 * t83, -t5 * qJD(5) - t35 * t15 - t21 * t29 + (qJD(1) * t81 + t27) * qJD(3), t7 * t15 - t8 * t16 - t6 * t27 - t5 * t71 - t56, t1 * t8 + t2 * t7 + t3 * t6 + t4 * t5 + (qJD(1) * t35 + t21) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t17, t61 + t76, t1 * t81 - t2 * t32 - t4 * t29 - t3 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t82 * qJD(1), 0, 0, 0, 0, 0, 0, -t55 * t49, -t55 * t51, 0, (-t72 * qJD(4) - t34) * qJD(1), 0, 0, 0, 0, 0, 0, -qJD(1) * t71 - t17, -qJD(1) * t27 - t18, -t60 - t62, -t21 * qJD(1) + t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72 * t55, t63 * qJD(1) + t47, 0, 0, 0, 0, 0, 0, -t38 + (t27 + t65) * qJD(5), -0.2e1 * t71 * qJD(5), -t22 - t79, t3 * t27 + t4 * t71 + t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t22 + t79, 0, -t75, t38 + (t27 - t65) * qJD(5), 0, -t21 * t27 - t58, t21 * t71 + t32 * t67, 0, 0;];
tauc_reg = t9;
