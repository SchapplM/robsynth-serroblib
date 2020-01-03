% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPPR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:44:01
% EndTime: 2019-12-31 17:44:03
% DurationCPUTime: 0.41s
% Computational Cost: add. (555->92), mult. (1319->137), div. (0->0), fcn. (892->6), ass. (0->65)
t47 = sin(pkin(8));
t45 = t47 ^ 2;
t49 = cos(pkin(8));
t75 = t49 ^ 2 + t45;
t85 = t75 * qJD(1) * qJD(3);
t52 = cos(qJ(5));
t72 = qJD(1) * t47;
t63 = t52 * t72;
t51 = sin(qJ(5));
t71 = qJD(1) * t49;
t64 = t51 * t71;
t28 = t63 - t64;
t34 = t47 * t52 - t49 * t51;
t58 = cos(pkin(7)) * pkin(1) + t47 * qJ(4) + pkin(2);
t84 = (t49 * pkin(3) + t58) * qJD(1);
t33 = t47 * t51 + t49 * t52;
t73 = qJD(1) * t33;
t83 = t73 ^ 2;
t82 = t28 ^ 2;
t43 = sin(pkin(7)) * pkin(1) + qJ(3);
t81 = -pkin(6) + t43;
t40 = qJD(5) * t64;
t21 = qJD(5) * t63 - t40;
t24 = t33 * qJD(5);
t80 = -t34 * t21 + t24 * t73;
t38 = t43 * qJD(1);
t19 = t47 * qJD(2) + t49 * t38;
t79 = t19 * t49;
t78 = t28 * t73;
t17 = (pkin(3) + pkin(4)) * t49 + t58;
t74 = qJD(1) * t17;
t70 = qJD(3) * t47;
t69 = qJD(4) * t47;
t68 = t24 * qJD(5);
t66 = qJD(1) * qJD(4);
t65 = qJD(3) * t79 + t43 * t85;
t61 = t47 * t66;
t60 = 0.2e1 * t73;
t18 = t49 * qJD(2) - t47 * t38;
t16 = qJD(4) - t18;
t12 = -pkin(6) * t72 + t16;
t13 = -pkin(6) * t71 + t19;
t3 = t52 * t12 - t51 * t13;
t4 = t51 * t12 + t52 * t13;
t20 = qJD(1) * t24;
t25 = t34 * qJD(5);
t59 = -t33 * t20 + t28 * t25;
t29 = t81 * t47;
t30 = t81 * t49;
t7 = t52 * t29 - t51 * t30;
t8 = t51 * t29 + t52 * t30;
t57 = t28 * qJD(3);
t56 = t33 * qJD(3);
t54 = qJD(1) ^ 2;
t53 = qJD(5) ^ 2;
t37 = t75 * t54;
t23 = 0.2e1 * t85;
t22 = t25 * qJD(5);
t15 = qJD(3) - t84;
t11 = -qJD(3) + t74;
t6 = t34 * qJD(3) - t8 * qJD(5);
t5 = t7 * qJD(5) + t56;
t2 = -t4 * qJD(5) + t57;
t1 = qJD(1) * t56 + t3 * qJD(5);
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t18 * t70 + t65, 0, 0, 0, 0, 0, 0, 0.2e1 * t49 * t61, t23, 0.2e1 * t45 * t66, t16 * t70 + (-t15 + t84) * t69 + t65, -t20 * t34 - t28 * t24, -t59 + t80, -t68, t21 * t33 + t25 * t73, -t22, 0, t6 * qJD(5) + t11 * t25 + t17 * t21 + t60 * t69, -t5 * qJD(5) - t11 * t24 - t17 * t20 + (qJD(1) * t34 + t28) * t69, -t1 * t33 - t2 * t34 + t7 * t20 - t8 * t21 + t3 * t24 - t4 * t25 - t6 * t28 - t5 * t73, t1 * t8 + t2 * t7 + t3 * t6 + t4 * t5 + (t11 + t74) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t68, t59 + t80, t1 * t34 - t2 * t33 - t4 * t24 - t3 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, (t18 * t47 - t79) * qJD(1), 0, 0, 0, 0, 0, 0, 0, -t37, 0, (-t79 + (-qJD(4) - t16) * t47) * qJD(1), 0, 0, 0, 0, 0, 0, t40 + (-t28 - t63) * qJD(5), t60 * qJD(5), t82 + t83, -t3 * t28 - t4 * t73 - t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47 * t54 * t49, 0, -t45 * t54, (qJD(3) + t15) * t72, 0, 0, 0, 0, 0, 0, -t53 * t51 - t72 * t73, -t28 * t72 - t53 * t52, t52 * t20 - t51 * t21 + (t28 * t51 - t52 * t73) * qJD(5), -t11 * t72 + t1 * t51 + t2 * t52 + (-t3 * t51 + t4 * t52) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t82 - t83, 0, -t78, t40 + (t28 - t63) * qJD(5), 0, -t11 * t28 + t57, -(qJD(3) - t11) * t73, 0, 0;];
tauc_reg = t9;
