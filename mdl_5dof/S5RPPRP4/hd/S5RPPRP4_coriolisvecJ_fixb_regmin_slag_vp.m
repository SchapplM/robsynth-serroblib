% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% tauc_reg [5x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:20
% EndTime: 2019-12-31 17:52:21
% DurationCPUTime: 0.29s
% Computational Cost: add. (395->77), mult. (764->122), div. (0->0), fcn. (374->4), ass. (0->66)
t31 = sin(pkin(7));
t37 = qJD(4) ^ 2;
t38 = qJD(1) ^ 2;
t79 = t31 * (t37 + t38);
t36 = -pkin(1) - pkin(2);
t19 = t36 * qJD(1) + qJD(2);
t32 = cos(pkin(7));
t58 = qJD(1) * qJ(2);
t14 = t31 * t19 + t32 * t58;
t10 = -qJD(1) * pkin(6) + t14;
t34 = sin(qJ(4));
t35 = cos(qJ(4));
t78 = t35 * qJD(3) - t34 * t10;
t61 = qJD(2) * t32;
t49 = -qJD(5) + t61;
t41 = t49 * t35;
t59 = qJ(5) * qJD(4);
t1 = t78 * qJD(4) + (t34 * t59 + t41) * qJD(1);
t42 = t49 * t34;
t44 = -t34 * qJD(3) - t35 * t10;
t2 = t44 * qJD(4) + (t35 * t59 - t42) * qJD(1);
t60 = qJ(5) * qJD(1);
t6 = t34 * t60 + t78;
t64 = qJD(4) * pkin(4);
t3 = t6 + t64;
t7 = -t35 * t60 - t44;
t77 = -t1 * t35 + t2 * t34 + (t3 * t35 + t34 * t7) * qJD(4);
t57 = qJD(1) * qJD(2);
t20 = t31 * t57;
t76 = 0.2e1 * t20;
t75 = t3 - t6;
t74 = t35 * t7;
t73 = t32 * t38;
t71 = t37 * t34;
t70 = t37 * t35;
t69 = t32 * qJ(2) + t31 * t36;
t29 = t34 ^ 2;
t30 = t35 ^ 2;
t68 = t29 - t30;
t67 = t29 + t30;
t62 = qJD(1) * t35;
t13 = t32 * t19 - t31 * t58;
t9 = qJD(1) * pkin(3) - t13;
t8 = pkin(4) * t62 + qJD(5) + t9;
t65 = qJD(1) * t8;
t17 = -pkin(6) + t69;
t63 = qJ(5) - t17;
t56 = qJD(1) * qJD(4);
t55 = 0.2e1 * t57;
t54 = 0.2e1 * t56;
t53 = t34 * t56;
t52 = -t31 * qJ(2) + t32 * t36;
t51 = qJD(4) * t63;
t50 = t35 * t54;
t16 = pkin(3) - t52;
t46 = t3 * t34 - t74;
t45 = t13 * t31 - t14 * t32;
t43 = (t9 - t61) * qJD(1);
t40 = -t17 * t37 + t76;
t39 = qJD(4) * (-qJD(1) * t16 - t61 - t9);
t15 = -pkin(4) * t53 + t20;
t12 = t63 * t35;
t11 = t63 * t34;
t5 = t35 * t51 - t42;
t4 = t34 * t51 + t41;
t18 = [0, 0, 0, 0, t55, qJ(2) * t55, t76, t32 * t55, ((-t31 * t52 + t32 * t69) * qJD(1) - t45) * qJD(2), t34 * t50, -t68 * t54, -t70, t71, 0, t34 * t39 + t40 * t35, -t40 * t34 + t35 * t39, (t34 * t5 - t35 * t4 + (t11 * t35 - t12 * t34) * qJD(4)) * qJD(1) + t77, -t1 * t12 + t7 * t4 + t2 * t11 + t3 * t5 + t15 * (t35 * pkin(4) + t16) + t8 * (t31 * qJD(2) - t34 * t64); 0, 0, 0, 0, -t38, -t38 * qJ(2), -t31 * t38, -t73, t45 * qJD(1), 0, 0, 0, 0, 0, 0.2e1 * t32 * t53 - t35 * t79, t32 * t50 + t34 * t79, t67 * t73, (t46 * qJD(1) - t15) * t32 + (-t65 - t77) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t70, 0, -t46 * qJD(4) + t1 * t34 + t2 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34 * t38 * t35, t68 * t38, 0, 0, 0, t34 * t43, t35 * t43, (t64 - t75) * t62, t75 * t7 + (t34 * t65 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67 * t38, t20 + (t74 + (-t3 - t64) * t34) * qJD(1);];
tauc_reg = t18;
