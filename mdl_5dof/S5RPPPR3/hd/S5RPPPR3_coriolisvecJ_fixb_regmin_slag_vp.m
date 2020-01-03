% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% tauc_reg [5x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:44:01
% EndTime: 2019-12-31 17:44:02
% DurationCPUTime: 0.24s
% Computational Cost: add. (195->57), mult. (515->93), div. (0->0), fcn. (336->6), ass. (0->48)
t36 = sin(pkin(8));
t34 = t36 ^ 2;
t38 = cos(pkin(8));
t60 = t38 ^ 2 + t34;
t67 = t60 * qJD(1) * qJD(3);
t41 = cos(qJ(5));
t40 = sin(qJ(5));
t61 = t38 * t40;
t24 = t36 * t41 - t61;
t47 = cos(pkin(7)) * pkin(1) + t36 * qJ(4) + pkin(2);
t66 = (t38 * pkin(3) + t47) * qJD(1);
t32 = sin(pkin(7)) * pkin(1) + qJ(3);
t28 = t32 * qJD(1);
t65 = t38 * (t36 * qJD(2) + t38 * t28);
t64 = -pkin(6) + t32;
t63 = t36 * t38;
t23 = t36 * t40 + t38 * t41;
t59 = qJD(1) * t23;
t58 = qJD(1) * t36;
t57 = qJD(3) * t36;
t56 = qJD(4) * t36;
t14 = t23 * qJD(5);
t55 = t14 * qJD(5);
t53 = qJD(3) * t65 + t32 * t67;
t52 = qJD(1) * t61;
t51 = t41 * t58;
t50 = 0.2e1 * qJD(1) * qJD(4);
t48 = 0.2e1 * t59;
t8 = t38 * qJD(2) - t36 * t28;
t46 = t24 * qJD(3);
t45 = t23 * qJD(3);
t7 = (pkin(3) + pkin(4)) * t38 + t47;
t43 = qJD(1) ^ 2;
t42 = qJD(5) ^ 2;
t29 = qJD(5) * t52;
t27 = t60 * t43;
t20 = t64 * t38;
t19 = t64 * t36;
t18 = t51 - t52;
t15 = t24 * qJD(5);
t13 = 0.2e1 * t67;
t12 = t15 * qJD(5);
t11 = qJD(5) * t51 - t29;
t10 = qJD(1) * t14;
t6 = qJD(4) - t8;
t5 = qJD(3) - t66;
t1 = t7 * qJD(1) - qJD(3);
t2 = [0, 0, 0, 0, 0, 0, t13, -t8 * t57 + t53, t50 * t63, t13, t34 * t50, t6 * t57 + (-t5 + t66) * t56 + t53, -t10 * t24 - t18 * t14, t10 * t23 - t24 * t11 + t14 * t59 - t18 * t15, -t55, -t12, 0, t1 * t15 + t7 * t11 + ((-t19 * t40 - t20 * t41) * qJD(5) + t46) * qJD(5) + t48 * t56, -t1 * t14 - t7 * t10 + ((-t19 * t41 + t20 * t40) * qJD(5) - t45) * qJD(5) + (qJD(1) * t24 + t18) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t55; 0, 0, 0, 0, 0, 0, -t27, (t36 * t8 - t65) * qJD(1), 0, -t27, 0, (-t65 + (-qJD(4) - t6) * t36) * qJD(1), 0, 0, 0, 0, 0, t29 + (-t18 - t51) * qJD(5), t48 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, -t43 * t63, 0, -t34 * t43, (qJD(3) + t5) * t58, 0, 0, 0, 0, 0, -t42 * t40 - t58 * t59, -t18 * t58 - t42 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 * t59, t18 ^ 2 - t59 ^ 2, 0, t29 + (t18 - t51) * qJD(5), 0, qJD(1) * t46 - t1 * t18, -qJD(1) * t45 + t1 * t59;];
tauc_reg = t2;
