% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x20]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPPR5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:32
% EndTime: 2019-12-31 17:46:33
% DurationCPUTime: 0.38s
% Computational Cost: add. (242->65), mult. (502->102), div. (0->0), fcn. (487->6), ass. (0->62)
t36 = sin(pkin(8));
t34 = t36 ^ 2;
t38 = cos(pkin(8));
t35 = t38 ^ 2;
t26 = t34 + t35;
t37 = sin(pkin(7));
t39 = cos(pkin(7));
t42 = -pkin(1) - pkin(2);
t69 = t39 * qJ(2) + t37 * t42;
t21 = -qJ(4) + t69;
t74 = pkin(6) - t21;
t40 = sin(qJ(5));
t73 = t40 * t36;
t72 = t40 * t38;
t41 = cos(qJ(5));
t71 = t41 * t36;
t70 = t41 * t38;
t18 = t26 * t39;
t47 = -t37 * qJ(2) + t39 * t42;
t45 = pkin(3) - t47;
t1 = t18 * t21 + t37 * t45;
t68 = t1 * qJD(1);
t19 = -t70 + t73;
t20 = t71 + t72;
t2 = t19 ^ 2 - t20 ^ 2;
t67 = t2 * qJD(1);
t44 = -t72 / 0.2e1 - t71 / 0.2e1;
t3 = (-t20 / 0.2e1 + t44) * t39;
t66 = t3 * qJD(1);
t43 = -t70 / 0.2e1 + t73 / 0.2e1;
t4 = (t19 / 0.2e1 + t43) * t39;
t65 = t4 * qJD(1);
t7 = t26 * t21;
t64 = t7 * qJD(1);
t8 = -t47 * t37 + t69 * t39;
t63 = t8 * qJD(1);
t11 = (0.1e1 / 0.2e1 + t35 / 0.2e1 + t34 / 0.2e1) * t37;
t62 = t11 * qJD(1);
t61 = t18 * qJD(1);
t60 = t19 * qJD(1);
t16 = t19 * qJD(5);
t59 = t20 * qJD(1);
t17 = t20 * qJD(5);
t58 = t26 * qJD(1);
t57 = t37 * qJD(1);
t56 = t37 * qJD(2);
t55 = t39 * qJD(1);
t54 = t39 * qJD(2);
t53 = qJ(2) * qJD(1);
t52 = t19 * t59;
t51 = t19 * t57;
t50 = t20 * t57;
t49 = t36 * t57;
t48 = t38 * t57;
t15 = t38 * pkin(4) + t45;
t46 = qJD(1) * t15 + qJD(4);
t12 = (0.1e1 / 0.2e1 - t26 / 0.2e1) * t37;
t10 = t74 * t38;
t9 = t74 * t36;
t6 = (t20 / 0.2e1 + t44) * t39;
t5 = (-t19 / 0.2e1 + t43) * t39;
t13 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), t56, t54, t8 * qJD(2), t38 * t56, -t36 * t56, -t18 * qJD(2) + t26 * qJD(4), t1 * qJD(2) - t7 * qJD(4), -t19 * t17, t2 * qJD(5), 0, 0, 0, -t15 * t17 - t19 * t56, t15 * t16 - t20 * t56; 0, 0, 0, 0, qJD(1), t53, t57, t55, t63, t48, -t49, -t61, t68 + t12 * qJD(4) + (-0.1e1 + t26) * t37 * t54, 0, 0, 0, 0, 0, t6 * qJD(5) - t51, t5 * qJD(5) - t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t12 * qJD(2) - t64, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, t67, t16, t17, 0, -t15 * t59 + t6 * qJD(2) + (t41 * t10 - t40 * t9) * qJD(5), t15 * t60 + t5 * qJD(2) + (-t40 * t10 - t41 * t9) * qJD(5); 0, 0, 0, 0, -qJD(1), -t53, -t57, -t55, -t63, -t48, t49, t61, -t11 * qJD(4) - t68, 0, 0, 0, 0, 0, -t3 * qJD(5) + t51, -t4 * qJD(5) + t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t16 - t66, t37 * t17 - t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t11 * qJD(2) + t64, 0, 0, 0, 0, 0, -t17, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t67, 0, 0, 0, t3 * qJD(2) + t46 * t20, t4 * qJD(2) - t46 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t13;
