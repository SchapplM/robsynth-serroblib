% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPPR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:37
% EndTime: 2019-12-31 17:38:38
% DurationCPUTime: 0.27s
% Computational Cost: add. (230->68), mult. (498->105), div. (0->0), fcn. (307->6), ass. (0->50)
t32 = sin(pkin(8));
t37 = cos(qJ(2));
t54 = t37 * qJD(1);
t19 = (qJD(3) + t54) * qJD(2);
t33 = cos(pkin(8));
t35 = sin(qJ(2));
t56 = t35 * qJD(1);
t51 = qJD(2) * t56;
t5 = t32 * t19 - t33 * t51;
t8 = t32 * t54 - t33 * t56;
t72 = (t32 * qJD(3) - t8) * qJD(2) + t5;
t39 = qJD(5) ^ 2;
t40 = qJD(2) ^ 2;
t71 = t32 * (t39 + t40);
t38 = -pkin(2) - pkin(3);
t61 = t33 * qJ(3) + t32 * t38;
t70 = -(-pkin(6) + t61) * t39 + t72;
t14 = -t37 * t32 + t35 * t33;
t69 = 0.2e1 * qJD(2);
t23 = qJD(2) * qJ(3) + t56;
t68 = t23 * t37;
t34 = sin(qJ(5));
t65 = t39 * t34;
t36 = cos(qJ(5));
t64 = t39 * t36;
t63 = t40 * t35;
t62 = t40 * t37;
t6 = t33 * t19 + t32 * t51;
t60 = t34 ^ 2 - t36 ^ 2;
t58 = qJD(2) * pkin(2);
t57 = t34 * qJD(5);
t55 = t36 * qJD(5);
t53 = qJD(5) * t69;
t46 = qJD(3) - t54;
t15 = t38 * qJD(2) + t46;
t3 = t33 * t15 - t32 * t23;
t1 = qJD(2) * pkin(4) - t3;
t50 = t1 * qJD(2) - t6;
t44 = t35 * t32 + t37 * t33;
t9 = t44 * qJD(1);
t48 = t33 * qJD(3) - t9;
t47 = t36 * t53;
t4 = t32 * t15 + t33 * t23;
t45 = t3 * t32 - t4 * t33;
t43 = -t32 * qJ(3) + t33 * t38;
t41 = qJD(5) * (-qJD(2) * (pkin(4) - t43) - t1 - t48);
t20 = t46 - t58;
t11 = t44 * qJD(2);
t10 = t14 * qJD(2);
t2 = [0, 0, -t63, -t62, -t63, t62, t19 * t35 + (t68 + (t20 - t54) * t35) * qJD(2), -t10 * qJD(2), t11 * qJD(2), t3 * t10 + t4 * t11 + t6 * t14 + t44 * t5, 0, 0, 0, 0, 0, -t11 * t57 - t14 * t64 + (-t10 * t36 - t44 * t57) * qJD(2), -t11 * t55 + t14 * t65 + (t10 * t34 - t44 * t55) * qJD(2); 0, 0, 0, 0, 0, qJD(3) * t69, t19 * qJ(3) + t23 * qJD(3) + (-t68 + (-t20 - t58) * t35) * qJD(1), t72, t48 * qJD(2) + t6, -t45 * qJD(3) + t3 * t8 - t4 * t9 - t5 * t43 + t6 * t61, t34 * t47, -t60 * t53, -t64, t65, 0, t34 * t41 + t70 * t36, -t70 * t34 + t36 * t41; 0, 0, 0, 0, 0, -t40, (-t23 + t56) * qJD(2), -t32 * t40, -t33 * t40, t45 * qJD(2) + t6 * t32 - t5 * t33, 0, 0, 0, 0, 0, t33 * t34 * t53 - t36 * t71, t33 * t47 + t34 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34 * t40 * t36, t60 * t40, 0, 0, 0, t50 * t34, t50 * t36;];
tauc_reg = t2;
