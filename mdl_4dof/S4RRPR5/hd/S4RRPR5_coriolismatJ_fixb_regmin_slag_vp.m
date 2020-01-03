% Calculate minimal parameter regressor of coriolis matrix for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x16]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPR5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:34
% EndTime: 2019-12-31 17:03:35
% DurationCPUTime: 0.27s
% Computational Cost: add. (151->70), mult. (337->69), div. (0->0), fcn. (213->4), ass. (0->54)
t34 = sin(qJ(2));
t63 = t34 * pkin(1);
t24 = qJ(3) + t63;
t67 = -qJ(3) / 0.2e1 - t24 / 0.2e1;
t33 = sin(qJ(4));
t35 = cos(qJ(4));
t22 = t33 ^ 2 - t35 ^ 2;
t31 = qJD(1) + qJD(2);
t66 = t31 * t22;
t36 = cos(qJ(2));
t62 = t36 * pkin(1);
t57 = pkin(1) * qJD(2);
t28 = t36 * t57;
t29 = qJD(3) * t33;
t60 = t33 * t28 + t29;
t30 = qJD(3) * t35;
t59 = t35 * t28 + t30;
t58 = pkin(1) * qJD(1);
t43 = -pkin(2) - t62;
t5 = -t24 * t62 - t43 * t63;
t56 = t5 * qJD(1);
t55 = t24 * qJD(1);
t54 = t33 * qJD(4);
t53 = t35 * qJD(4);
t52 = t28 + qJD(3);
t51 = qJ(3) * qJD(2);
t50 = qJ(3) * qJD(4);
t49 = t34 * t57;
t48 = t34 * t58;
t27 = t36 * t58;
t47 = t63 / 0.2e1;
t46 = t33 * t55;
t45 = t35 * t55;
t44 = t33 * t53;
t42 = t33 * t27;
t41 = t35 * t27;
t13 = t31 * t35;
t40 = t47 + t67;
t1 = t40 * t33;
t39 = -t1 * qJD(1) + t33 * t51;
t2 = t40 * t35;
t38 = t2 * qJD(1) - t35 * t51;
t37 = -pkin(2) - pkin(6);
t32 = qJ(3) * qJD(3);
t23 = -pkin(6) + t43;
t21 = t31 * qJ(3);
t20 = t24 * qJD(3);
t15 = t22 * qJD(4);
t12 = t31 * t33;
t9 = t31 * t63;
t8 = t33 * t13;
t4 = (t47 - t67) * t35;
t3 = (-t63 / 0.2e1 + t67) * t33;
t6 = [0, 0, 0, 0, -t49, -t28, t49, t52, -t5 * qJD(2) + t20, -t44, t15, 0, 0, 0, t24 * t53 + t60, -t24 * t54 + t59; 0, 0, 0, 0, -t9, -t27 - t28, t9, t27 + t52, -t56 + t20 + (-pkin(2) * t34 + qJ(3) * t36) * t57, -t44, t15, 0, 0, 0, t4 * qJD(4) + t42 + t60, t3 * qJD(4) + t41 + t59; 0, 0, 0, 0, 0, 0, 0, t31, t31 * t24, 0, 0, 0, 0, 0, t12, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t66, -t54, -t53, 0, t4 * qJD(2) - t23 * t54 + t45, t3 * qJD(2) - t23 * t53 - t46; 0, 0, 0, 0, t48, t27, -t48, -t27 + qJD(3), t32 + t56, -t44, t15, 0, 0, 0, -t2 * qJD(4) + t29 - t42, t1 * qJD(4) + t30 - t41; 0, 0, 0, 0, 0, 0, 0, qJD(3), t32, -t44, t15, 0, 0, 0, t35 * t50 + t29, -t33 * t50 + t30; 0, 0, 0, 0, 0, 0, 0, t31, t21, 0, 0, 0, 0, 0, t12, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t66, -t54, -t53, 0, -t37 * t54 - t38, -t37 * t53 - t39; 0, 0, 0, 0, 0, 0, 0, -t31, -t51 - t55, 0, 0, 0, 0, 0, -t12, -t13; 0, 0, 0, 0, 0, 0, 0, -t31, -t21, 0, 0, 0, 0, 0, -t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t66, 0, 0, 0, t2 * qJD(2) - t45, -t1 * qJD(2) + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t66, 0, 0, 0, t38, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
