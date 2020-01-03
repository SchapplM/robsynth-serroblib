% Calculate minimal parameter regressor of coriolis matrix for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x20]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRR8_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:17
% EndTime: 2019-12-31 16:55:18
% DurationCPUTime: 0.26s
% Computational Cost: add. (242->52), mult. (448->70), div. (0->0), fcn. (444->4), ass. (0->45)
t37 = qJD(3) + qJD(4);
t29 = -pkin(1) - pkin(5);
t61 = -pkin(6) + t29;
t27 = sin(qJ(3));
t20 = t61 * t27;
t28 = cos(qJ(3));
t21 = t61 * t28;
t26 = sin(qJ(4));
t56 = cos(qJ(4));
t60 = t37 * (-t20 * t56 - t26 * t21);
t59 = t37 * (t26 * t20 - t21 * t56);
t58 = pkin(3) * t26;
t57 = pkin(3) * t28;
t17 = t26 * t27 - t28 * t56;
t18 = t26 * t28 + t27 * t56;
t5 = -t17 ^ 2 + t18 ^ 2;
t51 = t5 * qJD(1);
t23 = pkin(3) * t27 + qJ(2);
t6 = -t17 * t23 + t18 * t57;
t50 = t6 * qJD(1);
t7 = -t17 * t57 - t18 * t23;
t49 = t7 * qJD(1);
t48 = qJD(3) * t27;
t47 = qJD(3) * t28;
t46 = qJD(3) * t29;
t45 = qJD(4) * t23;
t44 = t17 * qJD(1);
t43 = t18 * qJD(1);
t22 = t27 ^ 2 - t28 ^ 2;
t42 = t22 * qJD(1);
t41 = t27 * qJD(1);
t40 = t28 * qJD(1);
t39 = qJ(2) * qJD(3);
t38 = qJD(1) * qJ(2);
t36 = t23 * t44;
t35 = t23 * t43;
t34 = t27 * t40;
t33 = t27 * t38;
t32 = t28 * t38;
t31 = t56 * qJD(3);
t30 = t56 * qJD(4);
t10 = t37 * t17;
t12 = t17 * t43;
t11 = t37 * t18;
t1 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t27 * t47, t22 * qJD(3), 0, 0, 0, qJD(2) * t27 + t28 * t39, qJD(2) * t28 - t27 * t39, t18 * t10, t37 * t5, 0, 0, 0, qJD(2) * t18 + qJD(3) * t6 - t17 * t45, -qJD(2) * t17 + qJD(3) * t7 - t18 * t45; 0, 0, 0, 0, qJD(1), t38, 0, 0, 0, 0, 0, t41, t40, 0, 0, 0, 0, 0, t43, -t44; 0, 0, 0, 0, 0, 0, -t34, t42, -t48, -t47, 0, -t27 * t46 + t32, -t28 * t46 - t33, t12, t51, -t11, t10, 0, t50 + t60, t49 + t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t51, -t11, t10, 0, -t36 + t60, -t35 + t59; 0, 0, 0, 0, -qJD(1), -t38, 0, 0, 0, 0, 0, -t41, -t40, 0, 0, 0, 0, 0, -t43, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t47, 0, 0, 0, 0, 0, -t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10; 0, 0, 0, 0, 0, 0, t34, -t42, 0, 0, 0, -t32, t33, -t12, -t51, 0, 0, 0, -t50, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t58, -pkin(3) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37 * t58, (-t31 - t30) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t51, 0, 0, 0, t36, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t58, pkin(3) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
