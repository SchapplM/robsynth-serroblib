% Calculate minimal parameter regressor of coriolis matrix for
% S4RPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:19
% EndTime: 2019-12-31 16:49:20
% DurationCPUTime: 0.24s
% Computational Cost: add. (245->46), mult. (477->62), div. (0->0), fcn. (470->6), ass. (0->41)
t39 = qJD(3) + qJD(4);
t29 = sin(qJ(4));
t30 = sin(qJ(3));
t31 = cos(qJ(3));
t53 = cos(qJ(4));
t20 = t29 * t30 - t53 * t31;
t59 = t39 * t20;
t25 = sin(pkin(7)) * pkin(1) + pkin(5);
t54 = pkin(6) + t25;
t18 = t54 * t30;
t19 = t54 * t31;
t58 = t39 * (t29 * t18 - t53 * t19);
t57 = t39 * (t53 * t18 + t29 * t19);
t56 = pkin(3) * t29;
t55 = pkin(3) * t30;
t22 = t29 * t31 + t53 * t30;
t26 = -cos(pkin(7)) * pkin(1) - pkin(2);
t23 = -t31 * pkin(3) + t26;
t5 = t20 * t55 + t23 * t22;
t48 = t5 * qJD(1);
t6 = -t23 * t20 + t22 * t55;
t47 = t6 * qJD(1);
t7 = t20 ^ 2 - t22 ^ 2;
t46 = t7 * qJD(1);
t45 = qJD(1) * t23;
t44 = qJD(1) * t31;
t43 = qJD(3) * t30;
t42 = qJD(3) * t31;
t41 = qJD(4) * t23;
t24 = -t30 ^ 2 + t31 ^ 2;
t40 = t24 * qJD(1);
t38 = t20 * t45;
t37 = t22 * t45;
t36 = t26 * t30 * qJD(1);
t35 = t26 * t44;
t34 = t30 * t44;
t33 = t53 * qJD(3);
t32 = t53 * qJD(4);
t12 = t39 * t22;
t13 = t22 * t20 * qJD(1);
t1 = [0, 0, 0, 0, t30 * t42, t24 * qJD(3), 0, 0, 0, t26 * t43, t26 * t42, -t20 * t12, t39 * t7, 0, 0, 0, t5 * qJD(3) + t22 * t41, t6 * qJD(3) - t20 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t34, t40, t42, -t43, 0, -t25 * t42 + t36, t25 * t43 + t35, -t13, t46, -t59, -t12, 0, t48 + t58, t47 + t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t46, -t59, -t12, 0, t37 + t58, -t38 + t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t42, 0, 0, 0, 0, 0, -t12, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t59; 0, 0, 0, 0, -t34, -t40, 0, 0, 0, -t36, -t35, t13, -t46, 0, 0, 0, -t48, -t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t56, -pkin(3) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39 * t56, (-t33 - t32) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t46, 0, 0, 0, -t37, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t56, pkin(3) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
