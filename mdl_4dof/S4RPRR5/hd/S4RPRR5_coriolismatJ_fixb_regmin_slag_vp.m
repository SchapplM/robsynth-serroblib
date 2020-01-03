% Calculate minimal parameter regressor of coriolis matrix for
% S4RPRR5
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
% cmat_reg [(4*%NQJ)%x16]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRR5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:41
% EndTime: 2019-12-31 16:51:42
% DurationCPUTime: 0.25s
% Computational Cost: add. (222->64), mult. (379->83), div. (0->0), fcn. (284->4), ass. (0->46)
t29 = sin(qJ(3));
t31 = cos(qJ(3));
t54 = -pkin(1) - pkin(2);
t17 = t31 * qJ(2) + t29 * t54;
t41 = qJD(1) - qJD(3);
t59 = t41 * t17;
t28 = sin(qJ(4));
t30 = cos(qJ(4));
t22 = -t28 ^ 2 + t30 ^ 2;
t58 = t41 * t22;
t16 = t29 * qJ(2) - t31 * t54;
t14 = pkin(3) + t16;
t57 = t14 + t16;
t53 = pkin(3) * t28;
t52 = qJD(1) * t14;
t51 = qJD(1) * t30;
t50 = qJD(3) * t30;
t49 = t22 * qJD(4);
t48 = t28 * qJD(4);
t47 = t29 * qJD(1);
t46 = t29 * qJD(2);
t45 = t30 * qJD(4);
t44 = t31 * qJD(1);
t43 = t31 * qJD(2);
t42 = qJ(2) * qJD(1);
t40 = t28 * t47;
t39 = t30 * t47;
t18 = t41 * t29;
t38 = t41 * t30;
t19 = t41 * t31;
t37 = t16 / 0.2e1 - pkin(3) / 0.2e1 - t14 / 0.2e1;
t36 = -t43 + t52;
t35 = t17 * qJD(1) + t46;
t34 = t17 * qJD(3) + t46;
t1 = t37 * t28;
t33 = t1 * qJD(1) + qJD(3) * t53;
t2 = t37 * t30;
t32 = pkin(3) * t50 + t2 * qJD(1);
t24 = t28 * t45;
t15 = -pkin(6) + t17;
t13 = (-t50 + t51) * t28;
t6 = t29 * t38 - t31 * t48;
t5 = t28 * t18 + t31 * t45;
t4 = (pkin(3) + t57) * t30 / 0.2e1;
t3 = t53 / 0.2e1 + t57 * t28 / 0.2e1;
t7 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), 0, t34, -t16 * qJD(3) + t43, t24, t49, 0, 0, 0, -t14 * t48 + t34 * t30, -t14 * t45 - t34 * t28; 0, 0, 0, 0, qJD(1), t42, 0, t47, t44, 0, 0, 0, 0, 0, t39, -t40; 0, 0, 0, 0, 0, 0, 0, t59, -t41 * t16, -t24, -t49, 0, 0, 0, t3 * qJD(4) + t17 * t38, t4 * qJD(4) - t28 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t58, -t45, t48, 0, t3 * qJD(3) - t15 * t45 - t28 * t52, t4 * qJD(3) - t14 * t51 + t15 * t48; 0, 0, 0, 0, -qJD(1), -t42, 0, -t18, -t19, 0, 0, 0, 0, 0, -t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t18, t19, 0, 0, 0, 0, 0, t6, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t19 - t29 * t45, t30 * t19 + t29 * t48; 0, 0, 0, 0, 0, 0, 0, -t35, t16 * qJD(1) - t43, -t24, -t49, 0, 0, 0, -t1 * qJD(4) - t35 * t30, -t2 * qJD(4) + t35 * t28; 0, 0, 0, 0, 0, 0, 0, -t47, -t44, 0, 0, 0, 0, 0, -t39, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t49, 0, 0, 0, -pkin(3) * t48, -pkin(3) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t58, t45, -t48, 0, -pkin(6) * t45 - t33, pkin(6) * t48 - t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t58, 0, 0, 0, t1 * qJD(3) + t36 * t28, t2 * qJD(3) + t36 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28 * t44, -t30 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t58, 0, 0, 0, t33, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t7;
