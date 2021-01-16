% Calculate minimal parameter regressor of coriolis matrix for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x16]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:56
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PPRRP4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:56:31
% EndTime: 2021-01-15 14:56:32
% DurationCPUTime: 0.27s
% Computational Cost: add. (164->45), mult. (439->78), div. (0->0), fcn. (331->4), ass. (0->48)
t35 = sin(qJ(4));
t33 = t35 ^ 2;
t37 = cos(qJ(4));
t34 = t37 ^ 2;
t24 = t33 + t34;
t59 = t37 * pkin(4);
t28 = -pkin(3) - t59;
t58 = t28 * t35;
t57 = qJ(5) + pkin(6);
t1 = pkin(4) * t58;
t56 = t1 * qJD(3);
t36 = sin(qJ(3));
t38 = cos(qJ(3));
t8 = (-0.1e1 + t24) * t38 * t36;
t55 = t8 * qJD(2);
t54 = qJD(3) * t38;
t30 = qJD(4) * t35;
t31 = qJD(4) * t37;
t53 = qJD(4) * t38;
t13 = t35 * t59 - t58;
t52 = t13 * qJD(3);
t20 = t33 * pkin(4) + t28 * t37;
t51 = t20 * qJD(3);
t22 = t57 * t37;
t50 = t22 * qJD(4);
t49 = t24 * qJD(3);
t25 = t34 - t33;
t48 = t25 * qJD(3);
t47 = t35 * qJD(3);
t46 = t35 * qJD(5);
t45 = t37 * qJD(3);
t44 = pkin(3) * t47;
t43 = pkin(3) * t45;
t42 = pkin(4) * t47;
t41 = t36 * t31;
t40 = t35 * t45;
t21 = t57 * t35;
t7 = t21 * t35 + t22 * t37;
t9 = (0.1e1 / 0.2e1 - t34 / 0.2e1 - t33 / 0.2e1) * t36;
t39 = -t9 * qJD(2) + t7 * qJD(3);
t29 = pkin(4) * t30;
t17 = -t35 * t53 - t36 * t45;
t16 = t36 * t47 - t37 * t53;
t15 = -t38 * t47 - t41;
t14 = t36 * t30 - t38 * t45;
t10 = (0.1e1 + t24) * t36 / 0.2e1;
t2 = t38 * t35 * pkin(4);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t31, t30, t31, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 * qJD(3); 0, 0, 0, -qJD(3) * t36, -t54, 0, 0, 0, 0, 0, t17, t16, t17, t16, t24 * t54, t55 + (t36 * t28 + t7 * t38) * qJD(3) - t2 * qJD(4) + t10 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t14, t15, t14, 0, -pkin(4) * t41 - t2 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9 * qJD(5) - t55; 0, 0, 0, 0, 0, t35 * t31, t25 * qJD(4), 0, 0, 0, -pkin(3) * t30, -pkin(3) * t31, -t13 * qJD(4), t20 * qJD(4), t24 * qJD(5), t1 * qJD(4) + t7 * qJD(5); 0, 0, 0, 0, 0, t40, t48, t31, -t30, 0, -pkin(6) * t31 - t44, pkin(6) * t30 - t43, -t50 - t52, t21 * qJD(4) + t51, -pkin(4) * t31, -pkin(4) * t50 + t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t40, -t48, 0, 0, 0, t44, t43, -t46 + t52, -t37 * qJD(5) - t51, 0, -pkin(4) * t46 - t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t45, 0, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t31, -t49, t29 - t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t45, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
