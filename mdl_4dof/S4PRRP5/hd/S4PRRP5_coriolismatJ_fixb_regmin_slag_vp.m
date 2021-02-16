% Calculate minimal parameter regressor of coriolis matrix for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x15]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 22:36
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRRP5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:36:10
% EndTime: 2021-01-14 22:36:11
% DurationCPUTime: 0.24s
% Computational Cost: add. (143->45), mult. (379->78), div. (0->0), fcn. (290->4), ass. (0->47)
t29 = sin(qJ(3));
t27 = t29 ^ 2;
t31 = cos(qJ(3));
t28 = t31 ^ 2;
t20 = t28 + t27;
t54 = t31 * pkin(3);
t24 = -pkin(2) - t54;
t53 = t24 * t29;
t52 = qJ(4) + pkin(5);
t1 = pkin(3) * t53;
t51 = t1 * qJD(2);
t30 = sin(qJ(2));
t32 = cos(qJ(2));
t6 = (-0.1e1 + t20) * t32 * t30;
t50 = t6 * qJD(1);
t49 = qJD(2) * t32;
t48 = qJD(3) * t29;
t25 = qJD(3) * t31;
t47 = qJD(3) * t32;
t11 = t29 * t54 - t53;
t46 = t11 * qJD(2);
t16 = t27 * pkin(3) + t24 * t31;
t45 = t16 * qJD(2);
t18 = t52 * t31;
t44 = t18 * qJD(3);
t43 = t20 * qJD(2);
t21 = t28 - t27;
t42 = t21 * qJD(2);
t41 = t29 * qJD(2);
t40 = t29 * qJD(4);
t39 = t31 * qJD(2);
t38 = pkin(2) * t41;
t37 = pkin(2) * t39;
t36 = pkin(3) * t41;
t35 = t30 * t25;
t34 = t29 * t39;
t17 = t52 * t29;
t5 = t17 * t29 + t18 * t31;
t7 = (0.1e1 / 0.2e1 - t28 / 0.2e1 - t27 / 0.2e1) * t30;
t33 = -t7 * qJD(1) + t5 * qJD(2);
t15 = -t29 * t47 - t30 * t39;
t14 = t30 * t41 - t31 * t47;
t13 = -t32 * t41 - t35;
t12 = t30 * t48 - t32 * t39;
t8 = (0.1e1 + t20) * t30 / 0.2e1;
t2 = t32 * t29 * pkin(3);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * qJD(2); 0, 0, -qJD(2) * t30, -t49, 0, 0, 0, 0, 0, t15, t14, t15, t14, t20 * t49, t50 + (t30 * t24 + t5 * t32) * qJD(2) - t2 * qJD(3) + t8 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t12, t13, t12, 0, -pkin(3) * t35 - t2 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * qJD(4) - t50; 0, 0, 0, 0, t29 * t25, t21 * qJD(3), 0, 0, 0, -pkin(2) * t48, -pkin(2) * t25, -t11 * qJD(3), t16 * qJD(3), t20 * qJD(4), t1 * qJD(3) + t5 * qJD(4); 0, 0, 0, 0, t34, t42, t25, -t48, 0, -pkin(5) * t25 - t38, pkin(5) * t48 - t37, -t44 - t46, t17 * qJD(3) + t45, -pkin(3) * t25, -pkin(3) * t44 + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t34, -t42, 0, 0, 0, t38, t37, -t40 + t46, -t31 * qJD(4) - t45, 0, -pkin(3) * t40 - t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t39, 0, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t25, -t43, pkin(3) * t48 - t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t39, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
