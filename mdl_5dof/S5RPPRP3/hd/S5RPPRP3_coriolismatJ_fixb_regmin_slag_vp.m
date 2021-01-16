% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRP3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 17:04:44
% EndTime: 2021-01-15 17:04:45
% DurationCPUTime: 0.30s
% Computational Cost: add. (287->59), mult. (429->59), div. (0->0), fcn. (317->4), ass. (0->44)
t21 = -cos(pkin(7)) * pkin(1) - pkin(2) - pkin(6);
t56 = -qJ(5) + t21;
t34 = sin(qJ(4));
t55 = t34 * pkin(4);
t22 = sin(pkin(7)) * pkin(1) + qJ(3);
t16 = t22 + t55;
t35 = cos(qJ(4));
t1 = t16 * t35 * pkin(4);
t54 = t1 * qJD(1);
t14 = t56 * t34;
t15 = t56 * t35;
t6 = t14 * t34 + t15 * t35;
t53 = t6 * qJD(1);
t8 = (t16 + t55) * t35;
t52 = t8 * qJD(1);
t32 = t35 ^ 2;
t13 = t32 * pkin(4) - t16 * t34;
t51 = t13 * qJD(1);
t50 = t14 * qJD(4);
t49 = t16 * qJD(1);
t31 = t34 ^ 2;
t36 = -t31 / 0.2e1 - t32 / 0.2e1;
t18 = -0.1e1 / 0.2e1 + t36;
t48 = t18 * qJD(1);
t19 = t31 - t32;
t47 = t19 * qJD(1);
t20 = -t31 - t32;
t46 = t20 * qJD(1);
t45 = t22 * qJD(1);
t26 = t34 * qJD(1);
t25 = t34 * qJD(4);
t27 = t35 * qJD(1);
t44 = t35 * qJD(4);
t43 = t35 * qJD(5);
t42 = pkin(4) * t25;
t41 = pkin(4) * t44;
t40 = pkin(4) * t27;
t39 = t22 * t26;
t38 = t22 * t27;
t37 = t34 * t27;
t29 = qJD(3) * t35;
t28 = qJD(3) * t34;
t17 = 0.1e1 / 0.2e1 + t36;
t2 = [0, 0, 0, 0, 0, qJD(3), t22 * qJD(3), -t34 * t44, t19 * qJD(4), 0, 0, 0, t22 * t44 + t28, -t22 * t25 + t29, t8 * qJD(4) + t28, t13 * qJD(4) + t29, -t20 * qJD(5), t16 * qJD(3) + t1 * qJD(4) - t6 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, qJD(1), t45, 0, 0, 0, 0, 0, t26, t27, t26, t27, 0, t17 * qJD(5) + t49; 0, 0, 0, 0, 0, 0, 0, -t37, t47, -t25, -t44, 0, -t21 * t25 + t38, -t21 * t44 - t39, -t50 + t52, -t15 * qJD(4) + t51, t42, -pkin(4) * t50 + t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t17 * qJD(3) - t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t25, -t44, t25, 0, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -qJD(1), -t45, 0, 0, 0, 0, 0, -t26, -t27, -t26, -t27, 0, t18 * qJD(5) - t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t44, -t25, -t44, 0, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48; 0, 0, 0, 0, 0, 0, 0, t37, -t47, 0, 0, 0, -t38, t39, -t43 - t52, t34 * qJD(5) - t51, 0, -pkin(4) * t43 - t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t26, 0, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t25, t46, -t18 * qJD(3) + t41 + t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t26, 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
