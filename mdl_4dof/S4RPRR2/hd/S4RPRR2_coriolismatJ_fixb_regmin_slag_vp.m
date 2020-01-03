% Calculate minimal parameter regressor of coriolis matrix for
% S4RPRR2
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
% cmat_reg [(4*%NQJ)%x14]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:14
% EndTime: 2019-12-31 16:48:14
% DurationCPUTime: 0.20s
% Computational Cost: add. (188->41), mult. (398->55), div. (0->0), fcn. (306->6), ass. (0->40)
t23 = cos(pkin(7)) * pkin(1) + pkin(2);
t29 = sin(qJ(3));
t31 = cos(qJ(3));
t44 = pkin(1) * sin(pkin(7));
t16 = -t31 * t23 + t29 * t44;
t53 = t16 / 0.2e1;
t14 = -pkin(3) + t16;
t52 = -pkin(3) / 0.2e1 + t14 / 0.2e1;
t51 = t53 + t52;
t30 = cos(qJ(4));
t38 = qJD(1) + qJD(3);
t50 = t30 * t38;
t28 = sin(qJ(4));
t21 = -t28 ^ 2 + t30 ^ 2;
t49 = t38 * t21;
t43 = pkin(3) * qJD(3);
t42 = qJD(1) * t14;
t41 = t16 * qJD(1);
t17 = t29 * t23 + t31 * t44;
t40 = t17 * qJD(1);
t13 = t17 * qJD(3);
t39 = t28 * qJD(4);
t26 = t30 * qJD(4);
t37 = t28 * t42;
t36 = t30 * t42;
t35 = t28 * t40;
t34 = t53 - t52;
t1 = t34 * t28;
t33 = t1 * qJD(1) + t28 * t43;
t2 = t34 * t30;
t32 = t2 * qJD(1) + t30 * t43;
t22 = t28 * t26;
t20 = t21 * qJD(4);
t18 = t28 * t50;
t15 = pkin(6) + t17;
t12 = t16 * qJD(3);
t9 = t28 * t13;
t4 = t51 * t30;
t3 = t51 * t28;
t5 = [0, 0, 0, 0, 0, -t13, t12, t22, t20, 0, 0, 0, -t30 * t13 + t14 * t39, t14 * t26 + t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t13 - t40, t12 + t41, t22, t20, 0, 0, 0, t3 * qJD(4) - t17 * t50, t4 * qJD(4) + t35 + t9; 0, 0, 0, 0, 0, 0, 0, t18, t49, t26, -t39, 0, t3 * qJD(3) - t15 * t26 + t37, t4 * qJD(3) + t15 * t39 + t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t26; 0, 0, 0, 0, 0, t40, -t41, t22, t20, 0, 0, 0, -t1 * qJD(4) + t30 * t40, -t2 * qJD(4) - t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t22, t20, 0, 0, 0, -pkin(3) * t39, -pkin(3) * t26; 0, 0, 0, 0, 0, 0, 0, t18, t49, t26, -t39, 0, -pkin(6) * t26 - t33, pkin(6) * t39 - t32; 0, 0, 0, 0, 0, 0, 0, -t18, -t49, 0, 0, 0, t1 * qJD(3) - t37, t2 * qJD(3) - t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t18, -t49, 0, 0, 0, t33, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
