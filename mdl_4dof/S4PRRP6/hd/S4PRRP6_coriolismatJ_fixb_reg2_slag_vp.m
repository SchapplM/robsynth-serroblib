% Calculate inertial parameters regressor of coriolis matrix for
% S4PRRP6
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
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRRP6_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:45
% EndTime: 2019-12-31 16:30:46
% DurationCPUTime: 0.38s
% Computational Cost: add. (130->57), mult. (356->75), div. (0->0), fcn. (280->4), ass. (0->51)
t30 = sin(qJ(3));
t28 = t30 ^ 2;
t32 = cos(qJ(3));
t29 = t32 ^ 2;
t58 = t28 + t29;
t60 = t30 * pkin(3);
t33 = cos(qJ(2));
t59 = t58 * pkin(5) * t33;
t37 = -t32 * pkin(3) - t30 * qJ(4);
t16 = -pkin(2) + t37;
t56 = t32 * qJ(4);
t17 = -t56 + t60;
t3 = t16 * t32 + t17 * t30;
t57 = t3 * qJD(2);
t4 = -t16 * t30 + t17 * t32;
t55 = t4 * qJD(2);
t31 = sin(qJ(2));
t8 = (-0.1e1 + t58) * t33 * t31;
t54 = t8 * qJD(1);
t53 = qJD(2) * t30;
t52 = qJD(2) * t32;
t51 = qJD(4) * t30;
t22 = t29 - t28;
t50 = t22 * qJD(2);
t49 = t22 * qJD(3);
t48 = t28 * qJD(2);
t47 = t30 * qJD(3);
t46 = t31 * qJD(2);
t27 = t32 * qJD(3);
t45 = t32 * qJD(4);
t44 = t33 * qJD(2);
t43 = qJD(3) * qJ(4);
t42 = pkin(2) * t53;
t41 = pkin(2) * t52;
t40 = pkin(5) * t47;
t39 = pkin(5) * t27;
t38 = t16 * t53;
t36 = -t60 / 0.2e1 + t56 / 0.2e1;
t1 = (t17 / 0.2e1 + t36) * t33;
t35 = -t16 * t17 * qJD(2) + t1 * qJD(1);
t34 = t37 * qJD(3) + t45;
t24 = t30 * t27;
t23 = t30 * t52;
t13 = t58 * t44;
t12 = -t32 * t46 - t33 * t47;
t11 = -t33 * t27 + t30 * t46;
t10 = t31 * t27 + t30 * t44;
t9 = t31 * t47 - t32 * t44;
t5 = t8 * qJD(2);
t2 = (-t17 / 0.2e1 + t36) * t33;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t44, 0, 0, 0, 0, 0, 0, 0, 0, t12, t11, t13, t54 + (-t31 * pkin(2) + t59) * qJD(2), 0, 0, 0, 0, 0, 0, t12, t13, -t11, t54 + (t31 * t16 + t59) * qJD(2) + t2 * qJD(3) + t33 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, -t9, t2 * qJD(2) + t34 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * qJD(3) - t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t49, 0, -t24, 0, 0, -pkin(2) * t47, -pkin(2) * t27, 0, 0, t24, 0, -t49, 0, 0, -t24, -t4 * qJD(3) + t30 * t45, 0, -t3 * qJD(3) + t28 * qJD(4), (qJD(3) * t17 - t51) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t50, t27, -t23, -t47, 0, -t39 - t42, t40 - t41, 0, 0, t23, t27, -t50, 0, t47, -t23, -t39 - t55, t34, -t40 - t57, t34 * pkin(5) - t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t27, t48, -t38 + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t50, 0, t23, 0, 0, t42, t41, 0, 0, -t23, 0, t50, 0, 0, t23, t55, 0, t57, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, -t48, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
