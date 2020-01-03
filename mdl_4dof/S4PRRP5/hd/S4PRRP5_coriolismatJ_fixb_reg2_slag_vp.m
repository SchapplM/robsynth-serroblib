% Calculate inertial parameters regressor of coriolis matrix for
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
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRRP5_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:29:25
% EndTime: 2019-12-31 16:29:25
% DurationCPUTime: 0.44s
% Computational Cost: add. (162->55), mult. (445->81), div. (0->0), fcn. (349->4), ass. (0->52)
t38 = sin(qJ(3));
t36 = t38 ^ 2;
t40 = cos(qJ(3));
t37 = t40 ^ 2;
t27 = t37 + t36;
t63 = t40 * pkin(3);
t33 = -pkin(2) - t63;
t62 = t33 * t38;
t61 = pkin(5) + qJ(4);
t1 = pkin(3) * t62;
t60 = t1 * qJD(2);
t39 = sin(qJ(2));
t41 = cos(qJ(2));
t9 = (-0.1e1 + t27) * t41 * t39;
t59 = t9 * qJD(1);
t14 = t38 * t63 - t62;
t58 = t14 * qJD(2);
t20 = t36 * pkin(3) + t33 * t40;
t57 = t20 * qJD(2);
t22 = t61 * t40;
t56 = t22 * qJD(3);
t55 = t27 * qJD(2);
t28 = t37 - t36;
t54 = t28 * qJD(2);
t53 = t38 * qJD(2);
t52 = t38 * qJD(3);
t51 = t38 * qJD(4);
t50 = t39 * qJD(2);
t49 = t40 * qJD(2);
t34 = t40 * qJD(3);
t48 = t41 * qJD(2);
t47 = pkin(2) * t53;
t46 = pkin(2) * t49;
t45 = pkin(3) * t53;
t44 = t39 * t34;
t43 = t27 * t41;
t21 = t61 * t38;
t5 = t21 * t38 + t22 * t40;
t10 = (0.1e1 / 0.2e1 - t37 / 0.2e1 - t36 / 0.2e1) * t39;
t42 = -t10 * qJD(1) + t5 * qJD(2);
t32 = t38 * t34;
t31 = t38 * t49;
t23 = t28 * qJD(3);
t19 = qJD(2) * t43;
t18 = -t39 * t49 - t41 * t52;
t17 = -t41 * t34 + t38 * t50;
t16 = -t38 * t48 - t44;
t15 = t39 * t52 - t40 * t48;
t11 = (0.1e1 + t27) * t39 / 0.2e1;
t6 = t9 * qJD(2);
t2 = t41 * t38 * pkin(3);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t48, 0, 0, 0, 0, 0, 0, 0, 0, t18, t17, t19, t59 + (-t39 * pkin(2) + pkin(5) * t43) * qJD(2), 0, 0, 0, 0, 0, 0, t18, t17, t19, t59 + (t39 * t33 + t5 * t41) * qJD(2) - t2 * qJD(3) + t11 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t15, 0, 0, 0, 0, 0, 0, 0, 0, t16, t15, 0, -pkin(3) * t44 - t2 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * qJD(4) - t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t23, 0, -t32, 0, 0, -pkin(2) * t52, -pkin(2) * t34, 0, 0, t32, t23, 0, -t32, 0, 0, -t14 * qJD(3), t20 * qJD(3), t27 * qJD(4), t1 * qJD(3) + t5 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t54, t34, -t31, -t52, 0, -pkin(5) * t34 - t47, pkin(5) * t52 - t46, 0, 0, t31, t54, t34, -t31, -t52, 0, -t56 - t58, t21 * qJD(3) + t57, -pkin(3) * t34, -pkin(3) * t56 + t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t54, 0, t31, 0, 0, t47, t46, 0, 0, -t31, -t54, 0, t31, 0, 0, -t51 + t58, -t40 * qJD(4) - t57, 0, -pkin(3) * t51 - t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t49, 0, -t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t34, -t55, pkin(3) * t52 - t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t49, 0, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
