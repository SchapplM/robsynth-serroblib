% Calculate inertial parameters regressor of coriolis matrix for
% S4RPRR1
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
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRR1_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:31:56
% EndTime: 2019-03-08 18:31:56
% DurationCPUTime: 0.24s
% Computational Cost: add. (550->48), mult. (1058->64), div. (0->0), fcn. (894->6), ass. (0->47)
t62 = -pkin(3) / 0.2e1;
t42 = cos(pkin(7)) * pkin(1) + pkin(2);
t60 = cos(qJ(3));
t34 = t60 * t42;
t44 = sin(pkin(7)) * pkin(1);
t59 = sin(qJ(3));
t41 = t59 * t44;
t32 = t41 - t34;
t61 = t32 / 0.2e1;
t37 = sin(qJ(4));
t58 = t37 * t32;
t33 = t59 * t42 + t60 * t44;
t57 = t37 * t33;
t38 = cos(qJ(4));
t56 = t38 * t32;
t29 = t38 * t33;
t55 = pkin(3) * qJD(3);
t54 = pkin(3) * qJD(4);
t40 = pkin(3) - t32;
t28 = t38 * t40;
t16 = -t28 + t57;
t39 = t37 * t40;
t17 = t39 + t29;
t18 = t29 - t58;
t19 = -t56 - t57;
t2 = t16 * t18 + t17 * t19;
t53 = t2 * qJD(1);
t5 = (t61 + pkin(3) - t41 / 0.2e1 + t34 / 0.2e1) * t37;
t52 = t5 * qJD(1);
t7 = t28 / 0.2e1 + (t61 + pkin(3) / 0.2e1) * t38;
t51 = t7 * qJD(1);
t50 = t16 * qJD(1);
t49 = t17 * qJD(1);
t48 = t18 * qJD(1);
t47 = t19 * qJD(1);
t46 = t32 * qJD(1);
t45 = t33 * qJD(1);
t43 = pkin(3) * (-qJD(3) - qJD(4));
t31 = t33 * qJD(3);
t30 = t32 * qJD(3);
t15 = t19 * qJD(3);
t14 = t18 * qJD(3);
t13 = t17 * qJD(4);
t12 = t16 * qJD(4);
t8 = t38 * t62 + t57 - t28 / 0.2e1 + t56 / 0.2e1;
t6 = t37 * t62 - t29 - t39 / 0.2e1 + t58 / 0.2e1;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t30, 0, 0, 0, 0, 0, 0, 0, 0, -t14 - t13, -t15 + t12, 0, t2 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31 - t45, t30 + t46, 0, 0, 0, 0, 0, 0, 0, 0, t6 * qJD(4) - t14 - t48, t8 * qJD(4) - t15 - t47, 0, t53 + (-t18 * t38 + t19 * t37) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * qJD(3) - t13 - t49, t8 * qJD(3) + t12 + t50, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t46, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * qJD(4) + t48, -t7 * qJD(4) + t47, 0, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37 * t54, -t38 * t54, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t43 - t52, t38 * t43 - t51, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * qJD(3) + t49, t7 * qJD(3) - t50, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t55 + t52, t38 * t55 + t51, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t1;
