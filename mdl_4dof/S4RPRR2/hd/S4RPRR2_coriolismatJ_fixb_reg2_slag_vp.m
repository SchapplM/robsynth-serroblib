% Calculate inertial parameters regressor of coriolis matrix for
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
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRR2_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:15
% EndTime: 2019-12-31 16:48:16
% DurationCPUTime: 0.43s
% Computational Cost: add. (304->54), mult. (647->65), div. (0->0), fcn. (510->6), ass. (0->47)
t34 = sin(qJ(3));
t38 = cos(pkin(7)) * pkin(1) + pkin(2);
t51 = cos(qJ(3));
t52 = sin(pkin(7)) * pkin(1);
t20 = t34 * t52 - t51 * t38;
t61 = t20 / 0.2e1;
t18 = -pkin(3) + t20;
t60 = -pkin(3) / 0.2e1 + t18 / 0.2e1;
t59 = t61 + t60;
t35 = cos(qJ(4));
t43 = qJD(1) + qJD(3);
t58 = t35 * t43;
t33 = sin(qJ(4));
t30 = t33 ^ 2;
t31 = t35 ^ 2;
t25 = t31 - t30;
t57 = t43 * t25;
t50 = pkin(3) * qJD(3);
t21 = t34 * t38 + t51 * t52;
t19 = pkin(6) + t21;
t8 = (-t30 - t31) * t20;
t1 = t18 * t21 + t19 * t8;
t49 = t1 * qJD(1);
t48 = t8 * qJD(1);
t47 = qJD(1) * t18;
t46 = t20 * qJD(1);
t45 = t21 * qJD(1);
t17 = t21 * qJD(3);
t44 = t33 * qJD(4);
t29 = t35 * qJD(4);
t42 = t33 * t47;
t41 = t35 * t47;
t40 = t33 * t45;
t39 = t61 - t60;
t3 = t39 * t33;
t37 = t3 * qJD(1) + t33 * t50;
t4 = t39 * t35;
t36 = t4 * qJD(1) + t35 * t50;
t26 = t33 * t29;
t24 = t25 * qJD(4);
t22 = t33 * t58;
t16 = t20 * qJD(3);
t13 = t33 * t17;
t7 = t8 * qJD(3);
t6 = t59 * t35;
t5 = t59 * t33;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t16, 0, 0, t26, t24, 0, -t26, 0, 0, -t35 * t17 + t18 * t44, t18 * t29 + t13, t7, t1 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17 - t45, t16 + t46, 0, 0, t26, t24, 0, -t26, 0, 0, t5 * qJD(4) - t21 * t58, t6 * qJD(4) + t13 + t40, t7 + t48, t49 + (-t21 * pkin(3) + pkin(6) * t8) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t57, t29, -t22, -t44, 0, t5 * qJD(3) - t19 * t29 + t42, t6 * qJD(3) + t19 * t44 + t41, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t29, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t46, 0, 0, t26, t24, 0, -t26, 0, 0, -t3 * qJD(4) + t35 * t45, -t4 * qJD(4) - t40, -t48, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t24, 0, -t26, 0, 0, -pkin(3) * t44, -pkin(3) * t29, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t57, t29, -t22, -t44, 0, -pkin(6) * t29 - t37, pkin(6) * t44 - t36, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t57, 0, t22, 0, 0, t3 * qJD(3) - t42, t4 * qJD(3) - t41, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t57, 0, t22, 0, 0, t37, t36, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
