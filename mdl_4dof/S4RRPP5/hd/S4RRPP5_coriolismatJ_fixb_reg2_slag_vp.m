% Calculate inertial parameters regressor of coriolis matrix for
% S4RRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPP5_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:40
% EndTime: 2019-12-31 17:00:41
% DurationCPUTime: 0.47s
% Computational Cost: add. (217->74), mult. (411->74), div. (0->0), fcn. (316->2), ass. (0->56)
t34 = cos(qJ(2));
t25 = t34 * qJD(3);
t33 = sin(qJ(2));
t55 = t33 * qJ(3);
t36 = -t34 * pkin(2) - t55;
t62 = t36 * qJD(2) + t25;
t32 = pkin(2) + qJ(4);
t61 = -t32 * t34 - t55;
t60 = pkin(3) + pkin(5);
t8 = -pkin(1) + t61;
t12 = t33 * pkin(2) - t34 * qJ(3);
t9 = t33 * qJ(4) + t12;
t1 = t8 * t9;
t58 = t1 * qJD(1);
t2 = t9 * t33 + t8 * t34;
t57 = t2 * qJD(1);
t3 = -t8 * t33 + t9 * t34;
t56 = t3 * qJD(1);
t11 = -pkin(1) + t36;
t4 = t11 * t34 + t12 * t33;
t54 = t4 * qJD(1);
t5 = -t11 * t33 + t12 * t34;
t53 = t5 * qJD(1);
t52 = qJD(1) * t33;
t51 = qJD(1) * t34;
t50 = qJD(3) * t33;
t49 = qJD(4) * t34;
t17 = t60 * t34;
t48 = t17 * qJD(2);
t30 = t33 ^ 2;
t31 = t34 ^ 2;
t18 = t31 - t30;
t47 = t18 * qJD(1);
t13 = t18 * qJD(2);
t46 = t31 * qJD(1);
t45 = t32 * qJD(2);
t24 = t33 * qJD(2);
t26 = t34 * qJD(2);
t44 = pkin(1) * t52;
t43 = pkin(1) * t51;
t42 = pkin(5) * t24;
t41 = t8 * t52;
t40 = t8 * t51;
t39 = t11 * t12 * qJD(1);
t38 = t11 * t52;
t37 = t33 * t25;
t29 = qJ(3) * qJD(3);
t28 = qJD(2) * qJ(3);
t23 = t30 * qJD(1);
t22 = t30 * qJD(3);
t21 = pkin(5) * t26;
t20 = t33 * t26;
t19 = t33 * t51;
t16 = t60 * t33;
t10 = t16 * qJD(2);
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t13, 0, -t20, 0, 0, -pkin(1) * t24, -pkin(1) * t26, 0, 0, 0, 0, 0, t20, t13, -t20, 0, t5 * qJD(2) - t37, -t4 * qJD(2) + t22, (qJD(2) * t12 - t50) * t11, 0, 0, 0, -t20, -t13, t20, 0, -t2 * qJD(2) + t33 * t49 + t22, -t3 * qJD(2) + t31 * qJD(4) + t37, t1 * qJD(2) + (-t49 - t50) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t47, t26, -t19, -t24, 0, -t21 - t44, t42 - t43, 0, 0, 0, -t26, t24, t19, t47, -t19, t62, t21 + t53, -t42 - t54, t62 * pkin(5) + t39, 0, t24, t26, -t19, -t47, t19, t61 * qJD(2) - qJD(4) * t33 + t25, -t10 - t57, -t48 - t56, t58 + (-t16 * qJ(3) - t17 * t32) * qJD(2) + t17 * qJD(3) - t16 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t19, t23, t21 - t38, 0, 0, 0, 0, 0, 0, t26, t23, t19, -t41 + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t19, t46, -t10 - t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t47, 0, t19, 0, 0, t44, t43, 0, 0, 0, 0, 0, -t19, -t47, t19, 0, -t53, t54, -t39, 0, 0, 0, t19, t47, -t19, 0, t57, t56, -t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t29, 0, 0, 0, 0, 0, 0, 0, qJD(3), qJD(4), t32 * qJD(4) + t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t28, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t23, t38, 0, 0, 0, 0, 0, 0, 0, -t23, -t19, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t28, 0, 0, 0, 0, 0, 0, 0, -qJD(2), 0, -qJD(4) - t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t46, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), qJD(3) - t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
