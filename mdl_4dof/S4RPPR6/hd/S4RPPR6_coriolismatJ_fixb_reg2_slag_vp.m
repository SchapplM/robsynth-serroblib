% Calculate inertial parameters regressor of coriolis matrix for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPPR6_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:46
% EndTime: 2019-12-31 16:40:47
% DurationCPUTime: 0.36s
% Computational Cost: add. (351->52), mult. (764->73), div. (0->0), fcn. (778->4), ass. (0->49)
t40 = sin(pkin(6));
t41 = cos(pkin(6));
t65 = sin(qJ(4));
t66 = cos(qJ(4));
t23 = t40 * t65 + t41 * t66;
t69 = t23 ^ 2;
t25 = t40 * t66 - t41 * t65;
t68 = t25 ^ 2;
t67 = -t40 / 0.2e1;
t64 = -pkin(5) + qJ(2);
t38 = t40 ^ 2;
t35 = t41 ^ 2 + t38;
t31 = t64 * t41;
t45 = t64 * t40;
t14 = t65 * t31 - t66 * t45;
t15 = t66 * t31 + t65 * t45;
t3 = -t14 * t25 + t15 * t23;
t63 = t3 * qJD(1);
t6 = -t68 - t69;
t62 = t6 * qJD(1);
t7 = -t68 + t69;
t61 = t7 * qJD(1);
t42 = -t23 * t65 / 0.2e1 - t25 * t66 / 0.2e1;
t9 = t67 + t42;
t60 = t9 * qJD(1);
t59 = t23 * qJD(1);
t58 = t23 * qJD(4);
t57 = t25 * qJD(1);
t21 = t25 * qJD(4);
t30 = t35 * qJ(2);
t56 = t30 * qJD(1);
t55 = t35 * qJD(1);
t54 = t38 * qJD(1);
t53 = t40 * qJD(1);
t52 = t40 * qJD(3);
t44 = t40 * qJ(3) + pkin(1);
t20 = (pkin(2) + pkin(3)) * t41 + t44;
t51 = t20 * t53;
t50 = t23 * t57;
t49 = t23 * t21;
t48 = t23 * t53;
t47 = t25 * t53;
t46 = t41 * t53;
t43 = qJD(1) * t20 - qJD(2);
t32 = t35 * qJD(2);
t29 = -t41 * pkin(2) - t44;
t26 = t30 * qJD(2);
t8 = t67 - t42;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t26, 0, 0, 0, 0, 0, 0, t41 * t52, t32, t38 * qJD(3), -t29 * t52 + t26, -t49, t7 * qJD(4), 0, t49, 0, 0, t20 * t21 + t23 * t52, -t20 * t58 + t25 * t52, t6 * qJD(2), t3 * qJD(2) + t20 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t56, 0, 0, 0, 0, 0, 0, 0, t55, 0, t56, 0, 0, 0, 0, 0, 0, 0, 0, t62, t8 * qJD(3) + t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, t54, -t29 * t53, 0, 0, 0, 0, 0, 0, t48, t47, 0, t8 * qJD(2) + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t61, -t58, t50, -t21, 0, -t15 * qJD(4) + t20 * t57, t14 * qJD(4) - t20 * t59, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t56, 0, 0, 0, 0, 0, 0, 0, -t55, 0, -t56 - t52, 0, 0, 0, 0, 0, 0, -t21, t58, -t62, t9 * qJD(3) - t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t59, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, 0, -t54, (qJD(1) * t29 + qJD(2)) * t40, 0, 0, 0, 0, 0, 0, -t48, -t47, 0, -t9 * qJD(2) - t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65 * qJD(4), -t66 * qJD(4), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t61, 0, -t50, 0, 0, -t43 * t25, t43 * t23, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t59, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
