% Calculate minimal parameter regressor of coriolis matrix for
% S5PPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x16]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PPRRP2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:09:12
% EndTime: 2019-12-05 15:09:14
% DurationCPUTime: 0.29s
% Computational Cost: add. (271->52), mult. (659->80), div. (0->0), fcn. (674->6), ass. (0->52)
t38 = sin(qJ(4));
t71 = t38 * pkin(4);
t70 = cos(qJ(3));
t35 = t38 ^ 2;
t40 = cos(qJ(4));
t36 = t40 ^ 2;
t69 = -t35 - t36;
t37 = sin(pkin(8));
t39 = sin(qJ(3));
t66 = cos(pkin(8));
t26 = t39 * t37 - t70 * t66;
t27 = t70 * t37 + t39 * t66;
t1 = (0.1e1 + t69) * t27 * t26;
t68 = t1 * qJD(1);
t67 = t40 * qJ(5);
t65 = qJD(3) * t38;
t64 = qJD(3) * t40;
t63 = qJD(4) * t38;
t33 = qJD(4) * t40;
t48 = -t40 * pkin(4) - t38 * qJ(5);
t28 = -pkin(3) + t48;
t29 = t67 - t71;
t16 = t28 * t40 - t29 * t38;
t62 = t16 * qJD(3);
t17 = -t28 * t38 - t29 * t40;
t61 = t17 * qJD(3);
t60 = t27 * qJD(3);
t30 = t36 - t35;
t59 = t30 * qJD(3);
t58 = t35 * qJD(3);
t57 = t40 * qJD(5);
t56 = qJD(4) * qJ(5);
t55 = pkin(3) * t65;
t54 = pkin(3) * t64;
t53 = pkin(6) * t63;
t52 = pkin(6) * t33;
t51 = t28 * t65;
t50 = t69 * t26;
t47 = t29 * qJD(4) + t38 * qJD(5);
t46 = -t67 / 0.2e1 + t71 / 0.2e1;
t2 = (t29 / 0.2e1 + t46) * t26;
t45 = t28 * t29 * qJD(3) + t2 * qJD(1);
t13 = t26 * t38;
t44 = -t13 * qJD(3) + t27 * t33;
t15 = t26 * t40;
t43 = t15 * qJD(3) + t27 * t63;
t42 = t15 * qJD(4) + t38 * t60;
t41 = qJD(4) * t48 + t57;
t31 = t38 * t64;
t4 = t13 * qJD(4) - t40 * t60;
t3 = (-t29 / 0.2e1 + t46) * t26;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t60, t26 * qJD(3), 0, 0, 0, 0, 0, t4, t42, t4, qJD(3) * t50, -t42, t68 + (pkin(6) * t50 + t27 * t28) * qJD(3) + t3 * qJD(4) - t13 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t43, -t44, 0, -t43, t3 * qJD(3) + t27 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t33, -t63, 0, t33, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * qJD(4) - t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t38 * t33, t30 * qJD(4), 0, 0, 0, -pkin(3) * t63, -pkin(3) * t33, -t17 * qJD(4) + t38 * t57, 0, -t16 * qJD(4) + t35 * qJD(5), -t47 * t28; 0, 0, 0, 0, 0, t31, t59, t33, -t63, 0, -t52 - t55, t53 - t54, -t52 - t61, t41, -t53 - t62, pkin(6) * t41 - t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t33, t58, -t51 + t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t31, -t59, 0, 0, 0, t55, t54, t61, 0, t62, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, 0, -t58, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), -t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
