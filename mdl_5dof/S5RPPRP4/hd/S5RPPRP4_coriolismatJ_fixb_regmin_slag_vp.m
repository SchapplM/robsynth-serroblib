% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRP4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 17:13:18
% EndTime: 2021-01-15 17:13:20
% DurationCPUTime: 0.39s
% Computational Cost: add. (440->74), mult. (667->94), div. (0->0), fcn. (549->4), ass. (0->65)
t46 = sin(qJ(4));
t42 = t46 ^ 2;
t47 = cos(qJ(4));
t43 = t47 ^ 2;
t33 = t42 + t43;
t44 = sin(pkin(7));
t45 = cos(pkin(7));
t48 = -pkin(1) - pkin(2);
t79 = t45 * qJ(2) + t44 * t48;
t28 = -pkin(6) + t79;
t83 = -qJ(5) + t28;
t81 = pkin(4) * t46;
t80 = t47 * pkin(4);
t51 = -t44 * qJ(2) + t45 * t48;
t27 = pkin(3) - t51;
t21 = t27 + t80;
t1 = t21 * t81;
t78 = t1 * qJD(1);
t17 = t83 * t46;
t18 = t83 * t47;
t8 = t17 * t46 + t18 * t47;
t2 = t21 * t44 + t8 * t45;
t77 = t2 * qJD(1);
t3 = t45 * t81;
t76 = t3 * qJD(1);
t75 = t8 * qJD(1);
t9 = -t51 * t44 + t79 * t45;
t74 = t9 * qJD(1);
t73 = qJD(2) * t44;
t72 = qJD(2) * t45;
t12 = (t21 + t80) * t46;
t71 = t12 * qJD(1);
t15 = t42 * pkin(4) - t21 * t47;
t70 = t15 * qJD(1);
t69 = t18 * qJD(4);
t19 = (0.1e1 / 0.2e1 + t43 / 0.2e1 + t42 / 0.2e1) * t44;
t68 = t19 * qJD(1);
t26 = t33 * t45;
t67 = t26 * qJD(1);
t66 = t33 * qJD(1);
t34 = t43 - t42;
t65 = t34 * qJD(1);
t64 = t46 * qJD(1);
t63 = t46 * qJD(4);
t62 = t47 * qJD(1);
t61 = t47 * qJD(4);
t60 = qJ(2) * qJD(1);
t59 = pkin(4) * t63;
t58 = pkin(4) * t64;
t57 = t46 * t73;
t56 = t45 * t62;
t55 = t44 * t64;
t54 = t45 * t64;
t53 = t46 * t62;
t52 = t44 * t61;
t50 = qJD(5) - t72;
t49 = qJD(1) * t27 - t72;
t32 = t44 * t62;
t31 = t47 * t73;
t25 = t45 * t63 - t32;
t24 = t44 * t63 + t56;
t23 = -t52 + t54;
t22 = t45 * t61 + t55;
t20 = (0.1e1 / 0.2e1 - t33 / 0.2e1) * t44;
t4 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), t9 * qJD(2), t46 * t61, t34 * qJD(4), 0, 0, 0, -t27 * t63 + t31, -t27 * t61 - t57, -t12 * qJD(4) + t31, t15 * qJD(4) - t57, -t26 * qJD(2) + t33 * qJD(5), t2 * qJD(2) - t1 * qJD(4) - t8 * qJD(5); 0, 0, 0, 0, qJD(1), t60, t74, 0, 0, 0, 0, 0, t32, -t55, t32, -t55, -t67, t77 + t20 * qJD(5) + (-0.1e1 + t33) * t44 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t53, t65, -t61, t63, 0, -t27 * t64 - t28 * t61, -t27 * t62 + t28 * t63, -t69 - t71, t17 * qJD(4) + t70, pkin(4) * t61, -pkin(4) * t69 - t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t20 * qJD(2) - t75; 0, 0, 0, 0, -qJD(1), -t60, -t74, 0, 0, 0, 0, 0, t25, t22, t25, t22, t67, t3 * qJD(4) - t19 * qJD(5) - t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t24, t23, t24, 0, -pkin(4) * t52 + t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t61, -t63, -t61, 0, -t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t53, -t65, 0, 0, 0, t49 * t46, t49 * t47, t50 * t46 + t71, t50 * t47 - t70, 0, -t3 * qJD(2) + qJD(5) * t81 + t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t56, -t54, -t56, 0, -t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t62, 0, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t61, -t66, t19 * qJD(2) - t59 + t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t62, 0, -t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t4;
