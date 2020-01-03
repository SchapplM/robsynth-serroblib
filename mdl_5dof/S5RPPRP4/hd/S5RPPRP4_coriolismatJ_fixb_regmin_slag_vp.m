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
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2019-12-31 17:52:23
% EndTime: 2019-12-31 17:52:23
% DurationCPUTime: 0.33s
% Computational Cost: add. (383->58), mult. (579->85), div. (0->0), fcn. (472->4), ass. (0->53)
t37 = sin(qJ(4));
t33 = t37 ^ 2;
t38 = cos(qJ(4));
t34 = t38 ^ 2;
t25 = t33 + t34;
t35 = sin(pkin(7));
t36 = cos(pkin(7));
t39 = -pkin(1) - pkin(2);
t68 = t36 * qJ(2) + t35 * t39;
t22 = -pkin(6) + t68;
t71 = qJ(5) - t22;
t69 = pkin(4) * t37;
t41 = -t35 * qJ(2) + t36 * t39;
t21 = pkin(3) - t41;
t19 = t38 * pkin(4) + t21;
t1 = t19 * t69;
t67 = t1 * qJD(1);
t16 = t71 * t38;
t8 = -t16 * t38 - t71 * t33;
t2 = t19 * t35 + t36 * t8;
t66 = t2 * qJD(1);
t3 = t36 * t69;
t65 = t3 * qJD(1);
t64 = t8 * qJD(1);
t9 = -t41 * t35 + t68 * t36;
t63 = t9 * qJD(1);
t62 = qJD(1) * t37;
t61 = qJD(1) * t38;
t17 = (0.1e1 / 0.2e1 + t34 / 0.2e1 + t33 / 0.2e1) * t35;
t60 = t17 * qJD(1);
t20 = t25 * t36;
t59 = t20 * qJD(1);
t58 = t25 * qJD(1);
t26 = t34 - t33;
t57 = t26 * qJD(1);
t56 = t35 * qJD(1);
t55 = t35 * qJD(2);
t54 = t36 * qJD(1);
t53 = t36 * qJD(2);
t52 = t37 * qJD(4);
t51 = t38 * qJD(4);
t50 = qJ(2) * qJD(1);
t49 = pkin(4) * t52;
t48 = pkin(4) * t62;
t47 = t38 * t54;
t46 = t37 * t56;
t45 = t37 * t54;
t44 = t37 * t61;
t43 = t38 * t56;
t42 = t35 * t51;
t40 = qJD(1) * t21 - t53;
t18 = (0.1e1 / 0.2e1 - t25 / 0.2e1) * t35;
t4 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), t55, t53, t9 * qJD(2), t37 * t51, t26 * qJD(4), 0, 0, 0, -t21 * t52 + t38 * t55, -t21 * t51 - t37 * t55, -t20 * qJD(2) + t25 * qJD(5), t2 * qJD(2) - t1 * qJD(4) - t8 * qJD(5); 0, 0, 0, 0, qJD(1), t50, t56, t54, t63, 0, 0, 0, 0, 0, t43, -t46, -t59, t66 + t18 * qJD(5) + (-0.1e1 + t25) * t35 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t57, -t51, t52, 0, -t21 * t62 - t22 * t51, -t21 * t61 + t22 * t52, pkin(4) * t51, t16 * pkin(4) * qJD(4) - t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t18 * qJD(2) - t64; 0, 0, 0, 0, -qJD(1), -t50, -t56, -t54, -t63, 0, 0, 0, 0, 0, t36 * t52 - t43, t36 * t51 + t46, t59, t3 * qJD(4) - t17 * qJD(5) - t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42 + t45, t35 * t52 + t47, 0, -pkin(4) * t42 + t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t51, 0, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t57, 0, 0, 0, t40 * t37, t40 * t38, 0, -t3 * qJD(2) + qJD(5) * t69 + t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t47, 0, -t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t17 * qJD(2) - t49 + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t4;
