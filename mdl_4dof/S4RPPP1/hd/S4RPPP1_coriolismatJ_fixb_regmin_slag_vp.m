% Calculate minimal parameter regressor of coriolis matrix for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x15]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPPP1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:26
% EndTime: 2019-03-08 18:26:27
% DurationCPUTime: 0.24s
% Computational Cost: add. (164->68), mult. (550->84), div. (0->0), fcn. (478->4), ass. (0->59)
t34 = sin(pkin(6));
t37 = cos(pkin(4));
t35 = sin(pkin(4));
t53 = qJD(2) * t35;
t31 = t35 ^ 2;
t36 = cos(pkin(6));
t67 = t31 * t36;
t69 = (-qJD(3) * t67 + t37 * t53) * t34;
t68 = pkin(1) * t37;
t66 = t34 * t35;
t65 = t35 * t36;
t30 = t34 ^ 2;
t32 = t36 ^ 2;
t15 = (t30 + t32) * t31;
t12 = t15 * qJD(2);
t47 = t37 * qJD(3);
t64 = t47 * t65 + t12;
t33 = t37 ^ 2;
t19 = t31 * t30 + t33;
t42 = t36 * t53;
t25 = t37 * t42;
t63 = t19 * qJD(3) + t25;
t61 = qJ(2) * t35;
t62 = t34 * t68 + t36 * t61;
t11 = -t37 * qJ(3) - t62;
t38 = pkin(3) * t65 - t11;
t41 = -pkin(1) * t36 - pkin(2);
t7 = (pkin(3) + qJ(2)) * t66 + (-qJ(4) + t41) * t37;
t1 = (t7 * t34 + t38 * t36) * t35;
t60 = t1 * qJD(1);
t40 = -qJ(3) * t34 - pkin(1);
t9 = ((-pkin(2) - qJ(4)) * t36 + t40) * t35;
t2 = t7 * t37 + t9 * t65;
t59 = t2 * qJD(1);
t3 = -t38 * t37 + t9 * t66;
t58 = t3 * qJD(1);
t44 = t34 * t61;
t4 = t11 * t65 - (t41 * t37 + t44) * t66;
t57 = t4 * qJD(1);
t5 = t11 * t37 + (-pkin(2) * t36 + t40) * t34 * t31;
t56 = t5 * qJD(1);
t6 = (t62 * t36 + (-t36 * t68 + t44) * t34) * t35;
t55 = t6 * qJD(1);
t54 = qJD(1) * t35;
t52 = qJD(3) * t34;
t51 = t15 * qJD(1);
t50 = t19 * qJD(1);
t20 = t31 * t32 + t33;
t49 = t20 * qJD(1);
t48 = t37 * qJD(1);
t46 = t37 * qJD(4);
t45 = t34 * t67;
t26 = t34 * t53;
t27 = t34 * t54;
t43 = t36 * t54;
t24 = t37 * t43;
t22 = t37 * t27;
t21 = qJD(1) * t45;
t8 = [0, 0, 0, -t37 * t26, -t25, t12, t6 * qJD(2), t64, t69, t63, -t4 * qJD(2) - t5 * qJD(3), -t46 * t66 + t64, qJD(4) * t45 + t63, t20 * qJD(4) - t69, t1 * qJD(2) - t3 * qJD(3) - t2 * qJD(4); 0, 0, 0, -t22, -t24, t51, t55, t51, t22, t24, -t57, t51, t24, -t22, t60; 0, 0, 0, 0, 0, 0, 0, t24, -t21, t50, -t56, t24, t50, t21, -t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t21, t49, -t59; 0, 0, 0, t22, t24, -t51, -t55, -t51, -t22, -t24, -t35 * t52 + t57, -t51, -t24, t22, -t60 + (-qJD(4) * t36 - t52) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, 0, 0, 0, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43; 0, 0, 0, 0, 0, 0, 0, -t24, t21, -t50, t26 + t56, -t24, -t50, -t21, t26 - t46 + t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t21, -t49, t42 + t47 + t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t8;
