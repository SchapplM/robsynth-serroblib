% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:40
% EndTime: 2019-12-31 16:51:41
% DurationCPUTime: 0.32s
% Computational Cost: add. (451->60), mult. (812->97), div. (0->0), fcn. (339->4), ass. (0->59)
t24 = sin(qJ(3));
t51 = t24 * qJD(2);
t26 = cos(qJ(3));
t53 = qJD(3) * t26;
t27 = -pkin(1) - pkin(2);
t15 = qJD(1) * t27 + qJD(2);
t59 = t24 * t15;
t2 = (qJ(2) * t53 + t51) * qJD(1) + qJD(3) * t59;
t48 = qJD(1) - qJD(3);
t50 = qJD(1) * qJ(2);
t8 = t26 * t50 + t59;
t67 = -t8 * t48 - t2;
t28 = qJD(4) ^ 2;
t69 = t48 ^ 2;
t70 = t24 * (t28 + t69);
t63 = t48 * pkin(3);
t40 = t24 * t50;
t7 = t15 * t26 - t40;
t3 = -t7 + t63;
t68 = t48 * t3;
t52 = qJD(4) * t48;
t49 = qJD(1) * qJD(2);
t1 = -qJD(3) * t40 + t15 * t53 + t26 * t49;
t65 = t48 * t7 + t1;
t23 = sin(qJ(4));
t21 = t23 ^ 2;
t25 = cos(qJ(4));
t22 = t25 ^ 2;
t54 = t21 + t22;
t42 = t54 * t1;
t56 = t26 * qJ(2) + t24 * t27;
t6 = qJD(3) * t56 + t51;
t64 = t48 * t6 + t2;
t60 = t23 * t25;
t58 = t28 * t23;
t57 = t28 * t25;
t55 = t21 - t22;
t47 = t69 * t60;
t45 = 0.2e1 * t49;
t44 = -t1 + t68;
t43 = t54 * t7;
t37 = -qJ(2) * t24 + t26 * t27;
t5 = t26 * qJD(2) + qJD(3) * t37;
t41 = t54 * t5;
t38 = t52 * t60;
t36 = pkin(6) * t28 - t67;
t12 = -pkin(6) + t56;
t35 = t12 * t28 - t64;
t34 = qJD(4) * (t3 + t7 + t63);
t11 = pkin(3) - t37;
t33 = qJD(4) * (-t11 * t48 - t3 - t5);
t32 = 0.2e1 * t26 * t52;
t30 = t48 * t54;
t29 = qJD(1) ^ 2;
t14 = 0.2e1 * t38;
t13 = -0.2e1 * t38;
t9 = t55 * t52;
t4 = -pkin(6) * t48 + t8;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, qJ(2) * t45, 0, 0, 0, 0, 0, 0, t64, t48 * t5 + t1, 0, t1 * t56 - t2 * t37 + t5 * t8 - t6 * t7, t14, -0.2e1 * t9, -t57, t13, t58, 0, t23 * t33 - t25 * t35, t23 * t35 + t25 * t33, -t41 * t48 - t42, t2 * t11 + t12 * t42 + t3 * t6 + t4 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t29 * qJ(2), 0, 0, 0, 0, 0, 0, -t24 * t69, -t26 * t69, 0, t65 * t24 + t26 * t67, 0, 0, 0, 0, 0, 0, t23 * t32 - t25 * t70, t23 * t70 + t25 * t32, t30 * t26 * t48, (t42 - t68) * t24 + (-t30 * t4 - t2) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t65, 0, 0, t13, 0.2e1 * t9, t57, t14, -t58, 0, t23 * t34 - t25 * t36, t23 * t36 + t25 * t34, t43 * t48 + t42, -t2 * pkin(3) + pkin(6) * t42 - t3 * t8 - t4 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t55 * t69, 0, t47, 0, 0, t44 * t23, t44 * t25, 0, 0;];
tauc_reg = t10;
