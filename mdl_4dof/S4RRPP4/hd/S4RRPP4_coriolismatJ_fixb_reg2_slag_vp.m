% Calculate inertial parameters regressor of coriolis matrix for
% S4RRPP4
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
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPP4_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:59:21
% EndTime: 2019-12-31 16:59:22
% DurationCPUTime: 0.49s
% Computational Cost: add. (240->83), mult. (461->86), div. (0->0), fcn. (347->2), ass. (0->59)
t38 = sin(qJ(2));
t33 = t38 * qJ(3);
t39 = cos(qJ(2));
t65 = pkin(2) + pkin(3);
t68 = t65 * t39 + t33;
t41 = -t39 * pkin(2) - t33;
t49 = t39 * qJD(3);
t67 = t41 * qJD(2) + t49;
t66 = pkin(5) - qJ(4);
t10 = pkin(1) + t68;
t60 = t39 * qJ(3);
t11 = -t65 * t38 + t60;
t1 = t10 * t11;
t63 = t1 * qJD(1);
t2 = t10 * t39 + t11 * t38;
t62 = t2 * qJD(1);
t3 = -t10 * t38 + t11 * t39;
t61 = t3 * qJD(1);
t13 = -pkin(1) + t41;
t15 = t38 * pkin(2) - t60;
t4 = t13 * t39 + t15 * t38;
t59 = t4 * qJD(1);
t5 = -t13 * t38 + t15 * t39;
t58 = t5 * qJD(1);
t14 = t66 * t38;
t16 = t66 * t39;
t7 = t14 * t38 + t16 * t39;
t57 = t7 * qJD(1);
t48 = -pkin(2) / 0.2e1 - pkin(3) / 0.2e1;
t8 = t60 + (-t65 / 0.2e1 + t48) * t38;
t56 = t8 * qJD(1);
t55 = t16 * qJD(2);
t36 = t38 ^ 2;
t37 = t39 ^ 2;
t20 = t37 + t36;
t54 = t20 * qJD(1);
t21 = t37 - t36;
t53 = t21 * qJD(1);
t18 = t21 * qJD(2);
t30 = t38 * qJD(1);
t52 = t38 * qJD(2);
t51 = t38 * qJD(3);
t50 = t39 * qJD(1);
t31 = t39 * qJD(2);
t47 = pkin(1) * t30;
t46 = pkin(1) * t50;
t45 = pkin(5) * t52;
t44 = pkin(5) * t31;
t43 = t13 * t15 * qJD(1);
t42 = t13 * t30;
t35 = qJ(3) * qJD(3);
t34 = qJD(2) * qJ(3);
t29 = t36 * qJD(1);
t28 = t36 * qJD(3);
t24 = t38 * t31;
t23 = t38 * t50;
t22 = t38 * t49;
t9 = (t65 / 0.2e1 + t48) * t38;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t18, 0, -t24, 0, 0, -pkin(1) * t52, -pkin(1) * t31, 0, 0, t24, 0, -t18, 0, 0, -t24, -t5 * qJD(2) + t22, 0, -t4 * qJD(2) + t28, (qJD(2) * t15 - t51) * t13, t24, -t18, 0, -t24, 0, 0, t3 * qJD(2) + t22, t2 * qJD(2) + t28, t20 * qJD(4), t1 * qJD(2) - t7 * qJD(4) + t10 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t53, t31, -t23, -t52, 0, -t44 - t47, t45 - t46, 0, 0, t23, t31, -t53, 0, t52, -t23, -t44 - t58, t67, -t45 - t59, t67 * pkin(5) + t43, t23, -t53, -t31, -t23, -t52, 0, -t55 + t61, -t14 * qJD(2) + t62, t68 * qJD(2) - t49, t63 + (-t14 * qJ(3) - t16 * t65) * qJD(2) + t16 * qJD(3) + t9 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t31, t29, -t42 + t44, 0, 0, 0, 0, 0, 0, t23, t29, -t31, t10 * t30 + t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t9 * qJD(2) - t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t53, 0, t23, 0, 0, t47, t46, 0, 0, -t23, 0, t53, 0, 0, t23, t58, 0, t59, -t43, -t23, t53, 0, t23, 0, 0, t38 * qJD(4) - t61, -t39 * qJD(4) - t62, 0, -t8 * qJD(4) - t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t35, 0, 0, 0, 0, 0, 0, 0, qJD(3), 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t34, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t50, 0, -t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, -t29, t42, 0, 0, 0, 0, 0, 0, -t23, -t29, 0, (-qJD(1) * t10 - qJD(4)) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t34, 0, 0, 0, 0, 0, 0, 0, -qJD(2), 0, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, t31, -t54, t8 * qJD(2) + t51 + t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t50, 0, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
