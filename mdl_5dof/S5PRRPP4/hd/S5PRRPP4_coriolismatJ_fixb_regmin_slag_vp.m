% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRPP4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:41:14
% EndTime: 2019-12-31 17:41:16
% DurationCPUTime: 0.38s
% Computational Cost: add. (254->76), mult. (479->88), div. (0->0), fcn. (351->2), ass. (0->61)
t41 = sin(qJ(3));
t35 = t41 * qJ(4);
t42 = cos(qJ(3));
t68 = pkin(3) + pkin(4);
t73 = t68 * t42 + t35;
t44 = -t42 * pkin(3) - t35;
t51 = t42 * qJD(4);
t72 = t44 * qJD(3) + t51;
t71 = pkin(6) - qJ(5);
t36 = t41 * pkin(3);
t61 = t42 * qJ(4);
t19 = -t36 + t61;
t30 = t41 * qJD(4);
t70 = t19 * qJD(3) + t30;
t69 = -t36 / 0.2e1;
t67 = t41 * pkin(4);
t65 = t68 * t41;
t12 = pkin(2) + t73;
t15 = t19 - t67;
t1 = t12 * t15;
t64 = t1 * qJD(2);
t2 = t12 * t42 + t15 * t41;
t63 = t2 * qJD(2);
t3 = -t12 * t41 + t15 * t42;
t62 = t3 * qJD(2);
t17 = -pkin(2) + t44;
t6 = t17 * t42 - t19 * t41;
t60 = t6 * qJD(2);
t7 = -t17 * t41 - t19 * t42;
t59 = t7 * qJD(2);
t18 = t71 * t41;
t20 = t71 * t42;
t9 = t18 * t41 + t20 * t42;
t58 = t9 * qJD(2);
t10 = t61 + t69 + (-pkin(4) / 0.2e1 - t68 / 0.2e1) * t41;
t57 = t10 * qJD(2);
t55 = t20 * qJD(3);
t39 = t41 ^ 2;
t40 = t42 ^ 2;
t21 = t39 + t40;
t54 = t21 * qJD(2);
t22 = t40 - t39;
t53 = t22 * qJD(2);
t32 = t41 * qJD(2);
t31 = t41 * qJD(3);
t52 = t42 * qJD(2);
t33 = t42 * qJD(3);
t50 = pkin(2) * t32;
t49 = pkin(2) * t52;
t48 = pkin(6) * t31;
t47 = pkin(6) * t33;
t46 = t17 * t19 * qJD(2);
t45 = t17 * t32;
t38 = qJ(4) * qJD(4);
t37 = qJD(3) * qJ(4);
t29 = t39 * qJD(2);
t28 = t39 * qJD(4);
t24 = t41 * t52;
t23 = t41 * t51;
t11 = t65 / 0.2e1 + t69 - t67 / 0.2e1;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t33, -t31, 0, t33, t70, -t31, t33, 0, (t61 - t65) * qJD(3) + t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t41 * t33, t22 * qJD(3), 0, 0, 0, -pkin(2) * t31, -pkin(2) * t33, -t7 * qJD(3) + t23, 0, -t6 * qJD(3) + t28, -t70 * t17, t3 * qJD(3) + t23, t2 * qJD(3) + t28, t21 * qJD(5), t1 * qJD(3) - t9 * qJD(5) + t12 * t30; 0, 0, 0, 0, t24, t53, t33, -t31, 0, -t47 - t50, t48 - t49, -t47 - t59, t72, -t48 - t60, t72 * pkin(6) - t46, -t55 + t62, -t18 * qJD(3) + t63, t73 * qJD(3) - t51, t64 + (-t18 * qJ(4) - t20 * t68) * qJD(3) + t20 * qJD(4) + t11 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t33, t29, -t45 + t47, t24, t29, -t33, t12 * t32 + t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t11 * qJD(3) - t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t24, -t53, 0, 0, 0, t50, t49, t59, 0, t60, t46, t41 * qJD(5) - t62, -t42 * qJD(5) - t63, 0, -t10 * qJD(5) - t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t38, 0, qJD(4), 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t37, 0, qJD(3), 0, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t52, 0, -t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, -t29, t45, -t24, -t29, 0, (-qJD(2) * t12 - qJD(5)) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t37, 0, -qJD(3), 0, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t33, -t54, t10 * qJD(3) + t30 + t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t52, 0, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t4;
