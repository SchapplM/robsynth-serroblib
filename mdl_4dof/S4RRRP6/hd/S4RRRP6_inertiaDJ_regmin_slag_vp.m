% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRP6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:19:10
% EndTime: 2019-12-31 17:19:12
% DurationCPUTime: 0.40s
% Computational Cost: add. (264->87), mult. (743->184), div. (0->0), fcn. (503->4), ass. (0->63)
t28 = sin(qJ(2));
t69 = -0.4e1 * t28;
t30 = cos(qJ(2));
t35 = -t30 * pkin(2) - t28 * pkin(6);
t17 = -pkin(1) + t35;
t29 = cos(qJ(3));
t66 = t29 * t30;
t21 = pkin(5) * t66;
t27 = sin(qJ(3));
t62 = t27 * t17 + t21;
t25 = t29 ^ 2;
t61 = t27 ^ 2 - t25;
t39 = t61 * qJD(3);
t68 = pkin(5) * t27;
t67 = t28 * t29;
t65 = -qJ(4) - pkin(6);
t34 = pkin(2) * t28 - pkin(6) * t30;
t15 = t34 * qJD(2);
t53 = qJD(3) * t29;
t64 = -t27 * t15 - t17 * t53;
t57 = qJD(2) * t28;
t44 = t27 * t57;
t63 = pkin(5) * t44 + t29 * t15;
t24 = t28 ^ 2;
t60 = -t30 ^ 2 + t24;
t59 = qJ(4) * t28;
t58 = t29 * qJ(4);
t56 = qJD(2) * t29;
t55 = qJD(2) * t30;
t54 = qJD(3) * t27;
t52 = qJD(3) * t30;
t51 = t29 * qJD(4);
t50 = -0.2e1 * pkin(1) * qJD(2);
t49 = -0.2e1 * pkin(2) * qJD(3);
t48 = pkin(3) * t54;
t47 = pkin(5) * t55;
t46 = t27 * t52;
t45 = t29 * t52;
t43 = t27 * t53;
t42 = t28 * t55;
t41 = t29 * t55;
t40 = qJD(3) * t65;
t38 = t60 * qJD(2);
t37 = 0.2e1 * t42;
t36 = t27 * t41;
t33 = -t28 * t54 + t41;
t32 = t28 * t56 + t46;
t31 = t27 * t55 + t28 * t53;
t22 = -t29 * pkin(3) - pkin(2);
t19 = t65 * t29;
t18 = t65 * t27;
t16 = (pkin(3) * t27 + pkin(5)) * t28;
t14 = t29 * t17;
t9 = -t27 * qJD(4) + t29 * t40;
t8 = t27 * t40 + t51;
t7 = pkin(3) * t31 + t47;
t6 = -t27 * t59 + t62;
t5 = -t28 * t58 + t14 + (-pkin(3) - t68) * t30;
t4 = -t62 * qJD(3) + t63;
t3 = pkin(5) * t32 + t64;
t2 = (-pkin(5) * qJD(2) - qJ(4) * qJD(3)) * t67 + (-qJD(4) * t28 + (-pkin(5) * qJD(3) - qJ(4) * qJD(2)) * t30) * t27 - t64;
t1 = -t28 * t51 + (pkin(3) * t28 - t30 * t58) * qJD(2) + (-t21 + (-t17 + t59) * t27) * qJD(3) + t63;
t10 = [0, 0, 0, t37, -0.2e1 * t38, 0, 0, 0, t28 * t50, t30 * t50, -0.2e1 * t24 * t43 + 0.2e1 * t25 * t42, 0.2e1 * t24 * t39 + t36 * t69, 0.2e1 * t28 * t46 + 0.2e1 * t60 * t56, -0.2e1 * t27 * t38 + 0.2e1 * t28 * t45, -0.2e1 * t42, 0.2e1 * t14 * t57 - 0.2e1 * t4 * t30 + 0.2e1 * (t24 * t53 + t27 * t42) * pkin(5), -0.2e1 * t3 * t30 - 0.2e1 * t62 * t57 + 0.2e1 * (-t24 * t54 + t29 * t37) * pkin(5), 0.2e1 * (-t27 * t6 - t29 * t5) * t55 + 0.2e1 * (-t1 * t29 - t2 * t27 + (t27 * t5 - t29 * t6) * qJD(3)) * t28, 0.2e1 * t5 * t1 + 0.2e1 * t16 * t7 + 0.2e1 * t6 * t2; 0, 0, 0, 0, 0, t55, -t57, 0, -t47, pkin(5) * t57, -t28 * t39 + t36, t43 * t69 - t61 * t55, t44 - t45, t32, 0, (pkin(6) * t66 + (-pkin(2) * t29 + t68) * t28) * qJD(3) + (t27 * t35 - t21) * qJD(2), (pkin(5) * t67 + t27 * t34) * qJD(3) + (t29 * t35 + t30 * t68) * qJD(2), (-t18 * t55 - t28 * t9 + t2 + (t19 * t28 - t5) * qJD(3)) * t29 + (t19 * t55 - t28 * t8 - t1 + (t18 * t28 - t6) * qJD(3)) * t27, t1 * t18 + t16 * t48 - t2 * t19 + t7 * t22 + t5 * t9 + t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t43, -0.2e1 * t39, 0, 0, 0, t27 * t49, t29 * t49, -0.2e1 * t9 * t27 + 0.2e1 * t8 * t29 + 0.2e1 * (-t18 * t29 + t19 * t27) * qJD(3), 0.2e1 * t18 * t9 - 0.2e1 * t19 * t8 + 0.2e1 * t22 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t31, t57, t4, t3, -t33 * pkin(3), t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t54, 0, -pkin(6) * t53, pkin(6) * t54, -pkin(3) * t53, t9 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
