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
% MMD_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:39
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 14:39:15
% EndTime: 2021-01-15 14:39:18
% DurationCPUTime: 0.56s
% Computational Cost: add. (365->106), mult. (1023->217), div. (0->0), fcn. (696->4), ass. (0->76)
t32 = sin(qJ(2));
t85 = -0.4e1 * t32;
t34 = cos(qJ(2));
t41 = -t34 * pkin(2) - t32 * pkin(6);
t18 = -pkin(1) + t41;
t33 = cos(qJ(3));
t73 = t33 * t34;
t23 = pkin(5) * t73;
t31 = sin(qJ(3));
t69 = t31 * t18 + t23;
t84 = qJD(4) * t32 + (pkin(5) * qJD(3) + qJ(4) * qJD(2)) * t34;
t83 = 0.2e1 * qJD(3);
t82 = pkin(5) * t31;
t81 = t33 * pkin(3);
t26 = qJD(3) * t33;
t60 = t34 * qJD(2);
t11 = t32 * t26 + t31 * t60;
t55 = pkin(5) * t60;
t7 = t11 * pkin(3) + t55;
t80 = t7 * t31;
t79 = t7 * t33;
t17 = (pkin(3) * t31 + pkin(5)) * t32;
t78 = t17 * t31;
t72 = qJ(4) + pkin(6);
t19 = t72 * t31;
t77 = t19 * t32;
t20 = t72 * t33;
t76 = t20 * t32;
t75 = t31 * t34;
t74 = t32 * t33;
t40 = pkin(2) * t32 - pkin(6) * t34;
t37 = t40 * qJD(2);
t71 = -t18 * t26 - t31 * t37;
t61 = t32 * qJD(2);
t50 = t31 * t61;
t70 = pkin(5) * t50 + t33 * t37;
t27 = t31 ^ 2;
t29 = t33 ^ 2;
t68 = t27 - t29;
t28 = t32 ^ 2;
t67 = -t34 ^ 2 + t28;
t66 = qJ(4) * t32;
t65 = qJD(2) * t33;
t64 = qJD(3) * t31;
t63 = qJD(3) * t34;
t59 = -0.2e1 * pkin(1) * qJD(2);
t58 = -0.2e1 * pkin(2) * qJD(3);
t57 = pkin(3) * t61;
t56 = pkin(3) * t64;
t54 = t32 * t64;
t53 = t31 * t63;
t52 = t33 * t63;
t51 = t17 * t64;
t49 = t31 * t26;
t48 = t32 * t60;
t47 = t33 * t60;
t24 = -pkin(2) - t81;
t46 = -t24 + t81;
t45 = t68 * qJD(3);
t44 = t67 * qJD(2);
t43 = 0.2e1 * t48;
t42 = t31 * t47;
t39 = pkin(3) * t27 + t24 * t33;
t10 = t47 - t54;
t36 = t33 * t61 + t53;
t35 = qJ(4) * t54 - t18 * t64 - t84 * t33 + t70;
t16 = t33 * t18;
t9 = -t31 * qJD(4) - t72 * t26;
t8 = -t33 * qJD(4) + t72 * t64;
t6 = -t31 * t66 + t69;
t5 = -t33 * t66 + t16 + (-pkin(3) - t82) * t34;
t4 = -t69 * qJD(3) + t70;
t3 = t36 * pkin(5) + t71;
t2 = (pkin(5) * qJD(2) + qJ(4) * qJD(3)) * t74 + t84 * t31 + t71;
t1 = t35 + t57;
t12 = [0, 0, 0, t43, -0.2e1 * t44, 0, 0, 0, t32 * t59, t34 * t59, -0.2e1 * t28 * t49 + 0.2e1 * t29 * t48, t68 * t28 * t83 + t42 * t85, 0.2e1 * t32 * t53 + 0.2e1 * t67 * t65, -0.2e1 * t31 * t44 + 0.2e1 * t32 * t52, -0.2e1 * t48, 0.2e1 * t16 * t61 - 0.2e1 * t4 * t34 + 0.2e1 * (t28 * t26 + t31 * t48) * pkin(5), -0.2e1 * t3 * t34 - 0.2e1 * t69 * t61 + 0.2e1 * (-t28 * t64 + t33 * t43) * pkin(5), 0.2e1 * (qJD(2) * t78 - t1) * t34 + 0.2e1 * (qJD(2) * t5 + t17 * t26 + t80) * t32, 0.2e1 * (t17 * t65 - t2) * t34 + 0.2e1 * (-qJD(2) * t6 - t51 + t79) * t32, 0.2e1 * (-t31 * t6 - t33 * t5) * t60 + 0.2e1 * (-t1 * t33 + t2 * t31 + (t31 * t5 - t33 * t6) * qJD(3)) * t32, 0.2e1 * t5 * t1 + 0.2e1 * t17 * t7 - 0.2e1 * t6 * t2; 0, 0, 0, 0, 0, t60, -t61, 0, -t55, pkin(5) * t61, -t32 * t45 + t42, t49 * t85 - t68 * t60, t50 - t52, t36, 0, (pkin(6) * t73 + (-pkin(2) * t33 + t82) * t32) * qJD(3) + (t41 * t31 - t23) * qJD(2), (pkin(5) * t74 + t40 * t31) * qJD(3) + (pkin(5) * t75 + t41 * t33) * qJD(2), -t79 - t9 * t34 + (t24 * t75 - t77) * qJD(2) + (t39 * t32 + t78) * qJD(3), t80 - t8 * t34 + (t24 * t73 - t76) * qJD(2) + (t46 * t32 * t31 + t17 * t33) * qJD(3), (t19 * t60 - t32 * t9 - t2 + (-t5 - t76) * qJD(3)) * t33 + (-t20 * t60 + t32 * t8 - t1 + (-t6 - t77) * qJD(3)) * t31, pkin(3) * t51 - t1 * t19 - t2 * t20 + t7 * t24 + t5 * t9 - t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t49, -0.2e1 * t45, 0, 0, 0, t31 * t58, t33 * t58, -0.2e1 * t46 * t64, t39 * t83, -0.2e1 * t9 * t31 - 0.2e1 * t8 * t33 + 0.2e1 * (t19 * t33 - t20 * t31) * qJD(3), -0.2e1 * t19 * t9 - 0.2e1 * t20 * t8 + 0.2e1 * t24 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t11, t61, t4, t3, t35 + 0.2e1 * t57, t2, -t10 * pkin(3), t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t64, 0, -pkin(6) * t26, pkin(6) * t64, t9, t8, -pkin(3) * t26, t9 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t10, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t26, 0, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
