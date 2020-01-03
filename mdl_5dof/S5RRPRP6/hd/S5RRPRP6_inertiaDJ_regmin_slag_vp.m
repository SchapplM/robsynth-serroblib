% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:58:27
% EndTime: 2019-12-31 19:58:29
% DurationCPUTime: 0.56s
% Computational Cost: add. (897->117), mult. (2040->233), div. (0->0), fcn. (1800->6), ass. (0->77)
t45 = sin(pkin(8));
t46 = cos(pkin(8));
t48 = sin(qJ(2));
t50 = cos(qJ(2));
t32 = t45 * t48 - t46 * t50;
t33 = t45 * t50 + t46 * t48;
t66 = -t50 * pkin(2) - pkin(1);
t19 = t32 * pkin(3) - t33 * pkin(7) + t66;
t81 = -qJ(3) - pkin(6);
t36 = t81 * t48;
t37 = t81 * t50;
t22 = t45 * t36 - t46 * t37;
t49 = cos(qJ(4));
t20 = t49 * t22;
t47 = sin(qJ(4));
t80 = t47 * t19 + t20;
t27 = t33 * qJD(2);
t73 = t50 * qJD(2);
t74 = t48 * qJD(2);
t28 = -t45 * t74 + t46 * t73;
t55 = -qJ(5) * t28 - qJD(5) * t33;
t62 = qJD(2) * t81;
t26 = t50 * qJD(3) + t48 * t62;
t51 = -t48 * qJD(3) + t50 * t62;
t13 = t46 * t26 + t45 * t51;
t42 = pkin(2) * t74;
t14 = t27 * pkin(3) - t28 * pkin(7) + t42;
t64 = -t47 * t13 + t49 * t14;
t78 = qJ(5) * t33;
t1 = t27 * pkin(4) + t55 * t49 + (-t20 + (-t19 + t78) * t47) * qJD(4) + t64;
t75 = qJD(4) * t49;
t68 = t33 * t75;
t72 = t49 * t13 + t47 * t14 + t19 * t75;
t2 = -qJ(5) * t68 + (-qJD(4) * t22 + t55) * t47 + t72;
t63 = t49 * t19 - t47 * t22;
t5 = t32 * pkin(4) - t49 * t78 + t63;
t6 = -t47 * t78 + t80;
t88 = -t1 * t49 - t2 * t47 + (t47 * t5 - t49 * t6) * qJD(4);
t87 = 0.2e1 * qJD(4);
t86 = t33 * t47;
t85 = t33 * t49;
t84 = t47 * t27;
t83 = t49 * t27;
t82 = t49 * t28;
t43 = t47 ^ 2;
t44 = t49 ^ 2;
t79 = t43 - t44;
t40 = t45 * pkin(2) + pkin(7);
t77 = qJ(5) + t40;
t76 = qJD(4) * t47;
t71 = -0.2e1 * pkin(1) * qJD(2);
t41 = -t46 * pkin(2) - pkin(3);
t70 = t41 * t87;
t69 = pkin(4) * t76;
t67 = t47 * t75;
t65 = -0.4e1 * t47 * t85;
t12 = t45 * t26 - t46 * t51;
t21 = -t46 * t36 - t45 * t37;
t61 = t79 * qJD(4);
t60 = qJD(4) * t77;
t57 = -t27 * t40 + t28 * t41;
t56 = t32 * t40 - t33 * t41;
t54 = t32 * t75 + t84;
t53 = t47 * t28 + t68;
t52 = t33 * t76 - t82;
t35 = -t49 * pkin(4) + t41;
t31 = t33 ^ 2;
t30 = t77 * t49;
t29 = t77 * t47;
t24 = -t47 * qJD(5) - t49 * t60;
t23 = t49 * qJD(5) - t47 * t60;
t18 = -t32 * t76 + t83;
t11 = pkin(4) * t86 + t21;
t7 = t53 * pkin(4) + t12;
t4 = -t80 * qJD(4) + t64;
t3 = t22 * t76 - t72;
t8 = [0, 0, 0, 0.2e1 * t48 * t73, 0.2e1 * (-t48 ^ 2 + t50 ^ 2) * qJD(2), 0, 0, 0, t48 * t71, t50 * t71, 0.2e1 * t12 * t33 - 0.2e1 * t13 * t32 + 0.2e1 * t21 * t28 - 0.2e1 * t22 * t27, 0.2e1 * t21 * t12 + 0.2e1 * t22 * t13 + 0.2e1 * t66 * t42, 0.2e1 * t44 * t33 * t28 - 0.2e1 * t31 * t67, t79 * t31 * t87 + t28 * t65, -0.2e1 * t52 * t32 + 0.2e1 * t33 * t83, -0.2e1 * t53 * t32 - 0.2e1 * t33 * t84, 0.2e1 * t32 * t27, 0.2e1 * t12 * t86 + 0.2e1 * t53 * t21 + 0.2e1 * t63 * t27 + 0.2e1 * t4 * t32, 0.2e1 * t12 * t85 - 0.2e1 * t52 * t21 - 0.2e1 * t80 * t27 + 0.2e1 * t3 * t32, 0.2e1 * (-t47 * t6 - t49 * t5) * t28 + 0.2e1 * t88 * t33, 0.2e1 * t5 * t1 + 0.2e1 * t11 * t7 + 0.2e1 * t6 * t2; 0, 0, 0, 0, 0, t73, -t74, 0, -pkin(6) * t73, pkin(6) * t74, (-t27 * t45 - t28 * t46) * pkin(2), (-t12 * t46 + t13 * t45) * pkin(2), -t33 * t61 + t47 * t82, qJD(4) * t65 - t79 * t28, t54, t18, 0, -t12 * t49 + t57 * t47 + (t21 * t47 - t56 * t49) * qJD(4), t12 * t47 + t57 * t49 + (t21 * t49 + t56 * t47) * qJD(4), (-t24 * t33 + t28 * t29 + t2 + (-t30 * t33 - t5) * qJD(4)) * t49 + (-t23 * t33 - t28 * t30 - t1 + (-t29 * t33 - t6) * qJD(4)) * t47, -t1 * t29 + t11 * t69 + t2 * t30 + t6 * t23 + t5 * t24 + t7 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t67, -0.2e1 * t61, 0, 0, 0, t47 * t70, t49 * t70, 0.2e1 * t23 * t49 - 0.2e1 * t24 * t47 + 0.2e1 * (t29 * t49 - t30 * t47) * qJD(4), 0.2e1 * t30 * t23 - 0.2e1 * t29 * t24 + 0.2e1 * t35 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, t18, -t54, (-t43 - t44) * t28, -t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t47 + t24 * t49 + (t29 * t47 + t30 * t49) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t53, t27, t4, t3, t52 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t76, 0, -t40 * t75, t40 * t76, -pkin(4) * t75, t24 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t75, 0, -t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
