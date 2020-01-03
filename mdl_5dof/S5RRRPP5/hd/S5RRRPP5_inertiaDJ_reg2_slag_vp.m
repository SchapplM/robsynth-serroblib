% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPP5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:58:28
% EndTime: 2019-12-31 20:58:31
% DurationCPUTime: 0.67s
% Computational Cost: add. (806->108), mult. (1932->170), div. (0->0), fcn. (1579->4), ass. (0->62)
t52 = cos(qJ(2));
t76 = cos(qJ(3));
t63 = t76 * qJD(3);
t65 = t76 * t52;
t50 = sin(qJ(3));
t51 = sin(qJ(2));
t75 = t50 * t51;
t78 = qJD(2) + qJD(3);
t18 = -qJD(2) * t65 - t52 * t63 + t78 * t75;
t74 = t50 * t52;
t29 = t76 * t51 + t74;
t19 = t78 * t29;
t28 = -t65 + t75;
t4 = 0.2e1 * t18 * t28 - 0.2e1 * t29 * t19;
t54 = 2 * qJD(4);
t53 = -pkin(3) - pkin(4);
t77 = -pkin(7) - pkin(6);
t35 = t77 * t52;
t21 = -t76 * t35 + t77 * t75;
t48 = pkin(2) * t63;
t37 = t48 + qJD(4);
t43 = t50 * pkin(2) + qJ(4);
t73 = t37 * qJ(4) + t43 * qJD(4);
t72 = qJD(3) * t50;
t71 = t51 * qJD(2);
t70 = t52 * qJD(2);
t12 = 0.2e1 * t28 * t19;
t69 = -0.2e1 * pkin(1) * qJD(2);
t68 = pkin(2) * t71;
t47 = pkin(2) * t72;
t67 = t51 * t70;
t46 = -t52 * pkin(2) - pkin(1);
t62 = t77 * t76;
t31 = t51 * t62;
t20 = -t50 * t35 - t31;
t59 = qJD(2) * t62;
t64 = t77 * qJD(2);
t8 = -qJD(3) * t31 - t35 * t72 - t51 * t59 - t64 * t74;
t9 = -t35 * t63 - t52 * t59 + (qJD(3) * t77 + t64) * t75;
t66 = t20 * t9 - t21 * t8;
t45 = -t76 * pkin(2) - pkin(3);
t60 = -qJ(4) * t19 - qJD(4) * t28;
t58 = t29 * qJ(4) - t46;
t57 = -0.2e1 * t20 * t18 - 0.2e1 * t21 * t19 + 0.2e1 * t8 * t28 + 0.2e1 * t9 * t29;
t56 = -t18 * qJ(4) + t29 * qJD(4) - t68;
t55 = t43 * t19 + t37 * t28 - t29 * t47;
t49 = qJ(4) * t54;
t42 = -0.2e1 * t47;
t41 = -pkin(4) + t45;
t38 = t54 + t48;
t33 = 0.2e1 * t37;
t26 = t43 * t37;
t16 = t28 * pkin(3) - t58;
t15 = t28 * qJ(5) + t21;
t14 = -t29 * qJ(5) + t20;
t13 = -0.2e1 * t29 * t18;
t11 = t53 * t28 + t58;
t6 = t19 * pkin(3) - t56;
t3 = t18 * qJ(5) - t29 * qJD(5) + t9;
t2 = t19 * qJ(5) + t28 * qJD(5) - t8;
t1 = t53 * t19 + t56;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t67, 0.2e1 * (-t51 ^ 2 + t52 ^ 2) * qJD(2), 0, -0.2e1 * t67, 0, 0, t51 * t69, t52 * t69, 0, 0, t13, t4, 0, t12, 0, 0, 0.2e1 * t46 * t19 + 0.2e1 * t28 * t68, -0.2e1 * t46 * t18 + 0.2e1 * t29 * t68, t57, 0.2e1 * t46 * t68 + 0.2e1 * t66, t13, 0, -t4, 0, 0, t12, 0.2e1 * t16 * t19 + 0.2e1 * t6 * t28, t57, 0.2e1 * t16 * t18 - 0.2e1 * t6 * t29, 0.2e1 * t16 * t6 + 0.2e1 * t66, t13, -t4, 0, t12, 0, 0, -0.2e1 * t1 * t28 - 0.2e1 * t11 * t19, 0.2e1 * t1 * t29 - 0.2e1 * t11 * t18, 0.2e1 * t14 * t18 + 0.2e1 * t15 * t19 + 0.2e1 * t2 * t28 - 0.2e1 * t3 * t29, 0.2e1 * t11 * t1 + 0.2e1 * t14 * t3 + 0.2e1 * t15 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, -t71, 0, -pkin(6) * t70, pkin(6) * t71, 0, 0, 0, 0, -t18, 0, -t19, 0, -t9, t8, (t76 * t18 - t19 * t50 + (-t76 * t28 + t29 * t50) * qJD(3)) * pkin(2), (-t76 * t9 - t50 * t8 + (t20 * t50 + t76 * t21) * qJD(3)) * pkin(2), 0, -t18, 0, 0, t19, 0, -t9, -t45 * t18 - t55, -t8, t20 * t47 + t21 * t37 - t8 * t43 + t9 * t45, 0, 0, t18, 0, -t19, 0, -t3, t2, t41 * t18 + t55, t14 * t47 + t15 * t37 + t2 * t43 + t3 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -0.2e1 * t48, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t33, 0.2e1 * t45 * t47 + 0.2e1 * t26, 0, 0, 0, 0, 0, 0, t42, t33, 0, 0.2e1 * t41 * t47 + 0.2e1 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, -t19, 0, -t9, t8, 0, 0, 0, -t18, 0, 0, t19, 0, -t9, pkin(3) * t18 + t60, -t8, -t9 * pkin(3) - t8 * qJ(4) + t21 * qJD(4), 0, 0, t18, 0, -t19, 0, -t3, t2, t53 * t18 - t60, t2 * qJ(4) + t15 * qJD(4) + t3 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t48, 0, 0, 0, 0, 0, 0, 0, 0, -t47, 0, t38, -pkin(3) * t47 + t73, 0, 0, 0, 0, 0, 0, -t47, t38, 0, t53 * t47 + t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t49, 0, 0, 0, 0, 0, 0, 0, t54, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, t18, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t18, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
