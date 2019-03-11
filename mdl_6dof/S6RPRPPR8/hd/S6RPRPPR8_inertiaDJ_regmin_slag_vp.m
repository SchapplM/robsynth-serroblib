% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPPR8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:59:59
% EndTime: 2019-03-09 03:00:01
% DurationCPUTime: 0.53s
% Computational Cost: add. (323->95), mult. (647->166), div. (0->0), fcn. (424->4), ass. (0->70)
t36 = sin(qJ(3));
t39 = -pkin(3) - pkin(4);
t74 = t36 * t39;
t29 = -pkin(8) + t39;
t34 = qJ(4) + pkin(5);
t38 = cos(qJ(3));
t43 = t29 * t38 - t34 * t36;
t37 = cos(qJ(6));
t32 = t37 ^ 2;
t35 = sin(qJ(6));
t67 = t35 ^ 2 - t32;
t48 = qJD(6) * t67;
t25 = t38 * qJD(3);
t40 = -pkin(1) - pkin(7);
t62 = qJ(5) + t40;
t8 = t36 * qJD(5) + t62 * t25;
t73 = t43 * qJD(6) - t8;
t72 = 2 * qJD(2);
t41 = 2 * qJD(4);
t71 = t36 * pkin(3);
t70 = t29 * t36;
t31 = t36 ^ 2;
t33 = t38 ^ 2;
t66 = t31 - t33;
t65 = t31 + t33;
t64 = qJ(4) * t36;
t27 = t38 * qJ(4);
t63 = -qJ(2) + t27;
t61 = qJD(3) * t37;
t26 = qJD(6) * t35;
t60 = qJD(6) * t37;
t59 = qJD(6) * t38;
t17 = t62 * t36;
t58 = t17 * qJD(6);
t23 = t36 * qJD(3);
t57 = t36 * qJD(4);
t56 = t38 * qJD(4) - qJD(2);
t55 = qJ(2) * qJD(3);
t54 = t35 * t59;
t53 = t37 * t59;
t52 = t35 * t60;
t51 = t37 * t25;
t50 = t36 * t25;
t49 = t40 * t25;
t47 = t66 * qJD(3);
t46 = t35 * t51;
t18 = t62 * t38;
t6 = t38 * pkin(5) + t63 + t70;
t45 = -t37 * t18 + t35 * t6;
t44 = -t35 * t18 - t37 * t6;
t11 = t57 + (t27 - t71) * qJD(3);
t42 = t57 - t58 + (t34 * t38 + t70) * qJD(3);
t28 = qJ(4) * t41;
t22 = t40 * t23;
t20 = -0.2e1 * t50;
t19 = -t63 + t71;
t16 = t63 + t74;
t15 = -t35 * t23 + t53;
t14 = -t35 * t25 - t36 * t60;
t13 = t37 * t23 + t54;
t12 = -t36 * t26 + t51;
t10 = (pkin(3) * t38 + t64) * qJD(3) - t56;
t9 = t57 + (t27 + t74) * qJD(3);
t7 = qJ(5) * t23 - t38 * qJD(5) + t22;
t5 = (t39 * t38 - t64) * qJD(3) + t56;
t4 = t43 * qJD(3) + t56;
t3 = t8 * t36 - t7 * t38 + (t17 * t38 - t18 * t36) * qJD(3);
t2 = -t45 * qJD(6) - t35 * t7 + t37 * t4;
t1 = t44 * qJD(6) - t35 * t4 - t37 * t7;
t21 = [0, 0, 0, 0, t72, qJ(2) * t72, t20, 0.2e1 * t47, 0, 0, 0, 0.2e1 * qJD(2) * t36 + 0.2e1 * t38 * t55, 0.2e1 * qJD(2) * t38 - 0.2e1 * t36 * t55, 0.2e1 * t10 * t36 + 0.2e1 * t19 * t25, 0, -0.2e1 * t10 * t38 + 0.2e1 * t19 * t23, 0.2e1 * t19 * t10, -0.2e1 * t16 * t23 + 0.2e1 * t5 * t38, 0.2e1 * t16 * t25 + 0.2e1 * t5 * t36, 0.2e1 * t3, 0.2e1 * t16 * t5 + 0.2e1 * t17 * t8 - 0.2e1 * t18 * t7, -0.2e1 * t31 * t52 + 0.2e1 * t32 * t50, 0.2e1 * t31 * t48 - 0.4e1 * t36 * t46, -0.2e1 * t36 * t54 - 0.2e1 * t66 * t61, 0.2e1 * t35 * t47 - 0.2e1 * t36 * t53, t20, 0.2e1 * (t17 * t35 * qJD(3) + t2) * t38 + 0.2e1 * (t44 * qJD(3) + t8 * t35 + t37 * t58) * t36, 0.2e1 * (t17 * t61 + t1) * t38 + 0.2e1 * (t45 * qJD(3) - t35 * t58 + t8 * t37) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, t65 * t60, -t65 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t25, 0, -t22, -t49, -t22, -t11, t49, t11 * t40, t8, t7, t9, t8 * qJ(4) + t17 * qJD(4) + t7 * t39, t36 * t48 - t46, t67 * t25 + 0.4e1 * t36 * t52, -t15, t13, 0, t42 * t35 - t37 * t73, t35 * t73 + t42 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t25, -t23, 0, t25, t11, t25, t23, 0, t9, 0, 0, 0, 0, 0, t12, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t28, t41, 0, 0, t28, 0.2e1 * t52, -0.2e1 * t48, 0, 0, 0, 0.2e1 * qJD(4) * t37 - 0.2e1 * t34 * t26, -0.2e1 * qJD(4) * t35 - 0.2e1 * t34 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, t22, 0, 0, t23, t7, 0, 0, 0, 0, 0, -t15, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, 0, t23, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t25, 0, t5, 0, 0, 0, 0, 0, -t13, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t14, -t23, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, t26, 0, -t29 * t60, t29 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t21;
