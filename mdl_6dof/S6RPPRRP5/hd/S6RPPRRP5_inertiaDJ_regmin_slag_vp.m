% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRRP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:08:46
% EndTime: 2019-03-09 02:08:48
% DurationCPUTime: 0.71s
% Computational Cost: add. (538->119), mult. (1121->211), div. (0->0), fcn. (794->4), ass. (0->73)
t35 = -pkin(7) + qJ(2);
t38 = sin(qJ(4));
t40 = cos(qJ(4));
t70 = t40 * qJD(4);
t81 = t38 * qJD(2) + t35 * t70;
t32 = t38 ^ 2;
t34 = t40 ^ 2;
t57 = (t32 - t34) * qJD(4);
t37 = sin(qJ(5));
t31 = t37 ^ 2;
t39 = cos(qJ(5));
t33 = t39 ^ 2;
t78 = t31 - t33;
t58 = t78 * qJD(5);
t29 = t38 * qJD(4);
t47 = -t40 * qJD(2) + t35 * t29;
t53 = pkin(4) * t40 + pkin(8) * t38;
t80 = qJD(5) * t53 + t47;
t36 = pkin(1) + qJ(3);
t79 = -qJ(6) - pkin(8);
t77 = t31 + t33;
t75 = t32 + t34;
t74 = qJ(6) * t40;
t73 = qJD(5) * t37;
t30 = qJD(5) * t39;
t72 = qJD(5) * t40;
t69 = qJ(2) * qJD(2);
t68 = -0.2e1 * pkin(4) * qJD(5);
t25 = t39 * t38 * t35;
t28 = pkin(5) * t73;
t67 = t37 * t72;
t66 = t38 * t30;
t65 = t39 * t72;
t64 = t35 * t73;
t63 = t37 * t30;
t62 = t39 * t29;
t61 = t38 * t70;
t59 = qJD(5) * t79;
t56 = 0.2e1 * t61;
t20 = t53 * qJD(4) + qJD(3);
t52 = t38 * pkin(4) - t40 * pkin(8);
t21 = t52 + t36;
t55 = t37 * t20 + t21 * t30 + t81 * t39;
t54 = t37 * t62;
t18 = t39 * t21;
t5 = -t39 * t74 + t18 + (-t35 * t37 + pkin(5)) * t38;
t17 = t37 * t21;
t6 = -t37 * t74 + t17 + t25;
t51 = t37 * t6 + t39 * t5;
t50 = t37 * t5 - t39 * t6;
t23 = t79 * t37;
t24 = t79 * t39;
t49 = t23 * t37 + t24 * t39;
t46 = t62 + t67;
t15 = t37 * t70 + t66;
t16 = t37 * t29 - t65;
t45 = -qJD(6) * t40 + (qJ(6) * qJD(4) - qJD(5) * t35) * t38;
t44 = t52 * qJD(4) - t35 * t72;
t11 = t39 * t20;
t1 = pkin(5) * t70 + t11 + t45 * t39 + ((-t21 + t74) * qJD(5) - t81) * t37;
t2 = -qJ(6) * t65 + t45 * t37 + t55;
t43 = t50 * qJD(5) - t1 * t39 - t2 * t37;
t8 = t39 * qJD(6) + t37 * t59;
t9 = -t37 * qJD(6) + t39 * t59;
t42 = -t9 * t37 + t8 * t39 + (-t23 * t39 + t24 * t37) * qJD(5);
t41 = 0.2e1 * qJD(2);
t27 = -t39 * pkin(5) - pkin(4);
t19 = (pkin(5) * t37 - t35) * t40;
t13 = t38 * t73 - t39 * t70;
t7 = -t16 * pkin(5) + t47;
t4 = -t35 * t66 + t11 + (-qJD(5) * t21 - t81) * t37;
t3 = t38 * t64 - t55;
t10 = [0, 0, 0, 0, t41, 0.2e1 * t69, t41, 0.2e1 * qJD(3), 0.2e1 * t36 * qJD(3) + 0.2e1 * t69, -0.2e1 * t61, 0.2e1 * t57, 0, 0, 0, 0.2e1 * qJD(3) * t38 + 0.2e1 * t36 * t70, 0.2e1 * qJD(3) * t40 - 0.2e1 * t36 * t29, -0.2e1 * t33 * t61 - 0.2e1 * t34 * t63, 0.2e1 * t34 * t58 + 0.4e1 * t40 * t54, -0.2e1 * t38 * t67 - 0.2e1 * t39 * t57, 0.2e1 * t37 * t57 - 0.2e1 * t38 * t65, t56, -0.2e1 * t34 * qJD(2) * t37 + 0.2e1 * t18 * t70 + 0.2e1 * t4 * t38 + 0.2e1 * (-t34 * t30 + t37 * t61) * t35, 0.2e1 * t3 * t38 + 0.2e1 * (-qJD(2) * t39 + t64) * t34 + 0.2e1 * (-t17 + t25) * t70, 0.2e1 * t51 * t29 + 0.2e1 * t40 * t43, 0.2e1 * t5 * t1 + 0.2e1 * t19 * t7 + 0.2e1 * t6 * t2; 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), 0, 0, 0, 0, 0, -t70, t29, 0, 0, 0, 0, 0, t13, t15, -t77 * t29, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75 * t30, t75 * t73, 0 (-qJD(4) * t50 - t7) * t40 + (qJD(4) * t19 - qJD(5) * t51 - t1 * t37 + t2 * t39) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-0.1e1 + t77) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t70, 0, -t47, -t81, -t40 * t58 - t54, t78 * t29 - 0.4e1 * t40 * t63, t15, -t13, 0, t44 * t37 - t80 * t39, t80 * t37 + t44 * t39 (t23 * t29 - t40 * t9 + t2 + (t24 * t40 - t5) * qJD(5)) * t39 + (-t24 * t29 - t40 * t8 - t1 + (t23 * t40 - t6) * qJD(5)) * t37, t1 * t23 + t19 * t28 - t2 * t24 + t7 * t27 + t5 * t9 + t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5) * t49 - t37 * t8 - t39 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t70, 0, 0, 0, 0, 0, -t46, t16, t77 * t70 (-qJD(4) * t49 - t28) * t40 + (qJD(4) * t27 + t42) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t63, -0.2e1 * t58, 0, 0, 0, t37 * t68, t39 * t68, 0.2e1 * t42, 0.2e1 * t23 * t9 - 0.2e1 * t24 * t8 + 0.2e1 * t27 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t16, t70, t4, t3, t46 * pkin(5), t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t30, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t13, 0, -t15 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t73, 0, -pkin(8) * t30, pkin(8) * t73, -pkin(5) * t30, t9 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
