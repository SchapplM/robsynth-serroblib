% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:45:24
% EndTime: 2019-03-09 02:45:26
% DurationCPUTime: 0.51s
% Computational Cost: add. (319->93), mult. (688->171), div. (0->0), fcn. (474->6), ass. (0->66)
t41 = sin(qJ(3));
t27 = t41 * qJ(4);
t43 = cos(qJ(3));
t76 = -t43 * pkin(3) - t27;
t20 = sin(pkin(9)) * pkin(1) + pkin(7);
t75 = -qJ(5) + t20;
t42 = cos(qJ(6));
t36 = t42 ^ 2;
t40 = sin(qJ(6));
t71 = t40 ^ 2 - t36;
t54 = t71 * qJD(6);
t73 = -pkin(4) - pkin(8);
t33 = -pkin(3) + t73;
t24 = t41 * qJD(3);
t56 = t20 * t24;
t62 = qJ(5) * qJD(3);
t6 = t43 * qJD(5) - t41 * t62 + t56;
t39 = qJ(4) + pkin(5);
t72 = t39 * t43;
t74 = (t33 * t41 + t72) * qJD(6) + t6;
t45 = 2 * qJD(4);
t21 = -cos(pkin(9)) * pkin(1) - pkin(2);
t37 = t43 ^ 2;
t70 = t41 ^ 2 - t37;
t69 = qJ(4) * t43;
t68 = qJD(3) * t42;
t26 = qJD(6) * t40;
t67 = qJD(6) * t42;
t66 = qJD(6) * t43;
t16 = t75 * t43;
t65 = t16 * qJD(6);
t64 = t41 * qJD(4);
t25 = t43 * qJD(3);
t63 = t43 * qJD(4);
t61 = 0.2e1 * qJD(3) * t21;
t60 = t40 * t66;
t59 = t42 * t66;
t58 = t40 * t67;
t57 = t41 * t25;
t55 = t42 * t24;
t14 = t21 + t76;
t53 = t70 * qJD(3);
t52 = -pkin(3) * t24 + t64;
t51 = t40 * t55;
t8 = t43 * pkin(4) - t14;
t15 = t75 * t41;
t4 = t41 * pkin(5) + t43 * pkin(8) + t8;
t50 = t42 * t15 + t40 * t4;
t49 = t40 * t15 - t42 * t4;
t47 = t76 * qJD(3) + t63;
t46 = -t63 - t65 + (-t33 * t43 + t39 * t41) * qJD(3);
t44 = -pkin(3) - pkin(4);
t32 = qJ(4) * t45;
t18 = 0.2e1 * t57;
t17 = t20 * t25;
t13 = -t40 * t24 + t59;
t12 = -t40 * t25 - t41 * t67;
t11 = t55 + t60;
t10 = -t42 * t25 + t41 * t26;
t9 = qJ(4) * t25 + t52;
t7 = -t41 * qJD(5) - t43 * t62 + t17;
t5 = (-pkin(4) * t41 + t69) * qJD(3) + t52;
t3 = (t73 * t41 + t72) * qJD(3) + t52;
t2 = -t50 * qJD(6) + t42 * t3 - t40 * t7;
t1 = t49 * qJD(6) - t40 * t3 - t42 * t7;
t19 = [0, 0, 0, 0, t18, -0.2e1 * t53, 0, 0, 0, t41 * t61, t43 * t61, 0.2e1 * t14 * t24 + 0.2e1 * t9 * t43, 0, -0.2e1 * t14 * t25 + 0.2e1 * t9 * t41, -0.2e1 * t14 * t9, 0.2e1 * t8 * t25 + 0.2e1 * t5 * t41, 0.2e1 * t8 * t24 - 0.2e1 * t5 * t43, -0.2e1 * t7 * t41 + 0.2e1 * t6 * t43 + 0.2e1 * (-t15 * t43 + t16 * t41) * qJD(3), 0.2e1 * t15 * t7 - 0.2e1 * t16 * t6 + 0.2e1 * t8 * t5, -0.2e1 * t36 * t57 - 0.2e1 * t37 * t58, 0.2e1 * t37 * t54 + 0.4e1 * t43 * t51, 0.2e1 * t41 * t60 + 0.2e1 * t70 * t68, -0.2e1 * t40 * t53 + 0.2e1 * t41 * t59, t18, 0.2e1 * (t16 * t40 * qJD(3) + t2) * t41 + 0.2e1 * (-t49 * qJD(3) + t6 * t40 - t42 * t65) * t43, 0.2e1 * (t16 * t68 + t1) * t41 + 0.2e1 * (-t50 * qJD(3) + t40 * t65 + t6 * t42) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6 * t41 - t7 * t43 + (t15 * t41 + t16 * t43) * qJD(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t25, -t24, 0, -t17, t56, -t17, t47, -t56, t47 * t20, -t6, t7, -t63 + (-t43 * t44 + t27) * qJD(3), -t6 * qJ(4) + t16 * qJD(4) + t7 * t44, -t43 * t54 - t51, t71 * t24 - 0.4e1 * t43 * t58, t12, t10, 0, t46 * t40 - t74 * t42, t74 * t40 + t46 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t25, -t24, 0, t25, t9, t25, t24, 0, t64 + (t41 * t44 + t69) * qJD(3), 0, 0, 0, 0, 0, -t10, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t32, t45, 0, 0, t32, 0.2e1 * t58, -0.2e1 * t54, 0, 0, 0, 0.2e1 * qJD(4) * t42 - 0.2e1 * t39 * t26, -0.2e1 * qJD(4) * t40 - 0.2e1 * t39 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, t17, 0, 0, -t25, t7, 0, 0, 0, 0, 0, t12, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t24, 0, t5, 0, 0, 0, 0, 0, -t10, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t13, t25, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, t26, 0, -t33 * t67, t33 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t19;
