% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPPRR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:58
% EndTime: 2019-03-09 01:34:00
% DurationCPUTime: 0.51s
% Computational Cost: add. (568->85), mult. (1207->178), div. (0->0), fcn. (1167->8), ass. (0->72)
t44 = cos(pkin(9));
t67 = t44 * qJD(2);
t30 = -qJD(4) + t67;
t41 = sin(pkin(10));
t43 = cos(pkin(10));
t63 = (t41 ^ 2 + t43 ^ 2) * t30;
t47 = cos(qJ(6));
t40 = t47 ^ 2;
t45 = sin(qJ(6));
t70 = t45 ^ 2 - t40;
t59 = qJD(6) * t70;
t42 = sin(pkin(9));
t33 = t42 * qJD(2);
t85 = 0.2e1 * t33;
t84 = 0.2e1 * qJD(2);
t46 = sin(qJ(5));
t80 = cos(qJ(5));
t26 = t80 * t41 + t46 * t43;
t48 = -pkin(1) - pkin(2);
t72 = t44 * qJ(2) + t42 * t48;
t28 = -qJ(4) + t72;
t81 = pkin(7) - t28;
t16 = t81 * t41;
t17 = t81 * t43;
t8 = t46 * t16 - t80 * t17;
t4 = t8 * qJD(5) + t26 * t30;
t83 = t4 * t45;
t82 = t4 * t47;
t64 = t80 * t43;
t74 = t46 * t41;
t49 = t64 - t74;
t19 = t49 * t42;
t15 = qJD(5) * t19;
t79 = t15 * t45;
t60 = qJD(5) * t80;
t21 = qJD(5) * t74 - t43 * t60;
t78 = t26 * t21;
t77 = t26 * t47;
t76 = t45 * t21;
t22 = t26 * qJD(5);
t75 = t45 * t22;
t73 = t47 * t21;
t69 = qJD(6) * t45;
t68 = qJD(6) * t47;
t66 = -0.2e1 * pkin(5) * qJD(6);
t65 = t45 * t68;
t62 = 0.4e1 * t45 * t77;
t61 = -t42 * qJ(2) + t44 * t48;
t58 = pkin(3) - t61;
t57 = -pkin(5) * t21 + pkin(8) * t22;
t56 = pkin(5) * t26 - pkin(8) * t49;
t20 = t43 * pkin(4) + t58;
t9 = pkin(5) * t49 + t26 * pkin(8) + t20;
t55 = t45 * t9 + t47 * t8;
t54 = t45 * t8 - t47 * t9;
t53 = t47 * t19 - t45 * t44;
t52 = t45 * t19 + t47 * t44;
t51 = t26 * t68 - t76;
t50 = t26 * t69 + t73;
t12 = -t49 * t68 + t75;
t11 = -t22 * t47 - t49 * t69;
t23 = t26 ^ 2;
t18 = t26 * t42;
t14 = t42 * t22;
t10 = -t22 * pkin(5) - t21 * pkin(8) + t33;
t7 = -t80 * t16 - t46 * t17;
t6 = -t53 * qJD(6) + t45 * t14;
t5 = t52 * qJD(6) + t47 * t14;
t3 = -t16 * t60 - t30 * t64 + (-qJD(5) * t17 + t30 * t41) * t46;
t2 = -t55 * qJD(6) + t47 * t10 + t45 * t3;
t1 = t54 * qJD(6) - t45 * t10 + t47 * t3;
t13 = [0, 0, 0, 0, t84, qJ(2) * t84, t85, 0.2e1 * t67 (-t61 * t42 + t72 * t44) * t84, t43 * t85, -0.2e1 * t41 * t33, -0.2e1 * t63, 0.2e1 * t28 * t63 + 0.2e1 * t58 * t33, -0.2e1 * t78, -0.2e1 * t21 * t49 - 0.2e1 * t26 * t22, 0, 0, 0, -0.2e1 * t20 * t22 + 0.2e1 * t33 * t49, 0.2e1 * t20 * t21 - 0.2e1 * t26 * t33, -0.2e1 * t23 * t65 - 0.2e1 * t40 * t78, t21 * t62 + 0.2e1 * t23 * t59, 0.2e1 * t22 * t77 + 0.2e1 * t49 * t50, -0.2e1 * t26 * t75 + 0.2e1 * t49 * t51, -0.2e1 * t49 * t22, 0.2e1 * t2 * t49 + 0.2e1 * t54 * t22 + 0.2e1 * t7 * t76 + 0.2e1 * (-t7 * t68 - t83) * t26, 0.2e1 * t1 * t49 + 0.2e1 * t55 * t22 + 0.2e1 * t7 * t73 + 0.2e1 * (t7 * t69 - t82) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t63 - t67) * t42, 0, 0, 0, 0, 0, t44 * t22, -t44 * t21, 0, 0, 0, 0, 0, -t51 * t18 + t52 * t22 - t26 * t79 + t49 * t6, -t15 * t77 + t50 * t18 + t53 * t22 + t49 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, -t22, t21, 0, 0, 0, 0, 0, t11, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t22, 0, -t4, t3, t26 * t59 + t45 * t73, qJD(6) * t62 - t70 * t21, -t12, t11, 0, -t82 + t57 * t45 + (t45 * t7 + t56 * t47) * qJD(6), t83 + t57 * t47 + (-t56 * t45 + t47 * t7) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t14, 0, 0, 0, 0, 0, -t15 * t47 + t18 * t69, t18 * t68 + t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t21, 0, 0, 0, 0, 0, t11, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t65, -0.2e1 * t59, 0, 0, 0, t45 * t66, t47 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t51, -t22, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, -t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, -t69, 0, -pkin(8) * t68, pkin(8) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t13;
