% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPPR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:42:39
% EndTime: 2019-03-09 02:42:40
% DurationCPUTime: 0.61s
% Computational Cost: add. (763->103), mult. (1676->188), div. (0->0), fcn. (1512->8), ass. (0->72)
t42 = sin(pkin(10));
t43 = cos(pkin(10));
t48 = cos(qJ(3));
t71 = t48 * qJD(3);
t46 = sin(qJ(3));
t72 = t46 * qJD(3);
t26 = -t42 * t72 + t43 * t71;
t31 = t42 * t48 + t43 * t46;
t45 = sin(qJ(6));
t47 = cos(qJ(6));
t74 = qJD(6) * t47;
t16 = t26 * t45 + t31 * t74;
t40 = t45 ^ 2;
t77 = -t47 ^ 2 + t40;
t65 = t77 * qJD(6);
t85 = 2 * qJD(5);
t84 = pkin(4) + pkin(8);
t30 = t42 * t46 - t43 * t48;
t25 = t31 * qJD(3);
t36 = sin(pkin(9)) * pkin(1) + pkin(7);
t76 = qJ(4) + t36;
t64 = qJD(3) * t76;
t23 = t48 * qJD(4) - t46 * t64;
t52 = -t46 * qJD(4) - t48 * t64;
t9 = t43 * t23 + t42 * t52;
t5 = -t25 * pkin(5) + t9;
t83 = t5 * t30;
t82 = t30 * t25;
t81 = t31 * t26;
t80 = t45 * t25;
t78 = t47 * t25;
t75 = qJD(6) * t45;
t73 = t31 * qJD(5);
t70 = t45 * t78;
t69 = 0.2e1 * t71;
t39 = pkin(3) * t72;
t66 = t45 * t74;
t38 = -cos(pkin(9)) * pkin(1) - pkin(2);
t37 = -pkin(3) * t43 - pkin(4);
t8 = t23 * t42 - t43 * t52;
t27 = t76 * t46;
t28 = t76 * t48;
t17 = t27 * t43 + t28 * t42;
t18 = -t27 * t42 + t28 * t43;
t63 = t17 * t8 + t18 * t9;
t10 = pkin(5) * t31 + t17;
t59 = -t48 * pkin(3) + t38;
t54 = -t31 * qJ(5) + t59;
t7 = t30 * t84 + t54;
t62 = t10 * t47 - t45 * t7;
t61 = t10 * t45 + t47 * t7;
t34 = pkin(3) * t42 + qJ(5);
t60 = -qJD(5) * t30 - t25 * t34;
t15 = t30 * t74 + t80;
t58 = t30 * t75 - t78;
t57 = -t26 * t47 + t31 * t75;
t56 = -t26 * qJ(5) + t39 - t73;
t55 = 0.2e1 * t81 + 0.2e1 * t82;
t33 = -pkin(8) + t37;
t53 = t5 + (t30 * t34 - t31 * t33) * qJD(6);
t51 = t17 * t25 + t18 * t26 + t30 * t8 + t31 * t9;
t11 = -pkin(5) * t30 + t18;
t50 = -qJD(6) * t11 - t26 * t33 - t60;
t49 = 0.2e1 * t17 * t26 - 0.2e1 * t18 * t25 - 0.2e1 * t30 * t9 + 0.2e1 * t31 * t8;
t29 = t30 ^ 2;
t12 = pkin(4) * t30 + t54;
t6 = pkin(4) * t25 + t56;
t4 = pkin(5) * t26 + t8;
t3 = t25 * t84 + t56;
t2 = -qJD(6) * t61 - t45 * t3 + t47 * t4;
t1 = -qJD(6) * t62 - t47 * t3 - t45 * t4;
t13 = [0, 0, 0, 0, t46 * t69, 0.2e1 * (-t46 ^ 2 + t48 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * t38 * t72, t38 * t69, t49, 0.2e1 * t39 * t59 + 0.2e1 * t63, t49, -0.2e1 * t12 * t25 - 0.2e1 * t30 * t6, -0.2e1 * t12 * t26 - 0.2e1 * t31 * t6, 0.2e1 * t12 * t6 + 0.2e1 * t63, 0.2e1 * t29 * t66 + 0.2e1 * t40 * t82, -0.2e1 * t29 * t65 + 0.4e1 * t30 * t70, 0.2e1 * t16 * t30 + 0.2e1 * t31 * t80, -0.2e1 * t30 * t57 + 0.2e1 * t31 * t78, 0.2e1 * t81, 0.2e1 * t11 * t58 + 0.2e1 * t2 * t31 + 0.2e1 * t26 * t62 - 0.2e1 * t47 * t83, 0.2e1 * t1 * t31 + 0.2e1 * t11 * t15 - 0.2e1 * t26 * t61 + 0.2e1 * t45 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, t55, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t71, -t72, 0, -t36 * t71, t36 * t72 (-t25 * t42 - t26 * t43) * pkin(3) (t42 * t9 - t43 * t8) * pkin(3), t26 * t37 + t60, t8, t9, qJD(5) * t18 + t34 * t9 + t37 * t8, -t30 * t65 + t70, -t25 * t77 - 0.4e1 * t30 * t66, -t57, -t16, 0, t45 * t53 - t47 * t50, t45 * t50 + t47 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t71, 0 (-t25 * t43 + t26 * t42) * pkin(3), 0, t25, t26, t25 * t37 + t26 * t34 + t73, 0, 0, 0, 0, 0, t16, -t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, t34 * t85, -0.2e1 * t66, 0.2e1 * t65, 0, 0, 0, 0.2e1 * qJD(5) * t45 + 0.2e1 * t34 * t74, 0.2e1 * qJD(5) * t47 - 0.2e1 * t34 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, -t25, -t26, t6, 0, 0, 0, 0, 0, -t16, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, t8, 0, 0, 0, 0, 0, -t57, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t58, t26, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t74, 0, -t33 * t75, -t33 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t13;
