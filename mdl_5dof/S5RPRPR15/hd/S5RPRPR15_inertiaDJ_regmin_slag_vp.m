% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR15_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:37:24
% EndTime: 2019-12-31 18:37:26
% DurationCPUTime: 0.64s
% Computational Cost: add. (475->103), mult. (1164->208), div. (0->0), fcn. (981->6), ass. (0->73)
t43 = sin(pkin(8));
t44 = cos(pkin(8));
t45 = sin(qJ(5));
t47 = cos(qJ(5));
t84 = -t45 * t43 + t47 * t44;
t23 = t84 * qJD(5);
t48 = cos(qJ(3));
t87 = t48 * t23;
t46 = sin(qJ(3));
t86 = t23 * t46;
t29 = t43 * t47 + t44 * t45;
t85 = t29 * t48;
t40 = -pkin(4) * t44 - pkin(3);
t83 = 0.2e1 * t40;
t82 = 2 * qJD(2);
t81 = t44 * pkin(7);
t80 = t46 * pkin(3);
t24 = t29 * qJD(5);
t79 = t24 * t46;
t49 = -pkin(1) - pkin(6);
t77 = t46 * t49;
t75 = pkin(7) + qJ(4);
t21 = -t48 * qJD(4) + qJD(2) + (pkin(3) * t48 + qJ(4) * t46) * qJD(3);
t70 = t48 * qJD(3);
t64 = t49 * t70;
t13 = t21 * t43 + t44 * t64;
t58 = -qJ(4) * t48 + t80;
t30 = qJ(2) + t58;
t37 = t44 * t77;
t16 = t30 * t43 + t37;
t74 = t43 ^ 2 + t44 ^ 2;
t73 = qJD(4) * t46;
t71 = t46 * qJD(3);
t69 = qJ(2) * qJD(3);
t68 = t43 * t77;
t67 = t43 * t71;
t66 = t44 * t71;
t38 = t49 * t71;
t65 = t46 * t70;
t63 = -t43 * t49 + pkin(4);
t62 = t74 * t48;
t61 = t74 * qJD(4);
t60 = 0.2e1 * t65;
t59 = 0.2e1 * t61;
t26 = t44 * t30;
t11 = t46 * t63 - t48 * t81 + t26;
t14 = -pkin(7) * t43 * t48 + t16;
t57 = t11 * t47 - t14 * t45;
t56 = t11 * t45 + t14 * t47;
t18 = t44 * t21;
t12 = -t43 * t64 + t18;
t55 = -t12 * t43 + t13 * t44;
t15 = t26 - t68;
t54 = -t15 * t43 + t16 * t44;
t34 = t75 * t43;
t35 = t75 * t44;
t53 = -t34 * t47 - t35 * t45;
t52 = -t34 * t45 + t35 * t47;
t20 = t84 * t48;
t50 = qJD(3) * t84;
t27 = (pkin(4) * t43 - t49) * t48;
t22 = -pkin(4) * t67 + t38;
t10 = pkin(7) * t67 + t13;
t9 = -t45 * t66 - t47 * t67 + t87;
t8 = -qJD(3) * t85 - t86;
t7 = -qJD(5) * t85 - t46 * t50;
t6 = -t48 * t50 + t79;
t5 = -qJD(4) * t29 - qJD(5) * t52;
t4 = -qJD(4) * t84 - qJD(5) * t53;
t3 = t18 + (t46 * t81 + t48 * t63) * qJD(3);
t2 = -qJD(5) * t56 - t45 * t10 + t47 * t3;
t1 = -qJD(5) * t57 - t47 * t10 - t45 * t3;
t17 = [0, 0, 0, 0, t82, qJ(2) * t82, -0.2e1 * t65, 0.2e1 * (t46 ^ 2 - t48 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t46 + 0.2e1 * t48 * t69, 0.2e1 * qJD(2) * t48 - 0.2e1 * t46 * t69, 0.2e1 * t12 * t46 + 0.2e1 * (t15 + 0.2e1 * t68) * t70, -0.2e1 * t13 * t46 + 0.2e1 * (-t16 + 0.2e1 * t37) * t70, 0.2e1 * (-t12 * t44 - t13 * t43) * t48 + 0.2e1 * (t15 * t44 + t16 * t43) * t71, -0.2e1 * t49 ^ 2 * t65 + 0.2e1 * t12 * t15 + 0.2e1 * t13 * t16, 0.2e1 * t20 * t7, -0.2e1 * t20 * t9 - 0.2e1 * t7 * t85, 0.2e1 * t20 * t70 + 0.2e1 * t46 * t7, -0.2e1 * t46 * t9 - 0.2e1 * t70 * t85, t60, 0.2e1 * t2 * t46 + 0.2e1 * t22 * t85 + 0.2e1 * t27 * t9 + 0.2e1 * t57 * t70, 0.2e1 * t1 * t46 + 0.2e1 * t20 * t22 + 0.2e1 * t27 * t7 - 0.2e1 * t56 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 * t46 + (t54 - 0.2e1 * t77) * t70, 0, 0, 0, 0, 0, t8 * t46 - t48 * t9, t6 * t46 - t48 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-0.1e1 + t74) * t60, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t70, 0, -t38, -t64, -t43 * t73 + (t43 * t58 - t37) * qJD(3), -t44 * t73 + (t44 * t58 + t68) * qJD(3), t55, -pkin(3) * t38 + qJ(4) * t55 + qJD(4) * t54, t20 * t23 + t29 * t7, -t20 * t24 - t23 * t85 - t29 * t9 + t7 * t84, t29 * t70 + t86, t70 * t84 - t79, 0, -t22 * t84 + t24 * t27 + t40 * t9 + t46 * t5 + t53 * t70, t22 * t29 + t23 * t27 + t4 * t46 + t40 * t7 - t52 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t70, -t66, t67, qJD(3) * t62, t46 * t61 + (qJ(4) * t62 - t80) * qJD(3), 0, 0, 0, 0, 0, -t24 * t48 - t71 * t84, t29 * t71 - t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, qJ(4) * t59, 0.2e1 * t29 * t23, 0.2e1 * t23 * t84 - 0.2e1 * t24 * t29, 0, 0, 0, t24 * t83, t23 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t66, 0, t38, 0, 0, 0, 0, 0, t9, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t9, t70, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t24, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t17;
