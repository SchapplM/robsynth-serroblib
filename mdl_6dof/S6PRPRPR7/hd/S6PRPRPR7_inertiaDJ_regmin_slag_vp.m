% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRPR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:53:47
% EndTime: 2019-03-08 19:53:49
% DurationCPUTime: 0.72s
% Computational Cost: add. (311->112), mult. (846->206), div. (0->0), fcn. (703->8), ass. (0->83)
t43 = sin(qJ(4));
t46 = cos(qJ(4));
t85 = t46 * qJ(5);
t95 = -t43 * pkin(4) + t85;
t37 = t43 ^ 2;
t39 = t46 ^ 2;
t64 = (t37 - t39) * qJD(4);
t42 = sin(qJ(6));
t36 = t42 ^ 2;
t45 = cos(qJ(6));
t87 = -t45 ^ 2 + t36;
t63 = t87 * qJD(6);
t48 = -pkin(4) - pkin(9);
t49 = -pkin(2) - pkin(8);
t90 = pkin(5) - t49;
t23 = t90 * t43;
t75 = t23 * qJD(6);
t81 = qJD(5) * t43;
t94 = qJD(4) * (t43 * t48 + t85) + t75 + t81;
t41 = cos(pkin(6));
t40 = sin(pkin(6));
t47 = cos(qJ(2));
t88 = t40 * t47;
t13 = t41 * t43 + t46 * t88;
t14 = t41 * t46 - t43 * t88;
t44 = sin(qJ(2));
t89 = t40 * t44;
t28 = qJD(2) * t89;
t6 = t13 * qJD(4) - t43 * t28;
t7 = qJD(4) * t14 - t46 * t28;
t93 = -t6 * t43 - t7 * t46 + (t13 * t43 + t14 * t46) * qJD(4);
t92 = 2 * qJD(3);
t91 = 0.2e1 * qJD(5);
t84 = qJD(2) * t47;
t83 = qJD(4) * t42;
t82 = qJD(4) * t45;
t80 = qJD(6) * t42;
t79 = qJD(6) * t45;
t78 = qJD(6) * t46;
t77 = qJD(6) * t48;
t76 = t14 * qJD(6);
t34 = t43 * qJD(4);
t74 = t46 * qJD(4);
t73 = qJ(3) * qJD(4);
t72 = qJ(5) * qJD(6);
t71 = t42 * t78;
t70 = t45 * t78;
t69 = t40 * t84;
t68 = t45 * t74;
t67 = t42 * t79;
t66 = t43 * t74;
t65 = pkin(4) * t74 + qJ(5) * t34 + qJD(3);
t62 = qJD(6) * (t37 + t39);
t61 = t42 * t68;
t22 = qJ(3) - t95;
t19 = t43 * pkin(9) + t22;
t24 = t90 * t46;
t58 = t45 * t19 + t42 * t24;
t57 = t42 * t19 - t45 * t24;
t54 = t13 * t42 + t45 * t89;
t53 = t13 * t45 - t42 * t89;
t52 = -t6 * t42 + t45 * t76;
t51 = -t42 * t76 - t6 * t45;
t31 = t49 * t74;
t21 = -pkin(5) * t74 + t31;
t50 = t21 + (qJ(5) * t43 - t46 * t48) * qJD(6);
t12 = t95 * qJD(4) + t81;
t30 = t49 * t34;
t27 = -0.2e1 * t66;
t20 = -pkin(5) * t34 + t30;
t18 = -t42 * t34 + t70;
t17 = t42 * t74 + t43 * t79;
t16 = t45 * t34 + t71;
t15 = -t43 * t80 + t68;
t11 = -t46 * qJD(5) + t65;
t10 = (-t44 * t34 + t46 * t84) * t40;
t9 = (t43 * t84 + t44 * t74) * t40;
t8 = (qJD(4) * pkin(9) - qJD(5)) * t46 + t65;
t5 = t53 * qJD(6) + t7 * t42 + t45 * t69;
t4 = -t54 * qJD(6) - t42 * t69 + t7 * t45;
t3 = -t58 * qJD(6) + t45 * t20 - t42 * t8;
t2 = t57 * qJD(6) - t42 * t20 - t45 * t8;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t40 ^ 2 * t44 * t84 + 0.2e1 * t13 * t7 - 0.2e1 * t14 * t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t28, -t69, t28, t69 (qJD(3) * t44 + (-pkin(2) * t44 + qJ(3) * t47) * qJD(2)) * t40, 0, 0, 0, 0, 0, t9, t10, -t93, -t9, -t10 (t11 * t44 + t22 * t84) * t40 + t93 * t49, 0, 0, 0, 0, 0 (-t14 * t82 + t4) * t46 + (-qJD(4) * t53 - t51) * t43 (t14 * t83 - t5) * t46 + (qJD(4) * t54 + t52) * t43; 0, 0, 0, 0, 0, t92, qJ(3) * t92, t27, 0.2e1 * t64, 0, 0, 0, 0.2e1 * qJD(3) * t43 + 0.2e1 * t46 * t73, 0.2e1 * qJD(3) * t46 - 0.2e1 * t43 * t73, 0, -0.2e1 * t11 * t43 - 0.2e1 * t22 * t74, -0.2e1 * t11 * t46 + 0.2e1 * t22 * t34, 0.2e1 * t22 * t11, 0.2e1 * t36 * t66 + 0.2e1 * t37 * t67, -0.2e1 * t37 * t63 + 0.4e1 * t43 * t61, -0.2e1 * t42 * t64 + 0.2e1 * t43 * t70, -0.2e1 * t43 * t71 - 0.2e1 * t45 * t64, t27, 0.2e1 * (t23 * t82 + t3) * t46 + 0.2e1 * (qJD(4) * t57 - t21 * t45 - t42 * t75) * t43, 0.2e1 * (-t23 * t83 + t2) * t46 + 0.2e1 * (qJD(4) * t58 + t21 * t42 - t45 * t75) * t43; 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t62, t45 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t6, 0, t7, -t6, -t7 * pkin(4) - t6 * qJ(5) + t14 * qJD(5), 0, 0, 0, 0, 0, t52, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t74, 0, -t30, -t31, -t12, t30, t31, t12 * t49, -t43 * t63 + t61, -0.4e1 * t43 * t67 - t74 * t87, -t16, -t18, 0, t50 * t42 - t94 * t45, t94 * t42 + t50 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t74, 0, t34, t74, t12, 0, 0, 0, 0, 0, t17, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, qJ(5) * t91, -0.2e1 * t67, 0.2e1 * t63, 0, 0, 0, 0.2e1 * qJD(5) * t42 + 0.2e1 * t45 * t72, 0.2e1 * qJD(5) * t45 - 0.2e1 * t42 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, 0, 0, t30, 0, 0, 0, 0, 0, -t16, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t15, -t34, t3, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, -t79, 0, -t42 * t77, -t45 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, -t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t1;
