% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:44:58
% EndTime: 2019-03-09 01:44:59
% DurationCPUTime: 0.53s
% Computational Cost: add. (650->98), mult. (1328->175), div. (0->0), fcn. (1225->8), ass. (0->67)
t44 = sin(qJ(4));
t46 = cos(qJ(4));
t78 = sin(pkin(10));
t79 = cos(pkin(10));
t49 = t44 * t79 + t46 * t78;
t23 = t49 * qJD(4);
t28 = -t78 * t44 + t79 * t46;
t17 = t28 * t23;
t22 = t28 * qJD(4);
t84 = t49 * t22;
t93 = 0.2e1 * t17 - 0.2e1 * t84;
t92 = (-t22 * t78 + t23 * t79) * pkin(4);
t25 = t28 ^ 2;
t90 = 2 * qJD(3);
t89 = 2 * qJD(6);
t59 = -cos(pkin(9)) * pkin(1) - pkin(2) - pkin(7);
t55 = qJ(5) - t59;
t51 = t55 * t46;
t19 = -qJD(4) * t51 - t44 * qJD(5);
t75 = t44 * qJD(4);
t47 = -t46 * qJD(5) + t55 * t75;
t3 = t19 * t78 - t47 * t79;
t43 = sin(qJ(6));
t88 = t3 * t43;
t45 = cos(qJ(6));
t87 = t3 * t45;
t86 = t23 * t43;
t85 = t23 * t45;
t15 = t49 * t23;
t83 = t43 * t22;
t82 = t45 * t22;
t81 = -t45 * t15 + t28 * t82;
t41 = t45 ^ 2;
t80 = t43 ^ 2 - t41;
t31 = sin(pkin(9)) * pkin(1) + qJ(3);
t77 = qJD(6) * t43;
t76 = qJD(6) * t45;
t74 = t46 * qJD(4);
t29 = pkin(4) * t74 + qJD(3);
t73 = t43 * t85;
t36 = -pkin(4) * t79 - pkin(5);
t72 = t36 * t89;
t71 = t44 * pkin(4) + t31;
t70 = t28 * t77;
t69 = t43 * t76;
t68 = t49 ^ 2 + t25;
t65 = t80 * qJD(6);
t6 = pkin(5) * t49 - t28 * pkin(8) + t71;
t24 = t55 * t44;
t8 = -t24 * t79 - t51 * t78;
t64 = t43 * t8 - t45 * t6;
t63 = t43 * t6 + t45 * t8;
t35 = pkin(4) * t78 + pkin(8);
t61 = -t22 * t35 - t23 * t36;
t60 = t28 * t36 - t35 * t49;
t57 = -t49 * t76 - t83;
t56 = t28 * t76 - t86;
t10 = t70 + t85;
t52 = t59 * qJD(4);
t4 = t19 * t79 + t47 * t78;
t7 = -t24 * t78 + t51 * t79;
t48 = -t8 * t22 - t7 * t23 + t3 * t28 - t4 * t49;
t9 = t49 * t77 - t82;
t5 = t22 * pkin(5) + t23 * pkin(8) + t29;
t2 = -qJD(6) * t63 - t43 * t4 + t45 * t5;
t1 = qJD(6) * t64 - t45 * t4 - t43 * t5;
t11 = [0, 0, 0, 0, 0, t90, t31 * t90, -0.2e1 * t44 * t74, 0.2e1 * (t44 ^ 2 - t46 ^ 2) * qJD(4), 0, 0, 0, 0.2e1 * qJD(3) * t44 + 0.2e1 * t31 * t74, 0.2e1 * qJD(3) * t46 - 0.2e1 * t31 * t75, 0.2e1 * t48, 0.2e1 * t29 * t71 + 0.2e1 * t7 * t3 + 0.2e1 * t8 * t4, -0.2e1 * t17 * t41 - 0.2e1 * t25 * t69, t25 * t80 * t89 + 0.4e1 * t28 * t73, -0.2e1 * t49 * t70 + 0.2e1 * t81, -0.2e1 * t28 * t83 - 0.2e1 * t49 * t56, 0.2e1 * t84, 0.2e1 * t2 * t49 - 0.2e1 * t64 * t22 - 0.2e1 * t7 * t86 + 0.2e1 * (t7 * t76 + t88) * t28, 0.2e1 * t1 * t49 - 0.2e1 * t63 * t22 - 0.2e1 * t7 * t85 + 0.2e1 * (-t7 * t77 + t87) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * t22 - t8 * t23 + t4 * t28 + t3 * t49, 0, 0, 0, 0, 0, 0, t45 * (-t28 * t22 + t15) + t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, -t48, 0, 0, 0, 0, 0, t43 * t93 - t68 * t76, t45 * t93 + t68 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t74, 0, -t44 * t52, -t46 * t52, t92 (-t3 * t79 + t4 * t78) * pkin(4), -t28 * t65 - t73, t23 * t80 - 0.4e1 * t28 * t69, -t57, -t9, 0, -t87 + t61 * t43 + (t43 * t7 + t45 * t60) * qJD(6), t88 + t61 * t45 + (-t43 * t60 + t45 * t7) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, t75, 0 (-t22 * t79 - t23 * t78) * pkin(4), 0, 0, 0, 0, 0, t9, -t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t74, 0, -t92, 0, 0, 0, 0, 0, -t10, -t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t69, -0.2e1 * t65, 0, 0, 0, t43 * t72, t45 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, -t9, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t56, t22, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t77, 0, -t35 * t76, t35 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t11;
