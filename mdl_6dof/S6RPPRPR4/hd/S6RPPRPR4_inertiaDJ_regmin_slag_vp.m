% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRPR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:47:08
% EndTime: 2019-03-09 01:47:11
% DurationCPUTime: 0.82s
% Computational Cost: add. (812->109), mult. (1657->229), div. (0->0), fcn. (1559->8), ass. (0->77)
t44 = sin(qJ(4));
t46 = cos(qJ(4));
t78 = sin(pkin(10));
t79 = cos(pkin(10));
t25 = t44 * t79 + t46 * t78;
t91 = 0.4e1 * t25;
t24 = t44 * t78 - t46 * t79;
t22 = t24 * qJD(4);
t45 = cos(qJ(6));
t40 = t45 ^ 2;
t43 = sin(qJ(6));
t80 = t43 ^ 2 - t40;
t65 = qJD(6) * t80;
t21 = t25 * qJD(4);
t90 = 2 * qJD(2);
t41 = sin(pkin(9));
t42 = cos(pkin(9));
t47 = -pkin(1) - pkin(2);
t81 = qJ(2) * t42 + t41 * t47;
t72 = -pkin(7) + t81;
t64 = -qJ(5) + t72;
t51 = qJD(4) * t64;
t75 = t42 * qJD(2);
t63 = -qJD(5) + t75;
t14 = -t44 * t51 + t46 * t63;
t48 = -t44 * t63 - t46 * t51;
t3 = t14 * t78 - t48 * t79;
t89 = t3 * t43;
t88 = t3 * t45;
t87 = t24 * t21;
t86 = t25 * t22;
t85 = t25 * t43;
t84 = t25 * t45;
t83 = t43 * t22;
t82 = t45 * t22;
t77 = qJD(6) * t43;
t76 = qJD(6) * t45;
t35 = t41 * qJD(2);
t74 = t44 * qJD(4);
t73 = t46 * qJD(4);
t71 = t43 * t82;
t33 = -pkin(4) * t79 - pkin(5);
t70 = 0.2e1 * qJD(6) * t33;
t69 = t43 * t76;
t68 = -qJ(2) * t41 + t42 * t47;
t27 = pkin(3) - t68;
t62 = t72 * qJD(4);
t28 = -pkin(4) * t74 + t35;
t53 = pkin(4) * t46 + t27;
t10 = -pkin(5) * t24 + pkin(8) * t25 + t53;
t18 = t64 * t46;
t52 = t64 * t44;
t8 = t18 * t79 - t52 * t78;
t59 = t10 * t45 - t43 * t8;
t58 = t10 * t43 + t45 * t8;
t20 = t24 * t41;
t57 = -t20 * t45 - t42 * t43;
t56 = -t20 * t43 + t42 * t45;
t32 = pkin(4) * t78 + pkin(8);
t55 = t21 * t32 + t22 * t33;
t54 = t24 * t32 - t25 * t33;
t13 = t21 * t43 + t24 * t76;
t12 = -t21 * t45 + t24 * t77;
t50 = t25 * t76 - t83;
t49 = t25 * t77 + t82;
t23 = t25 ^ 2;
t19 = t25 * t41;
t17 = t41 * t21;
t16 = t22 * t41;
t9 = -pkin(5) * t21 - pkin(8) * t22 + t28;
t7 = t18 * t78 + t52 * t79;
t6 = -qJD(6) * t57 + t43 * t17;
t5 = qJD(6) * t56 + t45 * t17;
t4 = t14 * t79 + t48 * t78;
t2 = -qJD(6) * t58 - t43 * t4 + t45 * t9;
t1 = -qJD(6) * t59 - t45 * t4 - t43 * t9;
t11 = [0, 0, 0, 0, t90, qJ(2) * t90, 0.2e1 * t35, 0.2e1 * t75 (-t41 * t68 + t42 * t81) * t90, 0.2e1 * t44 * t73, 0.2e1 * (-t44 ^ 2 + t46 ^ 2) * qJD(4), 0, 0, 0, -0.2e1 * t27 * t74 + 0.2e1 * t35 * t46, -0.2e1 * t27 * t73 - 0.2e1 * t35 * t44, 0.2e1 * t21 * t8 + 0.2e1 * t22 * t7 + 0.2e1 * t24 * t4 - 0.2e1 * t25 * t3, 0.2e1 * t28 * t53 + 0.2e1 * t3 * t7 + 0.2e1 * t4 * t8, -0.2e1 * t23 * t69 - 0.2e1 * t40 * t86, 0.2e1 * t23 * t65 + t71 * t91, 0.2e1 * t21 * t84 - 0.2e1 * t24 * t49, -0.2e1 * t21 * t85 - 0.2e1 * t24 * t50, 0.2e1 * t87, -0.2e1 * t2 * t24 - 0.2e1 * t59 * t21 + 0.2e1 * t7 * t83 + 0.2e1 * (-t7 * t76 - t89) * t25, -0.2e1 * t1 * t24 + 0.2e1 * t58 * t21 + 0.2e1 * t7 * t82 + 0.2e1 * (t7 * t77 - t88) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t74, t42 * t73, t16 * t25 - t17 * t24 + t19 * t22 - t20 * t21, -t16 * t7 - t17 * t8 + t19 * t3 - t20 * t4 - t28 * t42, 0, 0, 0, 0, 0, t16 * t85 - t19 * t50 + t21 * t56 - t6 * t24, t16 * t84 + t19 * t49 + t21 * t57 - t5 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t16 * t19 + 0.2e1 * t17 * t20, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t7 - t22 * t8 + t24 * t3 + t25 * t4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16 * t24 - t17 * t25 + t19 * t21 + t20 * t22, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t86 + 0.2e1 * t87, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, t74, 0, -t44 * t75 - t46 * t62, t44 * t62 - t46 * t75 (t21 * t78 - t22 * t79) * pkin(4) (-t3 * t79 + t4 * t78) * pkin(4), t25 * t65 + t71, -t22 * t80 + t69 * t91, -t13, t12, 0, -t88 + t55 * t43 + (t43 * t7 + t45 * t54) * qJD(6), t89 + t55 * t45 + (-t43 * t54 + t45 * t7) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41 * t73, t41 * t74, 0 (t16 * t79 - t17 * t78) * pkin(4), 0, 0, 0, 0, 0, t16 * t45 + t19 * t77, -t16 * t43 + t19 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, -t73, 0 (-t21 * t79 - t22 * t78) * pkin(4), 0, 0, 0, 0, 0, t12, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t69, -0.2e1 * t65, 0, 0, 0, t43 * t70, t45 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, t12, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t50, -t21, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t77, 0, -t32 * t76, t32 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t11;
