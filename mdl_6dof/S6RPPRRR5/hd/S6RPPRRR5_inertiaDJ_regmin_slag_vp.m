% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRRR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:29:00
% EndTime: 2019-03-09 02:29:03
% DurationCPUTime: 0.65s
% Computational Cost: add. (721->107), mult. (1559->194), div. (0->0), fcn. (1390->6), ass. (0->76)
t93 = cos(qJ(5));
t70 = t93 * qJD(5);
t97 = t93 * qJD(4) + t70;
t92 = sin(qJ(5));
t72 = qJD(5) * t92;
t96 = -qJD(4) * t92 - t72;
t56 = cos(qJ(6));
t51 = t56 ^ 2;
t54 = sin(qJ(6));
t85 = t54 ^ 2 - t51;
t69 = qJD(6) * t85;
t55 = sin(qJ(4));
t57 = cos(qJ(4));
t27 = t92 * t55 - t93 * t57;
t25 = t27 ^ 2;
t52 = -pkin(7) + qJ(2);
t95 = pkin(8) - t52;
t30 = t95 * t55;
t31 = t95 * t57;
t20 = -t93 * t30 - t92 * t31;
t81 = t57 * qJD(2);
t82 = t55 * qJD(4);
t60 = t95 * t82 + t81;
t83 = t55 * qJD(2);
t61 = -qJD(4) * t31 + t83;
t6 = t20 * qJD(5) - t93 * t60 + t92 * t61;
t4 = t6 * t54;
t19 = -t92 * t30 + t93 * t31;
t48 = qJD(6) * t56;
t94 = t19 * t48 + t4;
t17 = -t97 * t55 + t96 * t57;
t91 = t27 * t17;
t90 = t27 * t56;
t18 = t96 * t55 + t97 * t57;
t28 = t93 * t55 + t92 * t57;
t89 = t28 * t18;
t88 = t54 * t17;
t87 = t56 * t17;
t53 = pkin(1) + qJ(3);
t46 = -t93 * pkin(4) - pkin(5);
t67 = pkin(4) * t72;
t86 = t46 * t48 + t54 * t67;
t84 = qJD(6) * t54;
t80 = t57 * qJD(4);
t36 = pkin(4) * t80 + qJD(3);
t79 = qJ(2) * qJD(2);
t43 = t55 * pkin(4) + t53;
t78 = pkin(5) * t84;
t77 = pkin(5) * t48;
t76 = t54 * t48;
t75 = t28 ^ 2 + t25;
t74 = 0.4e1 * t54 * t90;
t68 = pkin(4) * t70;
t13 = t28 * pkin(5) + t27 * pkin(9) + t43;
t66 = t56 * t13 - t54 * t20;
t65 = t54 * t13 + t56 * t20;
t45 = t92 * pkin(4) + pkin(9);
t64 = t27 * t46 + t28 * t45;
t11 = t27 * t48 - t88;
t12 = t27 * t84 + t87;
t10 = t54 * t18 + t28 * t48;
t9 = -t56 * t18 + t28 * t84;
t63 = t46 * t84 - t56 * t67;
t62 = -0.2e1 * t89 + 0.2e1 * t91;
t59 = t17 * t46 - t18 * t45 + (-t92 * t27 - t93 * t28) * qJD(5) * pkin(4);
t58 = 0.2e1 * qJD(2);
t35 = 0.2e1 * t76;
t26 = -0.2e1 * t69;
t14 = t19 * t84;
t8 = t18 * pkin(5) - t17 * pkin(9) + t36;
t7 = t27 * t69 + t54 * t87;
t5 = -t30 * t72 + t31 * t70 - t92 * t60 - t93 * t61;
t3 = qJD(6) * t74 - t85 * t17;
t2 = -t65 * qJD(6) + t54 * t5 + t56 * t8;
t1 = -t66 * qJD(6) + t56 * t5 - t54 * t8;
t15 = [0, 0, 0, 0, t58, 0.2e1 * t79, t58, 0.2e1 * qJD(3), 0.2e1 * t53 * qJD(3) + 0.2e1 * t79, -0.2e1 * t55 * t80, 0.2e1 * (t55 ^ 2 - t57 ^ 2) * qJD(4), 0, 0, 0, 0.2e1 * qJD(3) * t55 + 0.2e1 * t53 * t80, 0.2e1 * qJD(3) * t57 - 0.2e1 * t53 * t82, -0.2e1 * t91, -0.2e1 * t17 * t28 + 0.2e1 * t27 * t18, 0, 0, 0, 0.2e1 * t43 * t18 + 0.2e1 * t36 * t28, 0.2e1 * t43 * t17 - 0.2e1 * t36 * t27, -0.2e1 * t25 * t76 - 0.2e1 * t51 * t91, t17 * t74 + 0.2e1 * t25 * t69, 0.2e1 * t27 * t9 + 0.2e1 * t28 * t87, 0.2e1 * t10 * t27 - 0.2e1 * t28 * t88, 0.2e1 * t89, -0.2e1 * t11 * t19 + 0.2e1 * t18 * t66 + 0.2e1 * t2 * t28 - 0.2e1 * t27 * t4, 0.2e1 * t1 * t28 + 0.2e1 * t12 * t19 - 0.2e1 * t18 * t65 - 0.2e1 * t6 * t90; 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), 0, 0, 0, 0, 0, -t80, t82, 0, 0, 0, 0, 0, -t18, -t17, 0, 0, 0, 0, 0, t9, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75 * t48 + t54 * t62, t56 * t62 + t75 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, -t80, 0, -t52 * t82 + t81, -t52 * t80 - t83, 0, 0, t17, -t18, 0, -t6, t5, t7, t3, t10, -t9, 0, t14 + (-qJD(6) * t64 - t6) * t56 + t59 * t54, t56 * t59 + t64 * t84 + t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, -t80, 0, 0, 0, 0, 0, t17, -t18, 0, 0, 0, 0, 0, t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t67, -0.2e1 * t68, t35, t26, 0, 0, 0, 0.2e1 * t63, 0.2e1 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t18, 0, -t6, t5, t7, t3, t10, -t9, 0, t14 + (-pkin(5) * t17 - pkin(9) * t18) * t54 + (-t6 + (pkin(5) * t27 - pkin(9) * t28) * qJD(6)) * t56, -pkin(5) * t12 + pkin(9) * t9 + t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t18, 0, 0, 0, 0, 0, t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t68, t35, t26, 0, 0, 0, t63 - t78, -t77 + t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t26, 0, 0, 0, -0.2e1 * t78, -0.2e1 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t11, t18, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t84, 0, -t45 * t48 - t54 * t68, t45 * t84 - t56 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t84, 0, -pkin(9) * t48, pkin(9) * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t15;
