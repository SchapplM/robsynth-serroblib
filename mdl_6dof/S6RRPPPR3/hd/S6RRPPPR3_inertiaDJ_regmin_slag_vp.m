% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:15:38
% EndTime: 2019-03-09 08:15:41
% DurationCPUTime: 0.86s
% Computational Cost: add. (742->137), mult. (1634->256), div. (0->0), fcn. (1273->6), ass. (0->84)
t64 = sin(pkin(9));
t65 = cos(pkin(9));
t67 = sin(qJ(6));
t69 = cos(qJ(6));
t104 = -t67 * t64 + t69 * t65;
t33 = t104 * qJD(6);
t68 = sin(qJ(2));
t56 = t68 * qJ(3);
t70 = cos(qJ(2));
t105 = -t70 * pkin(2) - t56;
t94 = t64 ^ 2 + t65 ^ 2;
t41 = t94 * qJD(5);
t72 = 2 * qJD(3);
t71 = -pkin(2) - pkin(3);
t61 = -qJ(5) + t71;
t103 = pkin(8) - t61;
t89 = qJ(4) * qJD(2);
t81 = -t70 * qJD(4) + t68 * t89;
t91 = t68 * qJD(2);
t87 = pkin(7) * t91;
t29 = -t81 + t87;
t57 = t70 * qJ(4);
t47 = t70 * pkin(7) - t57;
t102 = t47 * t29;
t101 = t47 * t68;
t100 = t64 * t70;
t99 = t65 * t68;
t98 = t65 * t70;
t66 = qJ(3) + pkin(4);
t55 = t70 * qJD(2);
t95 = qJ(3) * t55 + t68 * qJD(3);
t18 = t70 * qJD(5) + (pkin(4) * t70 + t61 * t68) * qJD(2) + t95;
t53 = pkin(7) * t55;
t30 = -t68 * qJD(4) - t70 * t89 + t53;
t7 = t64 * t18 + t65 * t30;
t42 = -pkin(1) + t105;
t35 = t70 * pkin(3) - t42;
t24 = t68 * pkin(4) + t70 * qJ(5) + t35;
t46 = (pkin(7) - qJ(4)) * t68;
t17 = t64 * t24 + t65 * t46;
t92 = t47 * qJD(3);
t90 = t70 * qJD(3);
t88 = -0.2e1 * pkin(1) * qJD(2);
t86 = t64 * t91;
t85 = t64 * t55;
t84 = t65 * t91;
t83 = t65 * t55;
t82 = pkin(5) * t64 - pkin(7);
t6 = t65 * t18 - t64 * t30;
t16 = t65 * t24 - t64 * t46;
t80 = t6 * t65 + t7 * t64;
t3 = -t6 * t64 + t7 * t65;
t11 = pkin(8) * t100 + t17;
t8 = t68 * pkin(5) + pkin(8) * t98 + t16;
t79 = t69 * t11 + t67 * t8;
t78 = t67 * t11 - t69 * t8;
t39 = t103 * t64;
t40 = t103 * t65;
t77 = t69 * t39 + t67 * t40;
t76 = t67 * t39 - t69 * t40;
t37 = t69 * t64 + t67 * t65;
t19 = -t33 * t68 - t37 * t55;
t34 = t37 * qJD(6);
t75 = -t104 * t55 + t34 * t68;
t74 = t105 * qJD(2) + t90;
t73 = -t90 + qJD(5) * t68 + (-t61 * t70 + t66 * t68) * qJD(2);
t60 = qJ(3) * t72;
t49 = t65 * pkin(5) + t66;
t48 = 0.2e1 * t68 * t55;
t32 = -t82 * t70 - t57;
t31 = pkin(2) * t91 - t95;
t28 = t104 * t70;
t27 = t37 * t70;
t25 = t71 * t91 + t95;
t23 = t82 * t91 + t81;
t15 = t70 * t33 - t67 * t84 - t69 * t86;
t14 = t104 * t91 + t70 * t34;
t10 = t37 * qJD(5) - t76 * qJD(6);
t9 = qJD(5) * t104 - t77 * qJD(6);
t5 = -pkin(8) * t86 + t7;
t4 = (pkin(5) * t70 - pkin(8) * t99) * qJD(2) + t6;
t2 = -t79 * qJD(6) + t69 * t4 - t67 * t5;
t1 = t78 * qJD(6) - t67 * t4 - t69 * t5;
t12 = [0, 0, 0, t48, 0.2e1 * (-t68 ^ 2 + t70 ^ 2) * qJD(2), 0, 0, 0, t68 * t88, t70 * t88, -0.2e1 * t31 * t70 + 0.2e1 * t42 * t91, 0, -0.2e1 * t31 * t68 - 0.2e1 * t42 * t55, 0.2e1 * t42 * t31, 0.2e1 * t25 * t68 + 0.2e1 * t35 * t55, -0.2e1 * t25 * t70 + 0.2e1 * t35 * t91, 0.2e1 * t29 * t70 - 0.2e1 * t30 * t68 + 0.2e1 * (-t46 * t70 + t101) * qJD(2), 0.2e1 * t35 * t25 + 0.2e1 * t46 * t30 - 0.2e1 * t102, 0.2e1 * t29 * t100 + 0.2e1 * t6 * t68 + 0.2e1 * (t64 * t101 + t16 * t70) * qJD(2), 0.2e1 * t29 * t98 - 0.2e1 * t7 * t68 + 0.2e1 * (-t17 * t70 + t47 * t99) * qJD(2), 0.2e1 * t80 * t70 + 0.2e1 * (-t16 * t65 - t17 * t64) * t91, 0.2e1 * t16 * t6 + 0.2e1 * t17 * t7 - 0.2e1 * t102, -0.2e1 * t28 * t14, 0.2e1 * t14 * t27 - 0.2e1 * t28 * t15, 0.2e1 * t14 * t68 - 0.2e1 * t28 * t55, 0.2e1 * t15 * t68 + 0.2e1 * t27 * t55, t48, -0.2e1 * t32 * t15 + 0.2e1 * t2 * t68 - 0.2e1 * t23 * t27 - 0.2e1 * t78 * t55, 0.2e1 * t1 * t68 + 0.2e1 * t32 * t14 - 0.2e1 * t23 * t28 - 0.2e1 * t79 * t55; 0, 0, 0, 0, 0, t55, -t91, 0, -t53, t87, -t53, t74, -t87, t74 * pkin(7), -t29, t30, -t90 + (-t70 * t71 + t56) * qJD(2), -t29 * qJ(3) + t30 * t71 + t92, -t29 * t65 + t73 * t64, t29 * t64 + t73 * t65, -t3, t92 - t29 * t66 + t3 * t61 + (t16 * t64 - t17 * t65) * qJD(5), -t14 * t37 + t28 * t33, -t104 * t14 - t37 * t15 - t33 * t27 - t28 * t34, t19, t75, 0, -qJD(3) * t27 + t10 * t68 + t104 * t23 - t49 * t15 - t32 * t34 + t77 * t55, -qJD(3) * t28 + t49 * t14 - t23 * t37 - t32 * t33 - t76 * t55 + t9 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t60, t72, 0, 0, t60, t65 * t72, -0.2e1 * qJD(3) * t64, 0.2e1 * t41, 0.2e1 * t66 * qJD(3) - 0.2e1 * t61 * t41, 0.2e1 * t37 * t33, 0.2e1 * t104 * t33 - 0.2e1 * t37 * t34, 0, 0, 0, 0.2e1 * qJD(3) * t104 - 0.2e1 * t49 * t34, -0.2e1 * qJD(3) * t37 - 0.2e1 * t49 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, t53, 0, 0, -t55, t30, -t85, -t83, 0, t3, 0, 0, 0, 0, 0, t19, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t91, 0, t25, t83, -t85, -t94 * t91, t80, 0, 0, 0, 0, 0, -t75, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t84, 0, -t29, 0, 0, 0, 0, 0, -t15, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), 0, 0, 0, 0, 0, -t34, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t15, t55, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t34, 0, t10, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t12;
