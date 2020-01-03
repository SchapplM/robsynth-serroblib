% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x27]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRR12_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:13:08
% EndTime: 2019-12-31 19:13:10
% DurationCPUTime: 0.54s
% Computational Cost: add. (663->99), mult. (1477->187), div. (0->0), fcn. (1312->6), ass. (0->74)
t53 = sin(qJ(3));
t55 = cos(qJ(3));
t90 = cos(qJ(4));
t70 = t90 * qJD(3);
t89 = sin(qJ(4));
t72 = qJD(3) * t89;
t95 = -t53 * t70 - t55 * t72;
t69 = t90 * qJD(4);
t94 = t70 + t69;
t54 = cos(qJ(5));
t51 = t54 ^ 2;
t52 = sin(qJ(5));
t82 = t52 ^ 2 - t51;
t68 = qJD(5) * t82;
t27 = t53 * t89 - t55 * t90;
t25 = t27 ^ 2;
t93 = 2 * qJD(2);
t56 = -pkin(1) - pkin(6);
t92 = pkin(7) - t56;
t30 = t92 * t53;
t31 = t92 * t55;
t20 = -t30 * t90 - t31 * t89;
t8 = t20 * qJD(4) + t95 * t92;
t6 = t8 * t52;
t19 = -t89 * t30 + t90 * t31;
t48 = qJD(5) * t54;
t91 = t19 * t48 + t6;
t71 = qJD(4) * t89;
t17 = -t53 * t69 - t55 * t71 + t95;
t88 = t27 * t17;
t87 = t27 * t54;
t64 = t53 * t72;
t18 = -t53 * t71 + t94 * t55 - t64;
t28 = t53 * t90 + t55 * t89;
t86 = t28 * t18;
t85 = t52 * t17;
t84 = t54 * t17;
t46 = -pkin(3) * t90 - pkin(4);
t66 = pkin(3) * t71;
t83 = t46 * t48 + t52 * t66;
t43 = t53 * pkin(3) + qJ(2);
t81 = qJD(5) * t52;
t80 = t53 * qJD(3);
t79 = t55 * qJD(3);
t36 = pkin(3) * t79 + qJD(2);
t78 = qJ(2) * qJD(3);
t77 = pkin(4) * t81;
t76 = pkin(4) * t48;
t75 = t52 * t48;
t74 = t28 ^ 2 + t25;
t73 = 0.4e1 * t52 * t87;
t67 = pkin(3) * t69;
t13 = t28 * pkin(4) + t27 * pkin(8) + t43;
t62 = t54 * t13 - t52 * t20;
t61 = t52 * t13 + t54 * t20;
t45 = pkin(3) * t89 + pkin(8);
t60 = t27 * t46 + t28 * t45;
t11 = t27 * t48 - t85;
t12 = t27 * t81 + t84;
t10 = t52 * t18 + t28 * t48;
t9 = -t54 * t18 + t28 * t81;
t59 = t46 * t81 - t54 * t66;
t58 = -0.2e1 * t86 + 0.2e1 * t88;
t57 = t17 * t46 - t18 * t45 + (-t27 * t89 - t28 * t90) * qJD(4) * pkin(3);
t35 = 0.2e1 * t75;
t26 = -0.2e1 * t68;
t14 = t19 * t81;
t7 = -t30 * t71 + t94 * t31 - t64 * t92;
t5 = t18 * pkin(4) - t17 * pkin(8) + t36;
t4 = t27 * t68 + t52 * t84;
t3 = qJD(5) * t73 - t17 * t82;
t2 = -qJD(5) * t61 + t54 * t5 + t52 * t7;
t1 = -qJD(5) * t62 - t52 * t5 + t54 * t7;
t15 = [0, 0, 0, 0, t93, qJ(2) * t93, -0.2e1 * t53 * t79, 0.2e1 * (t53 ^ 2 - t55 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t53 + 0.2e1 * t55 * t78, 0.2e1 * qJD(2) * t55 - 0.2e1 * t53 * t78, -0.2e1 * t88, -0.2e1 * t17 * t28 + 0.2e1 * t27 * t18, 0, 0, 0, 0.2e1 * t43 * t18 + 0.2e1 * t36 * t28, 0.2e1 * t43 * t17 - 0.2e1 * t36 * t27, -0.2e1 * t25 * t75 - 0.2e1 * t51 * t88, t17 * t73 + 0.2e1 * t25 * t68, 0.2e1 * t27 * t9 + 0.2e1 * t28 * t84, 0.2e1 * t10 * t27 - 0.2e1 * t28 * t85, 0.2e1 * t86, -0.2e1 * t11 * t19 + 0.2e1 * t18 * t62 + 0.2e1 * t2 * t28 - 0.2e1 * t27 * t6, 0.2e1 * t1 * t28 + 0.2e1 * t12 * t19 - 0.2e1 * t18 * t61 - 0.2e1 * t8 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48 * t74 + t52 * t58, t54 * t58 + t74 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t80, -t79, 0, -t56 * t80, -t56 * t79, 0, 0, t17, -t18, 0, -t8, t7, t4, t3, t10, -t9, 0, t14 + (-qJD(5) * t60 - t8) * t54 + t57 * t52, t54 * t57 + t60 * t81 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, -t79, 0, 0, 0, 0, 0, t17, -t18, 0, 0, 0, 0, 0, t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t66, -0.2e1 * t67, t35, t26, 0, 0, 0, 0.2e1 * t59, 0.2e1 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t18, 0, -t8, t7, t4, t3, t10, -t9, 0, t14 + (-pkin(4) * t17 - pkin(8) * t18) * t52 + (-t8 + (pkin(4) * t27 - pkin(8) * t28) * qJD(5)) * t54, -pkin(4) * t12 + pkin(8) * t9 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t18, 0, 0, 0, 0, 0, t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t67, t35, t26, 0, 0, 0, t59 - t77, -t76 + t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t26, 0, 0, 0, -0.2e1 * t77, -0.2e1 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t11, t18, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t81, 0, -t45 * t48 - t52 * t67, t45 * t81 - t54 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t81, 0, -pkin(8) * t48, pkin(8) * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t15;
