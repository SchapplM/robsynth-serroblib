% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:12
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:10:44
% EndTime: 2021-01-15 23:10:50
% DurationCPUTime: 0.95s
% Computational Cost: add. (1711->123), mult. (3995->238), div. (0->0), fcn. (3714->8), ass. (0->90)
t110 = qJD(2) + qJD(3);
t109 = pkin(6) + pkin(7);
t62 = cos(qJ(2));
t111 = t109 * t62;
t58 = sin(qJ(3));
t61 = cos(qJ(3));
t59 = sin(qJ(2));
t90 = t109 * t59;
t66 = t111 * t58 + t61 * t90;
t21 = t110 * (-t111 * t61 + t58 * t90);
t71 = t58 * t59 - t61 * t62;
t29 = t110 * t71;
t43 = t58 * t62 + t61 * t59;
t112 = t29 * qJ(4) - t43 * qJD(4) + t21;
t60 = cos(qJ(5));
t55 = t60 ^ 2;
t57 = sin(qJ(5));
t100 = t57 ^ 2 - t55;
t81 = t100 * qJD(5);
t20 = t110 * t66;
t30 = t110 * t43;
t12 = -t30 * qJ(4) - t71 * qJD(4) - t20;
t56 = sin(pkin(9));
t98 = cos(pkin(9));
t4 = -t98 * t112 + t56 * t12;
t3 = t4 * t57;
t23 = t71 * (-qJ(4) - t109);
t63 = -t43 * qJ(4) - t66;
t16 = t56 * t23 - t98 * t63;
t53 = qJD(5) * t60;
t108 = t16 * t53 + t3;
t28 = t98 * t43 - t56 * t71;
t107 = t28 * t60;
t18 = -t56 * t29 + t98 * t30;
t106 = t57 * t18;
t104 = t60 * t18;
t19 = -t98 * t29 - t56 * t30;
t103 = t60 * t19;
t50 = t61 * pkin(2) + pkin(3);
t39 = -t56 * t58 * pkin(2) + t98 * t50;
t35 = -pkin(4) - t39;
t82 = t98 * t58;
t99 = pkin(2) * qJD(3);
t37 = (t56 * t61 + t82) * t99;
t101 = t35 * t53 + t37 * t57;
t40 = pkin(2) * t82 + t56 * t50;
t97 = qJD(5) * t57;
t96 = t59 * qJD(2);
t95 = t62 * qJD(2);
t94 = -0.2e1 * pkin(1) * qJD(2);
t52 = pkin(2) * t96;
t92 = t58 * t99;
t91 = t61 * t99;
t49 = -t98 * pkin(3) - pkin(4);
t87 = t49 * t97;
t86 = t49 * t53;
t85 = t57 * t53;
t51 = -t62 * pkin(2) - pkin(1);
t24 = t30 * pkin(3) + t52;
t84 = -0.4e1 * t57 * t107;
t83 = t35 * t97 - t37 * t60;
t27 = t56 * t43 + t98 * t71;
t34 = t71 * pkin(3) + t51;
t15 = t27 * pkin(4) - t28 * pkin(8) + t34;
t17 = t98 * t23 + t56 * t63;
t77 = t60 * t15 - t57 * t17;
t76 = t57 * t15 + t60 * t17;
t48 = t56 * pkin(3) + pkin(8);
t75 = -t18 * t48 + t19 * t49;
t36 = pkin(8) + t40;
t74 = t27 * t36 - t28 * t35;
t38 = -t56 * t92 + t98 * t91;
t73 = -t38 * t27 + t37 * t28;
t72 = t27 * t48 - t28 * t49;
t10 = t27 * t53 + t106;
t69 = t57 * t19 + t28 * t53;
t68 = -t28 * t97 + t103;
t65 = -t18 * t36 + t19 * t35 + t73;
t46 = 0.2e1 * t85;
t42 = -0.2e1 * t81;
t25 = t28 ^ 2;
t13 = t16 * t97;
t9 = -t27 * t97 + t104;
t8 = t57 * t103 - t28 * t81;
t7 = t18 * pkin(4) - t19 * pkin(8) + t24;
t6 = qJD(5) * t84 - t100 * t19;
t5 = t112 * t56 + t98 * t12;
t2 = -t76 * qJD(5) - t57 * t5 + t60 * t7;
t1 = -t77 * qJD(5) - t60 * t5 - t57 * t7;
t11 = [0, 0, 0, 0.2e1 * t59 * t95, 0.2e1 * (-t59 ^ 2 + t62 ^ 2) * qJD(2), 0, 0, 0, t59 * t94, t62 * t94, -0.2e1 * t43 * t29, 0.2e1 * t29 * t71 - 0.2e1 * t43 * t30, 0, 0, 0, 0.2e1 * t51 * t30 + 0.2e1 * t71 * t52, -0.2e1 * t51 * t29 + 0.2e1 * t43 * t52, 0.2e1 * t34 * t18 + 0.2e1 * t24 * t27, 0.2e1 * t34 * t19 + 0.2e1 * t24 * t28, 0.2e1 * t16 * t19 - 0.2e1 * t17 * t18 - 0.2e1 * t5 * t27 + 0.2e1 * t4 * t28, 0.2e1 * t16 * t4 + 0.2e1 * t17 * t5 + 0.2e1 * t34 * t24, 0.2e1 * t55 * t28 * t19 - 0.2e1 * t25 * t85, t19 * t84 + 0.2e1 * t25 * t81, 0.2e1 * t28 * t104 + 0.2e1 * t68 * t27, -0.2e1 * t28 * t106 - 0.2e1 * t69 * t27, 0.2e1 * t27 * t18, 0.2e1 * t69 * t16 + 0.2e1 * t77 * t18 + 0.2e1 * t2 * t27 + 0.2e1 * t28 * t3, 0.2e1 * t1 * t27 + 0.2e1 * t4 * t107 + 0.2e1 * t68 * t16 - 0.2e1 * t76 * t18; 0, 0, 0, 0, 0, t95, -t96, 0, -pkin(6) * t95, pkin(6) * t96, 0, 0, -t29, -t30, 0, t21, t20, -t4, -t5, -t40 * t18 - t39 * t19 + t73, t16 * t37 + t17 * t38 - t4 * t39 + t5 * t40, t8, t6, t10, t9, 0, t13 + (-t74 * qJD(5) - t4) * t60 + t65 * t57, t65 * t60 + t74 * t97 + t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t92, -0.2e1 * t91, -0.2e1 * t37, -0.2e1 * t38, 0, -0.2e1 * t39 * t37 + 0.2e1 * t40 * t38, t46, t42, 0, 0, 0, 0.2e1 * t83, 0.2e1 * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t30, 0, t21, t20, -t4, -t5, (-t18 * t56 - t98 * t19) * pkin(3), (-t98 * t4 + t5 * t56) * pkin(3), t8, t6, t10, t9, 0, t13 + t75 * t57 + (-t72 * qJD(5) - t4) * t60, t75 * t60 + t72 * t97 + t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, -t91, -t37, -t38, 0, (-t98 * t37 + t38 * t56) * pkin(3), t46, t42, 0, 0, 0, t83 + t87, t86 + t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t42, 0, 0, 0, 0.2e1 * t87, 0.2e1 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t19, 0, t24, 0, 0, 0, 0, 0, t9, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, -t69, t18, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t97, 0, -t36 * t53 - t57 * t38, t36 * t97 - t60 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t97, 0, -t48 * t53, t48 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
