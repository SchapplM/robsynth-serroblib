% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPR1
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
% MMD_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:38:51
% EndTime: 2019-12-05 18:38:53
% DurationCPUTime: 0.53s
% Computational Cost: add. (1519->100), mult. (3493->192), div. (0->0), fcn. (3265->8), ass. (0->77)
t90 = qJD(2) + qJD(3);
t66 = sin(qJ(2));
t89 = pkin(6) + pkin(7);
t77 = qJD(2) * t89;
t51 = t66 * t77;
t69 = cos(qJ(2));
t52 = t69 * t77;
t65 = sin(qJ(3));
t68 = cos(qJ(3));
t54 = t89 * t66;
t55 = t89 * t69;
t72 = -t68 * t54 - t65 * t55;
t26 = -t72 * qJD(3) + t68 * t51 + t65 * t52;
t62 = sin(pkin(9));
t88 = pkin(3) * t62;
t87 = t62 * t65;
t63 = cos(pkin(9));
t86 = t63 * t65;
t50 = t65 * t69 + t68 * t66;
t38 = t90 * t50;
t49 = t65 * t66 - t68 * t69;
t15 = -t38 * qJ(4) - t49 * qJD(4) - t26;
t71 = t65 * t54 - t68 * t55;
t27 = t71 * qJD(3) + t65 * t51 - t68 * t52;
t37 = t90 * t49;
t16 = t37 * qJ(4) - t50 * qJD(4) + t27;
t8 = t63 * t15 + t62 * t16;
t32 = -t50 * qJ(4) + t72;
t33 = -t49 * qJ(4) - t71;
t18 = t62 * t32 + t63 * t33;
t59 = t68 * pkin(2) + pkin(3);
t46 = -pkin(2) * t87 + t63 * t59;
t43 = pkin(4) + t46;
t57 = t63 * pkin(3) + pkin(4);
t85 = -t43 - t57;
t84 = pkin(2) * qJD(3);
t64 = sin(qJ(5));
t83 = qJD(5) * t64;
t82 = t66 * qJD(2);
t81 = t69 * qJD(2);
t80 = -0.2e1 * pkin(1) * qJD(2);
t61 = pkin(2) * t82;
t79 = t65 * t84;
t78 = t68 * t84;
t60 = -t69 * pkin(2) - pkin(1);
t34 = t38 * pkin(3) + t61;
t7 = -t62 * t15 + t63 * t16;
t17 = t63 * t32 - t62 * t33;
t44 = (-t62 * t68 - t86) * t84;
t47 = pkin(2) * t86 + t62 * t59;
t76 = -t64 * t44 + t47 * t83;
t45 = (t63 * t68 - t87) * t84;
t67 = cos(qJ(5));
t75 = t67 * t44 - t64 * t45;
t35 = -t63 * t49 - t62 * t50;
t36 = -t62 * t49 + t63 * t50;
t74 = t67 * t35 - t64 * t36;
t22 = t64 * t35 + t67 * t36;
t70 = t49 * pkin(3) + t60;
t56 = t83 * t88;
t42 = (-t57 * t64 - t67 * t88) * qJD(5);
t41 = -qJD(5) * t67 * t57 + t56;
t28 = -t35 * pkin(4) + t70;
t25 = -t63 * t37 - t62 * t38;
t24 = t62 * t37 - t63 * t38;
t20 = (-t43 * t64 - t47 * t67) * qJD(5) + t75;
t19 = (-qJD(5) * t43 - t45) * t67 + t76;
t14 = -t24 * pkin(4) + t34;
t10 = t35 * pkin(8) + t18;
t9 = -t36 * pkin(8) + t17;
t6 = t22 * qJD(5) - t67 * t24 + t64 * t25;
t5 = t74 * qJD(5) + t64 * t24 + t67 * t25;
t4 = t24 * pkin(8) + t8;
t3 = -t25 * pkin(8) + t7;
t2 = t67 * t3 - t64 * t4 + (-t10 * t67 - t64 * t9) * qJD(5);
t1 = -t64 * t3 - t67 * t4 + (t10 * t64 - t67 * t9) * qJD(5);
t11 = [0, 0, 0, 0.2e1 * t66 * t81, 0.2e1 * (-t66 ^ 2 + t69 ^ 2) * qJD(2), 0, 0, 0, t66 * t80, t69 * t80, -0.2e1 * t50 * t37, 0.2e1 * t37 * t49 - 0.2e1 * t50 * t38, 0, 0, 0, 0.2e1 * t60 * t38 + 0.2e1 * t49 * t61, -0.2e1 * t60 * t37 + 0.2e1 * t50 * t61, -0.2e1 * t17 * t25 + 0.2e1 * t18 * t24 + 0.2e1 * t8 * t35 - 0.2e1 * t7 * t36, 0.2e1 * t17 * t7 + 0.2e1 * t18 * t8 + 0.2e1 * t70 * t34, 0.2e1 * t22 * t5, -0.2e1 * t22 * t6 + 0.2e1 * t5 * t74, 0, 0, 0, -0.2e1 * t14 * t74 + 0.2e1 * t28 * t6, 0.2e1 * t14 * t22 + 0.2e1 * t28 * t5; 0, 0, 0, 0, 0, t81, -t82, 0, -pkin(6) * t81, pkin(6) * t82, 0, 0, -t37, -t38, 0, t27, t26, t47 * t24 - t46 * t25 + t45 * t35 - t44 * t36, t17 * t44 + t18 * t45 + t7 * t46 + t8 * t47, 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t79, -0.2e1 * t78, 0, 0.2e1 * t46 * t44 + 0.2e1 * t47 * t45, 0, 0, 0, 0, 0, 0.2e1 * t20, 0.2e1 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t38, 0, t27, t26, (t24 * t62 - t25 * t63) * pkin(3), (t62 * t8 + t63 * t7) * pkin(3), 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t78, 0, (t44 * t63 + t45 * t62) * pkin(3), 0, 0, 0, 0, 0, ((-t47 - t88) * t67 + t85 * t64) * qJD(5) + t75, t56 + (t85 * qJD(5) - t45) * t67 + t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t42, 0.2e1 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
