% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPPR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:30
% EndTime: 2019-12-31 17:46:32
% DurationCPUTime: 0.46s
% Computational Cost: add. (886->104), mult. (1785->164), div. (0->0), fcn. (1153->6), ass. (0->79)
t56 = sin(pkin(8));
t58 = cos(pkin(8));
t61 = sin(qJ(5));
t62 = cos(qJ(5));
t36 = t62 * t56 + t61 * t58;
t32 = qJD(1) * t36;
t99 = -0.2e1 * t32;
t91 = t61 * t56;
t35 = -t62 * t58 + t91;
t59 = cos(pkin(7));
t79 = t59 * qJD(2);
t44 = -qJD(4) + t79;
t39 = qJD(1) * t44;
t98 = t35 * t39;
t57 = sin(pkin(7));
t24 = t35 * t57;
t86 = t56 ^ 2 + t58 ^ 2;
t73 = t86 * t39;
t97 = t32 ^ 2;
t77 = qJD(1) * qJD(2);
t45 = t57 * t77;
t96 = 0.2e1 * t45;
t63 = -pkin(1) - pkin(2);
t87 = t59 * qJ(2) + t57 * t63;
t37 = -qJ(4) + t87;
t95 = pkin(6) - t37;
t84 = qJD(1) * t58;
t75 = t62 * t84;
t76 = qJD(1) * t91;
t31 = -t75 + t76;
t94 = t31 * t32;
t64 = qJD(1) ^ 2;
t93 = t57 * t64;
t92 = t59 * t64;
t41 = t63 * qJD(1) + qJD(2);
t78 = qJD(1) * qJ(2);
t28 = t57 * t41 + t59 * t78;
t16 = -qJD(1) * qJ(4) + t28;
t12 = t56 * qJD(3) + t58 * t16;
t89 = qJD(5) * t24 + t32 * t59;
t23 = t36 * t57;
t85 = qJD(1) * t35;
t88 = qJD(5) * t23 - t59 * t85;
t33 = t35 * qJD(5);
t83 = t33 * qJD(5);
t34 = t36 * qJD(5);
t82 = t34 * qJD(5);
t81 = t57 * qJD(1);
t80 = t57 * qJD(2);
t74 = 0.2e1 * t77;
t27 = t59 * t41 - t57 * t78;
t72 = -t57 * qJ(2) + t59 * t63;
t71 = pkin(3) - t72;
t10 = -pkin(6) * t84 + t12;
t51 = t58 * qJD(3);
t9 = t51 + (pkin(6) * qJD(1) - t16) * t56;
t4 = t62 * t10 + t61 * t9;
t3 = -t61 * t10 + t62 * t9;
t70 = -(-t56 * t16 + t51) * t56 + t12 * t58;
t21 = t95 * t56;
t22 = t95 * t58;
t7 = t62 * t21 + t61 * t22;
t8 = t61 * t21 - t62 * t22;
t40 = qJD(5) * t75;
t25 = qJD(5) * t76 - t40;
t69 = t25 * t35 - t32 * t34;
t26 = qJD(1) * t34;
t68 = -t36 * t26 + t33 * t31;
t67 = t27 * t57 - t28 * t59;
t15 = qJD(1) * pkin(3) + qJD(4) - t27;
t66 = t36 * t39;
t30 = t31 ^ 2;
t29 = t58 * pkin(4) + t71;
t13 = pkin(4) * t84 + t15;
t6 = -t8 * qJD(5) - t36 * t44;
t5 = t7 * qJD(5) - t35 * t44;
t2 = -t4 * qJD(5) - t66;
t1 = t3 * qJD(5) - t98;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, qJ(2) * t74, 0, 0, 0, 0, 0, 0, t96, t59 * t74, 0, ((-t57 * t72 + t59 * t87) * qJD(1) - t67) * qJD(2), 0, 0, 0, 0, 0, 0, t58 * t96, -0.2e1 * t56 * t45, -0.2e1 * t73, t70 * t44 + t37 * t73 + (qJD(1) * t71 + t15) * t80, -t25 * t36 - t32 * t33, t68 + t69, t83, t26 * t35 + t31 * t34, t82, 0, t6 * qJD(5) - t13 * t34 - t29 * t26 + (-t31 - t85) * t80, -t5 * qJD(5) + t13 * t33 + t29 * t25 + t80 * t99, t1 * t35 + t2 * t36 - t7 * t25 + t8 * t26 - t3 * t33 + t5 * t31 + t6 * t32 + t4 * t34, t1 * t8 + t2 * t7 + t3 * t6 + t4 * t5 + (qJD(1) * t29 + t13) * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t64 * qJ(2), 0, 0, 0, 0, 0, 0, -t93, -t92, 0, t67 * qJD(1), 0, 0, 0, 0, 0, 0, -t58 * t93, t56 * t93, t86 * t92, t57 * t73 + (-t15 * t57 + (-t70 - t80) * t59) * qJD(1), 0, 0, 0, 0, 0, 0, t89 * qJD(5) + t59 * t26 + t31 * t81, t88 * qJD(5) - t59 * t25 + t32 * t81, t23 * t25 - t24 * t26 - t88 * t31 + t89 * t32, -t1 * t24 - t2 * t23 - t88 * t4 + t89 * t3 + (-t13 - t79) * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, t83, -t68 + t69, t1 * t36 - t2 * t35 - t3 * t34 - t4 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86 * t64, t70 * qJD(1) + t45, 0, 0, 0, 0, 0, 0, qJD(5) * t99, -t40 + (t31 + t76) * qJD(5), -t30 - t97, -t3 * t32 - t4 * t31 + t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, -t30 + t97, -t40 + (-t31 + t76) * qJD(5), -t94, 0, 0, t13 * t32 - t66, -t13 * t31 + t98, 0, 0;];
tauc_reg = t11;
