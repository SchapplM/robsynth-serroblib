% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:18
% EndTime: 2019-12-31 16:49:20
% DurationCPUTime: 0.45s
% Computational Cost: add. (837->114), mult. (2050->171), div. (0->0), fcn. (1285->6), ass. (0->89)
t61 = sin(qJ(3));
t87 = qJD(3) * t61;
t82 = pkin(3) * t87;
t105 = 0.2e1 * t82;
t63 = cos(qJ(3));
t54 = t63 * qJD(2);
t52 = qJD(3) * t54;
t50 = sin(pkin(7)) * pkin(1) + pkin(5);
t43 = t50 * qJD(1);
t74 = pkin(6) * qJD(1) + t43;
t69 = t74 * t61;
t23 = -qJD(3) * t69 + t52;
t26 = t54 - t69;
t25 = qJD(3) * pkin(3) + t26;
t62 = cos(qJ(4));
t104 = (qJD(4) * t25 + t23) * t62;
t60 = sin(qJ(4));
t41 = t60 * t63 + t62 * t61;
t36 = qJD(1) * t41;
t55 = qJD(3) + qJD(4);
t86 = t61 * qJD(2);
t27 = t74 * t63 + t86;
t103 = pkin(6) + t50;
t98 = t60 * t61;
t71 = t55 * t98;
t95 = t62 * t63;
t20 = -t55 * t95 + t71;
t102 = t20 * t55;
t88 = qJD(1) * t63;
t80 = t62 * t88;
t89 = qJD(1) * t61;
t81 = t60 * t89;
t34 = -t80 + t81;
t101 = t36 * t34;
t51 = -cos(pkin(7)) * pkin(1) - pkin(2);
t42 = -t63 * pkin(3) + t51;
t37 = qJD(1) * t42;
t100 = t37 * t36;
t99 = t60 * t27;
t97 = t61 * t43;
t96 = t62 * t27;
t64 = qJD(3) ^ 2;
t94 = t64 * t61;
t93 = t64 * t63;
t21 = t55 * t41;
t17 = t21 * qJD(1);
t92 = -t41 * t17 + t20 * t34;
t85 = qJD(1) * qJD(3);
t78 = t63 * t85;
t91 = -qJD(4) * t80 - t62 * t78;
t90 = t61 ^ 2 - t63 ^ 2;
t44 = qJD(1) * t51;
t65 = qJD(1) ^ 2;
t84 = t61 * t65 * t63;
t83 = pkin(3) * t89;
t79 = -pkin(3) * t55 - t25;
t24 = t27 * qJD(3);
t77 = -t60 * t23 - t62 * t24;
t76 = -qJD(4) * t99 - t60 * t24;
t75 = qJD(3) * t103;
t72 = t61 * t78;
t16 = qJD(1) * t71 + t91;
t40 = -t95 + t98;
t70 = -t40 * t16 + t36 * t21;
t6 = t60 * t25 + t96;
t38 = t103 * t61;
t39 = t103 * t63;
t14 = -t62 * t38 - t60 * t39;
t15 = -t60 * t38 + t62 * t39;
t31 = t63 * t43 + t86;
t68 = 0.2e1 * qJD(3) * t44;
t67 = t37 * t34 - t76;
t2 = -qJD(4) * t6 + t77;
t28 = -t43 * t87 + t52;
t29 = t31 * qJD(3);
t30 = t54 - t97;
t66 = t28 * t63 + t29 * t61 + (-t30 * t63 - t31 * t61) * qJD(3);
t33 = t63 * t75;
t32 = t61 * t75;
t18 = t21 * t55;
t12 = -t34 ^ 2 + t36 ^ 2;
t10 = t62 * t26 - t99;
t9 = -t60 * t26 - t96;
t7 = -t91 + (t34 - t81) * t55;
t5 = t62 * t25 - t99;
t4 = -t15 * qJD(4) + t60 * t32 - t62 * t33;
t3 = t14 * qJD(4) - t62 * t32 - t60 * t33;
t1 = t76 + t104;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t72, -0.2e1 * t90 * t85, t93, -0.2e1 * t72, -t94, 0, -t50 * t93 + t61 * t68, t50 * t94 + t63 * t68, t66, t66 * t50, -t16 * t41 - t36 * t20, -t70 + t92, -t102, t17 * t40 + t34 * t21, -t18, 0, t42 * t17 + t37 * t21 + t4 * t55 + (qJD(1) * t40 + t34) * t82, t36 * t105 - t42 * t16 - t37 * t20 - t3 * t55, -t1 * t40 + t14 * t16 - t15 * t17 - t2 * t41 + t5 * t20 - t6 * t21 - t3 * t34 - t4 * t36, t1 * t15 + t37 * t105 + t2 * t14 + t6 * t3 + t5 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t93, 0, t28 * t61 - t29 * t63 + (-t30 * t61 + t31 * t63) * qJD(3), 0, 0, 0, 0, 0, 0, -t18, t102, t70 + t92, t1 * t41 - t2 * t40 - t6 * t20 - t5 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, t90 * t65, 0, t84, 0, 0, -t44 * t89, -t44 * t88 - t52 + (t30 + t97) * qJD(3), 0, 0, t101, t12, t7, -t101, 0, 0, -t34 * t83 - t100 - t9 * t55 + (t79 * t60 - t96) * qJD(4) + t77, -t36 * t83 + t10 * t55 + (t79 * qJD(4) - t23) * t62 + t67, (t6 + t9) * t36 + (t10 - t5) * t34 + (t16 * t62 - t17 * t60 + (-t34 * t62 + t36 * t60) * qJD(4)) * pkin(3), -t6 * t10 - t5 * t9 + (-t37 * t89 + t1 * t60 + t2 * t62 + (-t5 * t60 + t6 * t62) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t12, t7, -t101, 0, 0, t6 * t55 - t100 + t2, t5 * t55 - t104 + t67, 0, 0;];
tauc_reg = t8;
