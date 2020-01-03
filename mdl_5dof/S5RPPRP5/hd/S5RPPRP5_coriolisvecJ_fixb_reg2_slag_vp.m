% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRP5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:43
% EndTime: 2019-12-31 17:53:45
% DurationCPUTime: 0.62s
% Computational Cost: add. (832->151), mult. (2200->195), div. (0->0), fcn. (1498->4), ass. (0->92)
t67 = sin(pkin(7));
t107 = qJD(1) * t67;
t71 = cos(qJ(4));
t93 = t71 * t107;
t68 = cos(pkin(7));
t106 = qJD(1) * t68;
t70 = sin(qJ(4));
t94 = t70 * t106;
t37 = t93 - t94;
t42 = t67 * t70 + t68 * t71;
t108 = qJD(1) * t42;
t31 = -qJD(1) * pkin(1) - pkin(2) * t106 - qJ(3) * t107 + qJD(2);
t24 = pkin(3) * t106 - t31;
t6 = pkin(4) * t108 - t37 * qJ(5) + t24;
t98 = qJD(1) * qJD(2);
t90 = t67 * t98;
t91 = t68 * t98;
t81 = t70 * t91 - t71 * t90;
t123 = t6 * t37 + t81;
t122 = -t67 * t71 + t68 * t70;
t51 = qJD(4) * t94;
t121 = (-t37 + t93) * qJD(4) - t51;
t120 = t108 ^ 2;
t119 = t37 ^ 2;
t111 = -pkin(6) + qJ(2);
t46 = t111 * t67;
t47 = t111 * t68;
t82 = t71 * t46 - t70 * t47;
t52 = qJ(2) * t107 + qJD(3);
t40 = -pkin(6) * t107 + t52;
t44 = qJD(1) * t47;
t20 = t70 * t40 + t71 * t44;
t9 = t20 * qJD(4) + t81;
t117 = t9 * t82;
t116 = t9 * t71;
t115 = t37 * t108;
t112 = t70 * t44;
t15 = qJD(4) * qJ(5) + t20;
t110 = t15 - t20;
t65 = t67 ^ 2;
t66 = t68 ^ 2;
t109 = t65 + t66;
t103 = qJD(4) * t71;
t104 = qJD(4) * t70;
t34 = t67 * t103 - t68 * t104;
t105 = qJD(4) * t34;
t11 = t42 * qJD(2) + t82 * qJD(4);
t102 = t11 * qJD(4);
t22 = t70 * t46 + t71 * t47;
t12 = t122 * qJD(2) + t22 * qJD(4);
t101 = t12 * qJD(4);
t100 = t67 * qJD(3);
t19 = t71 * t40 - t112;
t99 = qJD(5) - t19;
t97 = qJD(1) * qJD(3);
t96 = -t68 * pkin(2) - t67 * qJ(3) - pkin(1);
t95 = t40 * t103 + t70 * t90 + t71 * t91;
t92 = -t119 + t120;
t73 = qJD(1) ^ 2;
t45 = t109 * t73;
t89 = t67 * t97;
t88 = t19 + t112;
t87 = qJ(2) * t98;
t39 = t68 * pkin(3) - t96;
t86 = 0.2e1 * t108;
t85 = t66 * t87;
t29 = qJD(4) * t93 - t51;
t83 = t108 * t34 + t29 * t42;
t72 = qJD(4) ^ 2;
t80 = t37 * t107 + t72 * t71;
t33 = t42 * qJD(4);
t79 = t51 + (-t37 - t93) * qJD(4);
t28 = qJD(1) * t33;
t78 = -t103 * t108 + t37 * t104 + t71 * t28 - t70 * t29;
t77 = -t29 * pkin(4) - t28 * qJ(5) - t89;
t76 = -t108 * t11 + t12 * t37 - t122 * t9 - t22 * t29 + t28 * t82;
t75 = t108 * t33 + t122 * t29 + t28 * t42 - t37 * t34;
t53 = t65 * t87;
t32 = 0.2e1 * t109 * t98;
t30 = t33 * qJD(4);
t23 = -t107 * t108 - t72 * t70;
t18 = t37 * pkin(4) + qJ(5) * t108;
t16 = t86 * qJD(4);
t14 = -qJD(4) * pkin(4) + t99;
t13 = t42 * pkin(4) + qJ(5) * t122 + t39;
t8 = -t44 * t104 + t95;
t7 = t34 * pkin(4) + t33 * qJ(5) + qJD(5) * t122 + t100;
t5 = (qJD(5) - t112) * qJD(4) + t95;
t3 = t119 + t120;
t2 = t122 * t28 - t37 * t33;
t1 = -t37 * qJD(5) - t77;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0.2e1 * t53 + 0.2e1 * t85, 0, 0, 0, 0, 0, 0, 0.2e1 * t68 * t89, t32, 0.2e1 * t65 * t97, 0.2e1 * t85 + t53 + (t52 * qJD(2) + (-qJD(1) * t96 - t31) * qJD(3)) * t67, t2, t75, -t30, t83, -t105, 0, t86 * t100 + t24 * t34 + t39 * t29 - t101, -t102 - t24 * t33 - t39 * t28 + (-qJD(1) * t122 + t37) * t100, t19 * t33 - t20 * t34 - t8 * t42 + t76, t20 * t11 - t19 * t12 - t117 + t8 * t22 + (qJD(1) * t39 + t24) * t100, t2, -t30, -t75, 0, t105, t83, t1 * t42 + t108 * t7 + t13 * t29 + t6 * t34 - t101, -t14 * t33 - t15 * t34 - t5 * t42 + t76, t1 * t122 + t13 * t28 + t6 * t33 - t7 * t37 + t102, t1 * t13 + t15 * t11 + t14 * t12 + t5 * t22 + t6 * t7 - t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -qJ(2) * t45, 0, 0, 0, 0, 0, 0, 0, -t45, 0, -t66 * t73 * qJ(2) + (-qJD(3) - t52) * t107, 0, 0, 0, 0, 0, 0, t79, t16, t3, -t108 * t20 - t19 * t37 - t89, 0, 0, 0, 0, 0, 0, t79, t3, -t16, -t15 * t108 + (qJD(5) + t14) * t37 + t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67 * t73 * t68, 0, -t65 * t73, (qJD(2) + t31) * t107, 0, 0, 0, 0, 0, 0, t23, -t80, t78, -t24 * t107 + t8 * t70 - t116 + (-t19 * t70 + t20 * t71) * qJD(4), 0, 0, 0, 0, 0, 0, t23, t78, t80, -t6 * t107 + t5 * t70 - t116 + (t14 * t70 + t15 * t71) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, -t92, 0, -t115, -t121, 0, -t24 * t37 - t81, t88 * qJD(4) + t108 * t24 - t95, 0, 0, t115, 0, t92, 0, t121, -t115, -t108 * t18 - t123, pkin(4) * t28 - t29 * qJ(5) + t110 * t37 + (t14 - t99) * t108, t18 * t37 - t6 * t108 + (0.2e1 * qJD(5) - t88) * qJD(4) + t95, -t9 * pkin(4) + t5 * qJ(5) - t14 * t20 + t99 * t15 - t6 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, 0, -t72 - t119, -t110 * qJD(4) + t123;];
tauc_reg = t4;
