% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:19
% EndTime: 2019-12-05 17:06:21
% DurationCPUTime: 0.49s
% Computational Cost: add. (972->102), mult. (1980->149), div. (0->0), fcn. (1056->6), ass. (0->82)
t41 = cos(qJ(5));
t35 = qJD(2) + qJD(3);
t43 = cos(qJ(3));
t79 = pkin(2) * qJD(2);
t72 = t43 * t79;
t24 = t35 * pkin(3) + t72;
t42 = cos(qJ(4));
t66 = qJD(3) * t72;
t39 = sin(qJ(4));
t40 = sin(qJ(3));
t74 = t40 * t79;
t68 = t39 * t74;
t83 = (qJD(3) + qJD(4)) * t68;
t5 = (qJD(4) * t24 + t66) * t42 - t83;
t16 = t39 * t24 + t42 * t74;
t33 = qJD(4) + t35;
t14 = t33 * pkin(8) + t16;
t38 = sin(qJ(5));
t7 = t41 * qJD(1) - t38 * t14;
t2 = t7 * qJD(5) + t41 * t5;
t8 = t38 * qJD(1) + t41 * t14;
t3 = -t8 * qJD(5) - t38 * t5;
t97 = t2 * t41 - t3 * t38 + (-t38 * t8 - t41 * t7) * qJD(5);
t86 = t39 * t40;
t58 = t42 * t43 - t86;
t21 = t58 * t79;
t77 = qJD(4) * t42;
t96 = pkin(3) * t77 - t21;
t85 = t40 * t42;
t59 = t39 * t43 + t85;
t20 = t59 * t79;
t29 = t39 * pkin(3) + pkin(8);
t44 = qJD(5) ^ 2;
t78 = qJD(4) * t39;
t95 = t33 * (pkin(3) * t78 - t20) + t29 * t44;
t94 = t59 * qJD(3) + t40 * t77;
t93 = t33 * pkin(4);
t31 = t43 * pkin(2) + pkin(3);
t9 = t31 * t77 + (t58 * qJD(3) - t40 * t78) * pkin(2);
t92 = t9 * t33;
t15 = t42 * t24 - t68;
t13 = -t15 - t93;
t47 = t94 * pkin(2);
t71 = t24 * t78;
t6 = qJD(2) * t47 + t71;
t76 = qJD(5) * t41;
t91 = t13 * t76 + t6 * t38;
t10 = t31 * t78 + t47;
t90 = t10 * t33;
t89 = t15 * t33;
t88 = t16 * t33;
t84 = t44 * t38;
t82 = pkin(2) * t85 + t39 * t31;
t36 = t38 ^ 2;
t37 = t41 ^ 2;
t81 = t36 - t37;
t80 = t36 + t37;
t32 = t33 ^ 2;
t75 = t38 * t32 * t41;
t69 = -t13 * t33 - t5;
t67 = t38 * t33 * t76;
t64 = t38 * t7 - t41 * t8;
t63 = (-qJD(3) + t35) * t79;
t62 = pkin(2) * qJD(3) * (-qJD(2) - t35);
t61 = pkin(8) * t44 - t88;
t19 = pkin(8) + t82;
t60 = t19 * t44 + t90;
t57 = qJD(5) * (t15 - t93);
t54 = -pkin(2) * t86 + t42 * t31;
t18 = -pkin(4) - t54;
t56 = qJD(5) * (t18 * t33 - t9);
t55 = (-pkin(3) * t33 - t24) * qJD(4);
t30 = -t42 * pkin(3) - pkin(4);
t52 = qJD(5) * (t30 * t33 - t96);
t46 = t94 * t79;
t45 = -t46 - t71;
t34 = t44 * t41;
t23 = -0.2e1 * t67;
t22 = 0.2e1 * t67;
t17 = -0.2e1 * t81 * t33 * qJD(5);
t11 = t13 * qJD(5) * t38;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, -t34, 0, -qJD(5) * t64 + t2 * t38 + t3 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40 * t62, t43 * t62, 0, 0, 0, 0, 0, 0, 0, 0, t45 - t90, -t5 - t92, 0, -t15 * t10 + t16 * t9 + t5 * t82 - t54 * t6, t22, t17, t34, t23, -t84, 0, t11 + t38 * t56 + (-t6 - t60) * t41, t38 * t60 + t41 * t56 + t91, t80 * t92 + t97, t13 * t10 + t6 * t18 + t19 * t97 - t64 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40 * t63, t43 * t63, 0, 0, 0, 0, 0, 0, 0, 0, t20 * t33 + t39 * t55 - t46, t21 * t33 + (t55 - t66) * t42 + t83, 0, t15 * t20 - t16 * t21 + (t39 * t5 - t42 * t6 + (-t15 * t39 + t16 * t42) * qJD(4)) * pkin(3), t22, t17, t34, t23, -t84, 0, t11 + t38 * t52 + (-t6 - t95) * t41, t95 * t38 + t41 * t52 + t91, t96 * t33 * t80 + t97, -t13 * t20 + t6 * t30 + t64 * t21 + (t13 * t39 - t42 * t64) * qJD(4) * pkin(3) + t97 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 + t88, -t5 + t89, 0, 0, t22, t17, t34, t23, -t84, 0, t11 + t38 * t57 + (-t6 - t61) * t41, t38 * t61 + t41 * t57 + t91, -t80 * t89 + t97, -t6 * pkin(4) + pkin(8) * t97 - t13 * t16 + t15 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, t81 * t32, 0, t75, 0, 0, t69 * t38, t69 * t41, 0, 0;];
tauc_reg = t1;
