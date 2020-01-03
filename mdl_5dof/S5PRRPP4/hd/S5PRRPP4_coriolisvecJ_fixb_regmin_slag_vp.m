% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% 
% Output:
% tauc_reg [5x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:41:12
% EndTime: 2019-12-31 17:41:14
% DurationCPUTime: 0.42s
% Computational Cost: add. (377->117), mult. (901->160), div. (0->0), fcn. (407->2), ass. (0->81)
t88 = pkin(3) + pkin(4);
t63 = t88 * qJD(3);
t44 = sin(qJ(3));
t75 = qJD(2) * t44;
t34 = pkin(6) * t75;
t45 = cos(qJ(3));
t22 = t45 * qJD(1) - t34;
t69 = qJ(5) * qJD(2);
t11 = t44 * t69 + t22;
t70 = qJD(4) - t11;
t4 = -t63 + t70;
t90 = qJD(4) - t22;
t60 = t44 * qJ(4) + pkin(2);
t18 = t88 * t45 + t60;
t76 = qJD(2) * t18;
t6 = qJD(5) + t76;
t77 = qJD(5) + t6;
t89 = t77 * t44;
t23 = t45 * qJD(2) * pkin(6) + t44 * qJD(1);
t13 = -t45 * t69 + t23;
t41 = qJD(3) * qJ(4);
t8 = t13 + t41;
t19 = t41 + t23;
t16 = -qJD(3) * pkin(3) + t90;
t87 = t45 * t8;
t48 = qJD(2) ^ 2;
t86 = t45 * t48;
t47 = qJD(3) ^ 2;
t85 = t47 * t44;
t38 = t47 * t45;
t84 = pkin(6) - qJ(5);
t66 = qJD(2) * qJD(3);
t61 = t45 * t66;
t73 = t44 * qJD(4);
t83 = qJ(4) * t61 + qJD(2) * t73;
t67 = qJD(1) * qJD(3);
t20 = pkin(6) * t61 + t44 * t67;
t33 = t45 * t67;
t40 = qJD(3) * qJD(4);
t82 = t33 + 0.2e1 * t40;
t81 = t33 + t40;
t42 = t44 ^ 2;
t43 = t45 ^ 2;
t80 = t42 - t43;
t79 = pkin(6) * qJD(3);
t78 = qJ(4) * t45;
t24 = -t45 * pkin(3) - t60;
t17 = qJD(2) * t24;
t74 = qJD(3) * t44;
t72 = t44 * qJD(5);
t71 = t45 * qJD(5);
t68 = qJ(5) * qJD(3);
t65 = t44 * t86;
t64 = t45 * t68;
t26 = t84 * t45;
t62 = t44 * t66;
t2 = -t88 * t62 + t83;
t50 = -t88 * t44 + t78;
t5 = t50 * qJD(3) + t73;
t59 = qJD(2) * t5 + t2;
t58 = t6 + t76;
t57 = 0.2e1 * t17;
t55 = -0.2e1 * t62;
t54 = 0.2e1 * t61;
t53 = t23 * qJD(3) - t20;
t52 = pkin(3) * t44 - t78;
t15 = t52 * qJD(3) - t73;
t7 = pkin(3) * t62 - t83;
t51 = -pkin(6) * t47 - qJD(2) * t15 - t7;
t10 = -pkin(6) * t62 + t81;
t49 = t10 * t45 + t20 * t44 + (t16 * t45 - t19 * t44) * qJD(3);
t28 = qJ(5) * t62;
t27 = -t42 * t48 - t47;
t25 = t84 * t44;
t21 = t52 * qJD(2);
t14 = qJD(3) * t26 - t72;
t12 = -t84 * t74 - t71;
t9 = t50 * qJD(2);
t3 = (-t64 - t72) * qJD(2) + t20;
t1 = t28 + (-pkin(6) * t74 - t71) * qJD(2) + t81;
t29 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, -t38, -t85, 0, t38, t10 * t44 - t20 * t45 + (t16 * t44 + t19 * t45) * qJD(3), -t85, t38, 0, t1 * t44 - t3 * t45 + (t4 * t44 + t87) * qJD(3); 0, 0, 0, 0, t44 * t54, -0.2e1 * t80 * t66, t38, -t85, 0, pkin(2) * t55 - pkin(6) * t38, -0.2e1 * pkin(2) * t61 + pkin(6) * t85, t51 * t45 + t57 * t74, t49, -t57 * t45 * qJD(3) + t51 * t44, t49 * pkin(6) + t17 * t15 + t7 * t24, t59 * t45 + (-t58 * t44 - t14) * qJD(3), t59 * t44 + (t58 * t45 + t12) * qJD(3), -t1 * t45 - t3 * t44 + (-t4 * t45 + t44 * t8) * qJD(3) + (-t12 * t45 - t14 * t44 + (-t25 * t45 + t26 * t44) * qJD(3)) * qJD(2), t1 * t26 + t8 * t12 + t4 * t14 + t2 * t18 + t3 * t25 + t6 * t5; 0, 0, 0, 0, -t65, t80 * t48, 0, 0, 0, t48 * pkin(2) * t44 + t53, pkin(2) * t86 - t33 + (t22 + t34) * qJD(3), (-t17 * t44 + t21 * t45) * qJD(2) + t53, 0, -t22 * qJD(3) + (t17 * t45 + (t21 - t79) * t44) * qJD(2) + t82, -t20 * pkin(3) + t10 * qJ(4) - t16 * t23 - t17 * t21 + t90 * t19, t13 * qJD(3) + ((-t9 + t68) * t45 + t89) * qJD(2) - t20, -t11 * qJD(3) + t28 + (-t77 * t45 + (-t9 - t79) * t44) * qJD(2) + t82, 0, t1 * qJ(4) - t4 * t13 - t3 * t88 - t6 * t9 + t70 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, 0, t27, -t19 * qJD(3) + t17 * t75 + t20, -t65, t27, 0, -t8 * qJD(3) + (-t64 - t89) * qJD(2) + t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t54, (-t42 - t43) * t48, (t87 + (t4 - t63) * t44) * qJD(2) + t83;];
tauc_reg = t29;
