% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRRR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:51
% EndTime: 2019-12-05 15:12:53
% DurationCPUTime: 0.56s
% Computational Cost: add. (1027->99), mult. (2580->147), div. (0->0), fcn. (2108->8), ass. (0->76)
t45 = sin(pkin(9));
t46 = cos(pkin(9));
t49 = sin(qJ(3));
t52 = cos(qJ(3));
t95 = -t49 * t45 + t52 * t46;
t29 = t95 * qJD(1);
t26 = qJD(3) * pkin(3) + t29;
t48 = sin(qJ(4));
t34 = t52 * t45 + t49 * t46;
t30 = t34 * qJD(1);
t51 = cos(qJ(4));
t81 = t51 * t30;
t16 = t48 * t26 + t81;
t42 = qJD(3) + qJD(4);
t14 = t42 * pkin(7) + t16;
t47 = sin(qJ(5));
t50 = cos(qJ(5));
t10 = t47 * qJD(2) + t50 * t14;
t31 = t95 * qJD(3);
t27 = qJD(1) * t31;
t32 = t34 * qJD(3);
t28 = qJD(1) * t32;
t75 = qJD(4) * t48;
t68 = -t48 * t28 - t30 * t75;
t5 = (qJD(4) * t26 + t27) * t51 + t68;
t9 = t50 * qJD(2) - t47 * t14;
t2 = t9 * qJD(5) + t50 * t5;
t3 = -t10 * qJD(5) - t47 * t5;
t94 = t2 * t50 - t3 * t47 + (-t10 * t47 - t50 * t9) * qJD(5);
t83 = t48 * t30;
t18 = t51 * t29 - t83;
t76 = pkin(3) * qJD(4);
t93 = t51 * t76 - t18;
t17 = t48 * t29 + t81;
t38 = t48 * pkin(3) + pkin(7);
t53 = qJD(5) ^ 2;
t92 = t42 * (pkin(3) * t75 - t17) + t38 * t53;
t69 = t48 * t27 + t51 * t28;
t6 = qJD(4) * t16 + t69;
t91 = t42 * pkin(4);
t61 = -t48 * t34 + t51 * t95;
t90 = t6 * t61;
t7 = qJD(4) * t61 + t51 * t31 - t48 * t32;
t89 = t7 * t42;
t20 = t51 * t34 + t48 * t95;
t8 = qJD(4) * t20 + t48 * t31 + t51 * t32;
t88 = t8 * t42;
t15 = t51 * t26 - t83;
t13 = -t15 - t91;
t74 = qJD(5) * t50;
t87 = t13 * t74 + t6 * t47;
t86 = t15 * t42;
t85 = t16 * t42;
t79 = t53 * t47;
t43 = t47 ^ 2;
t44 = t50 ^ 2;
t78 = t43 - t44;
t77 = t43 + t44;
t41 = t42 ^ 2;
t73 = t47 * t41 * t50;
t71 = -pkin(3) * t42 - t26;
t70 = -t13 * t42 - t5;
t66 = t47 * t42 * t74;
t64 = pkin(7) * t53 - t85;
t63 = t10 * t50 - t47 * t9;
t62 = t20 * t53 + t88;
t59 = qJD(5) * (t15 - t91);
t58 = qJD(5) * (-t42 * t61 - t7);
t39 = -t51 * pkin(3) - pkin(4);
t57 = qJD(5) * (t39 * t42 - t93);
t40 = t53 * t50;
t36 = -0.2e1 * t66;
t35 = 0.2e1 * t66;
t23 = -0.2e1 * t78 * t42 * qJD(5);
t11 = t13 * qJD(5) * t47;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32 * qJD(3), -t31 * qJD(3), 0, t27 * t34 - t28 * t95 - t29 * t32 + t30 * t31, 0, 0, 0, 0, 0, 0, -t88, -t89, 0, -t15 * t8 + t16 * t7 + t5 * t20 - t90, 0, 0, 0, 0, 0, 0, t47 * t58 - t50 * t62, t47 * t62 + t50 * t58, t77 * t89, t13 * t8 + t20 * t94 + t63 * t7 - t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t40, 0, qJD(5) * t63 + t2 * t47 + t3 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t42 + (t48 * t71 - t81) * qJD(4) - t69, t18 * t42 + (qJD(4) * t71 - t27) * t51 - t68, 0, t15 * t17 - t16 * t18 + (t48 * t5 - t51 * t6 + (-t15 * t48 + t16 * t51) * qJD(4)) * pkin(3), t35, t23, t40, t36, -t79, 0, t11 + t47 * t57 + (-t6 - t92) * t50, t92 * t47 + t50 * t57 + t87, t93 * t42 * t77 + t94, -t13 * t17 + t6 * t39 - t63 * t18 + (t13 * t48 + t51 * t63) * t76 + t94 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6 + t85, -t5 + t86, 0, 0, t35, t23, t40, t36, -t79, 0, t11 + t47 * t59 + (-t6 - t64) * t50, t47 * t64 + t50 * t59 + t87, -t77 * t86 + t94, -t6 * pkin(4) + pkin(7) * t94 - t13 * t16 - t15 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, t78 * t41, 0, t73, 0, 0, t70 * t47, t70 * t50, 0, 0;];
tauc_reg = t1;
