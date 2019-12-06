% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [5x16]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:01:44
% EndTime: 2019-12-05 18:01:46
% DurationCPUTime: 0.38s
% Computational Cost: add. (636->94), mult. (1286->133), div. (0->0), fcn. (688->6), ass. (0->76)
t92 = qJ(5) + pkin(7);
t39 = cos(pkin(8)) * pkin(1) + pkin(2);
t51 = sin(qJ(3));
t53 = cos(qJ(3));
t89 = pkin(1) * sin(pkin(8));
t91 = t53 * t39 - t51 * t89;
t50 = sin(qJ(4));
t52 = cos(qJ(4));
t33 = t39 * qJD(1);
t73 = qJD(1) * t89;
t21 = t51 * t33 + t53 * t73;
t45 = qJD(1) + qJD(3);
t67 = t92 * t45 + t21;
t7 = t52 * qJD(2) - t67 * t50;
t8 = t50 * qJD(2) + t67 * t52;
t78 = qJD(4) * pkin(4);
t6 = t7 + t78;
t90 = t6 - t7;
t88 = t45 * pkin(3);
t87 = t52 * pkin(4);
t86 = t21 * t45;
t55 = t51 * t39 + t53 * t89;
t23 = t55 * qJD(3);
t85 = t23 * t45;
t54 = qJD(4) ^ 2;
t83 = t54 * t50;
t20 = t53 * t33 - t51 * t73;
t14 = -t20 - t88;
t63 = qJD(3) * t73;
t76 = qJD(3) * t33;
t19 = t51 * t76 + t53 * t63;
t81 = t14 * qJD(4) * t52 + t19 * t50;
t46 = t50 ^ 2;
t47 = t52 ^ 2;
t80 = -t46 - t47;
t79 = t46 - t47;
t26 = pkin(7) + t55;
t77 = -qJ(5) - t26;
t74 = t50 * qJD(4);
t71 = t45 * t74;
t10 = pkin(4) * t71 + t19;
t72 = pkin(4) * t74;
t70 = -pkin(3) - t87;
t18 = -t51 * t63 + t53 * t76;
t65 = qJD(5) * t45 + t18;
t2 = t7 * qJD(4) + t65 * t52;
t3 = -qJD(4) * t8 - t65 * t50;
t69 = t2 * t52 - t3 * t50;
t68 = -t14 * t45 - t18;
t66 = qJD(4) * t92;
t64 = qJD(4) * t77;
t25 = -pkin(3) - t91;
t62 = -t50 * t8 - t52 * t6;
t61 = t6 * t50 - t8 * t52;
t60 = pkin(7) * t54 - t86;
t59 = t26 * t54 + t85;
t58 = qJD(4) * (t20 - t88);
t22 = t91 * qJD(3);
t56 = qJD(4) * (t25 * t45 - t22);
t44 = t45 ^ 2;
t43 = t52 * qJ(5);
t42 = t54 * t52;
t40 = t52 * qJD(5);
t35 = t52 * pkin(7) + t43;
t34 = t92 * t50;
t30 = 0.2e1 * t52 * t71;
t28 = -t50 * qJD(5) - t52 * t66;
t27 = -t50 * t66 + t40;
t24 = -0.2e1 * t79 * t45 * qJD(4);
t17 = t52 * t26 + t43;
t16 = t77 * t50;
t11 = t14 * t74;
t9 = t70 * t45 + qJD(5) - t20;
t5 = (-qJD(5) - t22) * t50 + t52 * t64;
t4 = t52 * t22 + t50 * t64 + t40;
t1 = [0, 0, 0, 0, 0, -t19 - t85, -t22 * t45 - t18, t30, t24, t42, -t83, 0, t11 + t50 * t56 + (-t19 - t59) * t52, t59 * t50 + t52 * t56 + t81, (t4 * t52 - t5 * t50) * t45 + ((-t16 * t52 - t17 * t50) * t45 + t62) * qJD(4) + t69, t2 * t17 + t8 * t4 + t3 * t16 + t6 * t5 + t10 * (t25 - t87) + t9 * (t23 + t72); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, -t42, 0, -t61 * qJD(4) + t2 * t50 + t3 * t52; 0, 0, 0, 0, 0, -t19 + t86, t20 * t45 - t18, t30, t24, t42, -t83, 0, t11 + t50 * t58 + (-t19 - t60) * t52, t60 * t50 + t52 * t58 + t81, t62 * qJD(4) + (t27 * t52 - t28 * t50 + t80 * t20 + (t34 * t52 - t35 * t50) * qJD(4)) * t45 + t69, t2 * t35 + t8 * t27 - t3 * t34 + t6 * t28 + t10 * t70 + (-t21 + t72) * t9 + t61 * t20; 0, 0, 0, 0, 0, 0, 0, -t50 * t44 * t52, t79 * t44, 0, 0, 0, t68 * t50, t68 * t52, (-t78 + t90) * t52 * t45, t90 * t8 + (-t45 * t50 * t9 + t3) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80 * t44, t61 * t45 + t10;];
tauc_reg = t1;
