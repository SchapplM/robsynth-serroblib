% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:43:04
% EndTime: 2019-12-05 15:43:07
% DurationCPUTime: 0.50s
% Computational Cost: add. (695->113), mult. (1861->170), div. (0->0), fcn. (1468->6), ass. (0->80)
t61 = cos(pkin(9));
t65 = cos(qJ(4));
t92 = t65 * t61;
t60 = sin(pkin(9));
t63 = sin(qJ(4));
t94 = t63 * t60;
t73 = -t92 + t94;
t39 = t73 * qJD(2);
t64 = cos(qJ(5));
t31 = t64 * t39;
t86 = t65 * qJD(4);
t89 = qJD(2) * t61;
t51 = t86 * t89;
t84 = qJD(2) * t94;
t34 = -qJD(4) * t84 + t51;
t45 = t65 * t60 + t63 * t61;
t42 = t45 * qJD(4);
t35 = qJD(2) * t42;
t40 = t45 * qJD(2);
t62 = sin(qJ(5));
t88 = qJD(5) * t62;
t2 = -qJD(5) * t31 + t64 * t34 - t62 * t35 - t40 * t88;
t18 = t62 * t40 + t31;
t59 = qJD(4) + qJD(5);
t96 = t18 * t59;
t109 = t2 + t96;
t77 = -t62 * t39 + t64 * t40;
t108 = t77 * t18;
t3 = t77 * qJD(5) + t62 * t34 + t64 * t35;
t97 = t77 * t59;
t107 = -t3 + t97;
t106 = -t18 ^ 2 + t77 ^ 2;
t91 = pkin(6) + qJ(3);
t49 = t91 * t60;
t56 = t61 * qJD(1);
t37 = -qJD(2) * t49 + t56;
t85 = qJ(3) * qJD(2);
t47 = t60 * qJD(1) + t61 * t85;
t38 = pkin(6) * t89 + t47;
t78 = -t63 * t37 - t65 * t38;
t14 = -t39 * pkin(7) - t78;
t54 = -t61 * pkin(3) - pkin(2);
t48 = t54 * qJD(2) + qJD(3);
t24 = t39 * pkin(4) + t48;
t71 = t45 * qJD(3);
t70 = qJD(2) * t71;
t5 = -t34 * pkin(7) + t78 * qJD(4) - t70;
t105 = t24 * t18 + t14 * t88 + (-t14 * t59 - t5) * t62;
t103 = t65 * t37 - t63 * t38;
t102 = qJD(3) * t39;
t101 = qJD(5) - t59;
t4 = -t35 * pkin(7) + t103 * qJD(4) - t102;
t100 = -t24 * t77 - t62 * t4 + t64 * t5;
t99 = pkin(4) * t40;
t41 = t73 * qJD(4);
t76 = -t62 * t45 - t64 * t73;
t7 = t76 * qJD(5) - t64 * t41 - t62 * t42;
t98 = t7 * t59;
t93 = t64 * t14;
t90 = t60 ^ 2 + t61 ^ 2;
t87 = t41 * qJD(4);
t13 = -t40 * pkin(7) + t103;
t10 = qJD(4) * pkin(4) + t13;
t82 = -pkin(4) * t59 - t10;
t80 = t90 * qJD(2);
t23 = t64 * t45 - t62 * t73;
t75 = (-t60 * t85 + t56) * t60 - t47 * t61;
t50 = t91 * t61;
t74 = t63 * t49 - t65 * t50;
t69 = -t49 * t86 + qJD(3) * t92 + (-qJD(3) * t60 - qJD(4) * t50) * t63;
t67 = t74 * qJD(4) - t71;
t36 = t42 * qJD(4);
t29 = pkin(4) * t73 + t54;
t16 = -pkin(7) * t73 - t74;
t15 = -t45 * pkin(7) - t65 * t49 - t63 * t50;
t12 = t41 * pkin(7) + t67;
t11 = -t42 * pkin(7) + t69;
t8 = t23 * qJD(5) - t62 * t41 + t64 * t42;
t6 = t8 * t59;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t87, 0, 0, 0, 0, 0, -t6, -t98; 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t80, (qJ(3) * t80 - t75) * qJD(3), t34 * t45 - t40 * t41, -t34 * t73 - t45 * t35 + t41 * t39 - t40 * t42, -t87, -t36, 0, t67 * qJD(4) + t54 * t35 + t48 * t42, -t69 * qJD(4) + t54 * t34 - t48 * t41, t2 * t23 + t7 * t77, -t7 * t18 + t2 * t76 - t23 * t3 - t77 * t8, t98, -t6, 0, t29 * t3 + t24 * t8 + (-t62 * t11 + t64 * t12 + (-t15 * t62 - t16 * t64) * qJD(5)) * t59 + (t42 * t18 - t35 * t76) * pkin(4), t29 * t2 + t24 * t7 - (t64 * t11 + t62 * t12 + (t15 * t64 - t16 * t62) * qJD(5)) * t59 + (t35 * t23 + t42 * t77) * pkin(4); 0, 0, 0, 0, 0, 0, -t90 * qJD(2) ^ 2, t75 * qJD(2), 0, 0, 0, 0, 0, 0.2e1 * t40 * qJD(4), t51 + (-t39 - t84) * qJD(4), 0, 0, 0, 0, 0, t3 + t97, t2 - t96; 0, 0, 0, 0, 0, 0, 0, 0, t40 * t39, -t39 ^ 2 + t40 ^ 2, t51 + (t39 - t84) * qJD(4), 0, 0, -t48 * t40 - t70, t48 * t39 + t102, t108, t106, t109, t107, 0, -t18 * t99 - (-t62 * t13 - t93) * t59 + (t82 * t62 - t93) * qJD(5) + t100, -t77 * t99 + (t82 * qJD(5) + t13 * t59 - t4) * t64 + t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, t106, t109, t107, 0, t101 * (-t62 * t10 - t93) + t100, (-t101 * t10 - t4) * t64 + t105;];
tauc_reg = t1;
