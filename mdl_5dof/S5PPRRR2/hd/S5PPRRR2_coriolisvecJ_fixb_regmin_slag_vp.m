% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PPRRR2
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
% tauc_reg [5x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:48
% EndTime: 2019-12-05 15:14:51
% DurationCPUTime: 0.55s
% Computational Cost: add. (433->93), mult. (1165->149), div. (0->0), fcn. (884->8), ass. (0->76)
t50 = sin(pkin(9));
t51 = cos(pkin(9));
t54 = sin(qJ(3));
t57 = cos(qJ(3));
t103 = -t54 * t50 + t57 * t51;
t22 = t103 * qJD(1);
t56 = cos(qJ(4));
t45 = -t56 * pkin(4) - pkin(3);
t17 = t45 * qJD(3) - t22;
t104 = t17 + t22;
t30 = t57 * t50 + t54 * t51;
t25 = t30 * qJD(3);
t52 = sin(qJ(5));
t53 = sin(qJ(4));
t55 = cos(qJ(5));
t32 = t52 * t56 + t55 * t53;
t47 = qJD(4) + qJD(5);
t102 = t32 * t47;
t10 = t102 * qJD(3);
t99 = pkin(6) + pkin(7);
t100 = qJD(5) - t47;
t31 = t52 * t53 - t55 * t56;
t62 = t31 * t47;
t23 = t30 * qJD(1);
t71 = t99 * qJD(3) + t23;
t12 = t56 * qJD(2) - t71 * t53;
t21 = qJD(1) * t25;
t101 = t23 * qJD(3) - t21;
t13 = t53 * qJD(2) + t71 * t56;
t98 = t62 * t47;
t85 = qJD(3) * t56;
t75 = t55 * t85;
t86 = qJD(3) * t53;
t76 = t52 * t86;
t26 = -t75 + t76;
t28 = -t52 * t85 - t55 * t86;
t97 = t28 * t26;
t96 = t30 * t47;
t93 = t55 * t13;
t58 = qJD(4) ^ 2;
t90 = t58 * t53;
t89 = t58 * t56;
t88 = t53 ^ 2 - t56 ^ 2;
t87 = qJD(3) * pkin(3);
t84 = qJD(4) * t53;
t83 = qJD(4) * t56;
t79 = qJD(3) * qJD(4);
t78 = pkin(4) * t86;
t77 = pkin(4) * t84;
t8 = qJD(4) * pkin(4) + t12;
t74 = -pkin(4) * t47 - t8;
t73 = qJD(4) * t99;
t72 = t56 * t79;
t18 = -t22 - t87;
t24 = t103 * qJD(3);
t20 = qJD(1) * t24;
t70 = -t18 * qJD(3) - t20;
t69 = -t23 + t77;
t2 = t12 * qJD(4) + t56 * t20;
t3 = -t13 * qJD(4) - t53 * t20;
t66 = t17 * t28 - t52 * t2 + t55 * t3;
t65 = pkin(6) * t58 - t101;
t64 = qJD(4) * (t18 + t22 - t87);
t9 = qJD(5) * t75 - t47 * t76 + t55 * t72;
t61 = t17 * t26 + (t100 * t13 - t3) * t52;
t59 = qJD(3) ^ 2;
t36 = t99 * t56;
t35 = t99 * t53;
t34 = t56 * t73;
t33 = t53 * t73;
t16 = qJD(3) * t77 + t21;
t11 = t102 * t47;
t6 = -t26 ^ 2 + t28 ^ 2;
t5 = -t28 * t47 - t10;
t4 = t26 * t47 + t9;
t1 = [0, 0, 0, -t25 * qJD(3), -t24 * qJD(3), 0, 0, 0, 0, 0, -t24 * t84 - t30 * t89 + (-t103 * t84 - t25 * t56) * qJD(3), -t24 * t83 + t30 * t90 + (-t103 * t83 + t25 * t53) * qJD(3), 0, 0, 0, 0, 0, -t10 * t103 - t102 * t24 + t25 * t26 + t96 * t62, t102 * t96 - t103 * t9 + t24 * t62 - t25 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90, -t89, 0, 0, 0, 0, 0, -t11, t98; 0, 0, 0, t101, 0, 0.2e1 * t53 * t72, -0.2e1 * t88 * t79, t89, -t90, 0, t53 * t64 - t65 * t56, t65 * t53 + t56 * t64, t28 * t62 + t9 * t32, -t32 * t10 + t102 * t28 + t26 * t62 - t9 * t31, -t98, -t11, 0, (t52 * t33 - t55 * t34 + (t35 * t52 - t36 * t55) * qJD(5)) * t47 + t45 * t10 + t16 * t31 + t69 * t26 + t104 * t102, -(-t55 * t33 - t52 * t34 + (-t35 * t55 - t36 * t52) * qJD(5)) * t47 + t45 * t9 + t16 * t32 - t69 * t28 - t104 * t62; 0, 0, 0, 0, 0, -t53 * t59 * t56, t88 * t59, 0, 0, 0, t70 * t53, t70 * t56, -t97, t6, t4, t5, 0, -(-t52 * t12 - t93) * t47 - t26 * t78 + (t74 * t52 - t93) * qJD(5) + t66, t28 * t78 + (t74 * qJD(5) + t12 * t47 - t2) * t55 + t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, t6, t4, t5, 0, t66 + t100 * (-t52 * t8 - t93), (-t100 * t8 - t2) * t55 + t61;];
tauc_reg = t1;
