% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% tauc_reg [5x16]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:09:10
% EndTime: 2019-12-05 15:09:12
% DurationCPUTime: 0.35s
% Computational Cost: add. (400->86), mult. (1045->114), div. (0->0), fcn. (725->6), ass. (0->64)
t37 = sin(pkin(8));
t38 = cos(pkin(8));
t40 = sin(qJ(3));
t42 = cos(qJ(3));
t82 = -t40 * t37 + t42 * t38;
t21 = t82 * qJD(1);
t26 = t42 * t37 + t40 * t38;
t24 = t26 * qJD(3);
t23 = t82 * qJD(3);
t83 = qJD(1) * t23 + qJD(2) * qJD(4);
t41 = cos(qJ(4));
t22 = t26 * qJD(1);
t17 = qJD(3) * pkin(6) + t22;
t39 = sin(qJ(4));
t79 = t39 * t17;
t9 = t41 * qJD(2) - t79;
t81 = qJD(5) - t9;
t6 = -qJD(4) * pkin(4) + t81;
t10 = t39 * qJD(2) + t41 * t17;
t7 = qJD(4) * qJ(5) + t10;
t43 = qJD(4) ^ 2;
t80 = pkin(6) * t43;
t78 = t39 * t41;
t73 = t43 * t39;
t34 = t43 * t41;
t72 = t83 * t41;
t61 = t22 * qJD(3);
t63 = qJD(4) * t39;
t71 = t21 * t63 + t41 * t61;
t47 = pkin(4) * t39 - qJ(5) * t41;
t20 = t47 * qJD(4) - t39 * qJD(5);
t70 = t20 - t22;
t19 = qJD(1) * t24;
t35 = t39 ^ 2;
t36 = t41 ^ 2;
t69 = t35 - t36;
t68 = t35 + t36;
t67 = qJD(3) * pkin(3);
t66 = qJD(3) * t20;
t28 = -t41 * pkin(4) - t39 * qJ(5) - pkin(3);
t65 = qJD(3) * t28;
t64 = qJD(3) * t39;
t62 = qJD(4) * t41;
t60 = t23 * qJD(3);
t44 = qJD(3) ^ 2;
t56 = t44 * t78;
t4 = t17 * t62 + t83 * t39;
t55 = 0.2e1 * qJD(3) * qJD(4);
t5 = t66 + t19;
t54 = -t5 - t80;
t53 = -t19 - t80;
t52 = t9 + t79;
t16 = -t21 - t67;
t51 = t16 - t67;
t8 = -t21 + t65;
t50 = t8 + t65;
t48 = t39 * t6 + t41 * t7;
t46 = t10 * qJD(4) - t4;
t3 = (qJD(5) - t79) * qJD(4) + t72;
t45 = t3 * t41 + t4 * t39 + (-t39 * t7 + t41 * t6) * qJD(4);
t27 = t47 * qJD(3);
t2 = -t23 * t63 - t26 * t34 + (-t24 * t41 - t63 * t82) * qJD(3);
t1 = -t23 * t62 + t26 * t73 + (t24 * t39 - t62 * t82) * qJD(3);
t11 = [0, 0, 0, -t24 * qJD(3), -t60, 0, 0, 0, 0, 0, t2, t1, t2, t68 * t60, -t1, t48 * t23 + t8 * t24 + t45 * t26 - t5 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, -t34, -t73, 0, t34, t48 * qJD(4) + t3 * t39 - t4 * t41; 0, 0, 0, t61 - t19, 0, t55 * t78, -t69 * t55, t34, -t73, 0, t53 * t41 + t51 * t63 + t71, (-t53 - t61) * t39 + (t21 + t51) * t62, t50 * t63 + (t54 - t66) * t41 + t71, -t68 * t21 * qJD(3) + t45, (-t21 - t50) * t62 + (-t70 * qJD(3) + t54) * t39, t45 * pkin(6) - t48 * t21 + t5 * t28 + t70 * t8; 0, 0, 0, 0, 0, -t56, t69 * t44, 0, 0, 0, -t16 * t64 + t46, -t16 * t41 * qJD(3) + t52 * qJD(4) - t72, (t27 * t41 - t39 * t8) * qJD(3) + t46, 0, (t27 * t39 + t41 * t8) * qJD(3) + (0.2e1 * qJD(5) - t52) * qJD(4) + t72, -t4 * pkin(4) + t3 * qJ(5) - t6 * t10 - t8 * t27 + t81 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, 0, -t35 * t44 - t43, -t7 * qJD(4) + t8 * t64 + t4;];
tauc_reg = t11;
