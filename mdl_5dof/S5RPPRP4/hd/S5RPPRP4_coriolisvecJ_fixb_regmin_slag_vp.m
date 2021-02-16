% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% tauc_reg [5x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 17:13:15
% EndTime: 2021-01-15 17:13:18
% DurationCPUTime: 0.35s
% Computational Cost: add. (489->99), mult. (944->145), div. (0->0), fcn. (465->4), ass. (0->76)
t38 = sin(pkin(7));
t44 = qJD(4) ^ 2;
t45 = qJD(1) ^ 2;
t90 = t38 * (t44 + t45);
t41 = sin(qJ(4));
t36 = t41 ^ 2;
t42 = cos(qJ(4));
t37 = t42 ^ 2;
t89 = (t36 + t37) * t45;
t39 = cos(pkin(7));
t72 = qJD(2) * t39;
t55 = -qJD(5) + t72;
t49 = t55 * t42;
t71 = t41 * qJD(4);
t64 = qJ(5) * t71;
t66 = qJD(4) * qJD(3);
t43 = -pkin(1) - pkin(2);
t24 = t43 * qJD(1) + qJD(2);
t69 = qJD(1) * qJ(2);
t17 = t38 * t24 + t39 * t69;
t13 = -qJD(1) * pkin(6) + t17;
t9 = t13 * t71;
t1 = t42 * t66 - t9 + (t49 + t64) * qJD(1);
t85 = t42 * t13;
t50 = -t41 * qJD(3) - t85;
t68 = qJD(1) * qJD(2);
t67 = qJD(1) * qJD(4);
t62 = t42 * t67;
t74 = qJD(1) * t41;
t81 = qJ(5) * t62 + qJD(5) * t74;
t2 = -t41 * t39 * t68 + t50 * qJD(4) + t81;
t35 = t42 * qJD(3);
t61 = -t41 * t13 + t35;
t70 = qJ(5) * qJD(1);
t6 = t41 * t70 + t61;
t76 = qJD(4) * pkin(4);
t3 = t6 + t76;
t7 = -t42 * t70 - t50;
t88 = -t1 * t42 + t2 * t41 + (t3 * t42 + t41 * t7) * qJD(4);
t87 = t3 - t6;
t86 = t42 * t7;
t84 = t42 * t45;
t83 = t44 * t41;
t82 = t44 * t42;
t80 = t39 * qJ(2) + t38 * t43;
t79 = t36 - t37;
t21 = -pkin(6) + t80;
t75 = qJ(5) - t21;
t73 = qJD(1) * t42;
t65 = 0.2e1 * t68;
t26 = t38 * t68;
t63 = t41 * t67;
t16 = t39 * t24 - t38 * t69;
t60 = -t38 * qJ(2) + t39 * t43;
t20 = pkin(3) - t60;
t18 = t42 * pkin(4) + t20;
t12 = qJD(1) * pkin(3) - t16;
t8 = pkin(4) * t73 + qJD(5) + t12;
t59 = -qJD(1) * t18 - t8;
t19 = -pkin(4) * t63 + t26;
t23 = -pkin(4) * t71 + t38 * qJD(2);
t58 = qJD(1) * t23 + t19;
t57 = qJD(4) * t75;
t56 = 0.2e1 * t62;
t52 = t3 * t41 - t86;
t51 = t16 * t38 - t17 * t39;
t48 = (t12 - t72) * qJD(1);
t47 = -t21 * t44 + 0.2e1 * t26;
t46 = qJD(4) * (-qJD(1) * t20 - t12 - t72);
t15 = t75 * t42;
t14 = t75 * t41;
t11 = t39 * t56 + t41 * t90;
t10 = 0.2e1 * t39 * t63 - t42 * t90;
t5 = -t55 * t41 + t42 * t57;
t4 = t41 * t57 + t49;
t22 = [0, 0, 0, 0, t65, qJ(2) * t65, ((-t38 * t60 + t39 * t80) * qJD(1) - t51) * qJD(2), t41 * t56, -0.2e1 * t79 * t67, -t82, t83, 0, t41 * t46 + t47 * t42, -t47 * t41 + t42 * t46, t58 * t42 + (t59 * t41 + t5) * qJD(4), -t58 * t41 + (t59 * t42 - t4) * qJD(4), (-t4 * t42 + t41 * t5 + (t14 * t42 - t15 * t41) * qJD(4)) * qJD(1) + t88, -t1 * t15 + t2 * t14 + t19 * t18 + t8 * t23 + t3 * t5 + t7 * t4; 0, 0, 0, 0, -t45, -t45 * qJ(2), t51 * qJD(1), 0, 0, 0, 0, 0, t10, t11, t10, t11, t39 * t89, (t52 * qJD(1) - t19) * t39 + (-qJD(1) * t8 - t88) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, -t82, -t83, -t82, 0, -t52 * qJD(4) + t1 * t41 + t2 * t42; 0, 0, 0, 0, 0, 0, 0, -t41 * t84, t79 * t45, 0, 0, 0, t41 * t48, t9 + t61 * qJD(4) + (t48 - t66) * t42, (t7 - t85) * qJD(4) + (pkin(4) * t84 - t66 + (t8 - t72) * qJD(1)) * t41 + t81, -t36 * t45 * pkin(4) + t9 + (t6 - t35) * qJD(4) + (-t64 + (-t55 + t8) * t42) * qJD(1), (t76 - t87) * t73, t87 * t7 + (t8 * t74 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t63, -0.2e1 * t62, -t89, t26 + (t86 + (-t3 - t76) * t41) * qJD(1);];
tauc_reg = t22;
