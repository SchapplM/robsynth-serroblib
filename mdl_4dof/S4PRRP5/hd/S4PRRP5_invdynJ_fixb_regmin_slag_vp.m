% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% tau_reg [4x15]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 22:36
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:36:08
% EndTime: 2021-01-14 22:36:11
% DurationCPUTime: 0.54s
% Computational Cost: add. (374->132), mult. (853->181), div. (0->0), fcn. (518->6), ass. (0->82)
t38 = qJ(4) + pkin(5);
t42 = cos(qJ(2));
t36 = sin(pkin(6));
t37 = cos(pkin(6));
t58 = g(1) * t37 + g(2) * t36;
t51 = t58 * t42;
t40 = sin(qJ(2));
t94 = g(3) * t40;
t102 = t51 + t94;
t39 = sin(qJ(3));
t61 = t40 * qJD(1) + t38 * qJD(2);
t7 = t61 * t39;
t41 = cos(qJ(3));
t8 = t61 * t41;
t44 = qJD(2) ^ 2;
t101 = qJDD(2) * t40 + t44 * t42;
t81 = qJD(3) * pkin(3);
t6 = -t7 + t81;
t57 = t39 * t6 - t41 * t8;
t100 = t57 * qJD(2);
t76 = qJDD(1) - g(3);
t99 = t58 * t40 + t76 * t42;
t98 = t7 + t6;
t34 = t39 ^ 2;
t97 = pkin(3) * t34;
t93 = g(3) * t42;
t92 = t41 * pkin(3);
t91 = t39 * t42;
t90 = t40 * t41;
t89 = t41 * t42;
t88 = t41 * t44;
t77 = t42 * qJD(1);
t67 = qJD(3) * t77;
t86 = g(3) * t91 + t41 * t67;
t35 = t41 ^ 2;
t85 = t34 - t35;
t84 = t34 + t35;
t43 = qJD(3) ^ 2;
t83 = t43 + t44;
t82 = qJD(2) * pkin(2);
t79 = qJDD(3) * pkin(3);
t32 = pkin(2) + t92;
t15 = -t32 * qJD(2) + qJD(4) - t77;
t78 = t15 * qJD(2);
t75 = qJDD(2) * t32;
t73 = qJDD(3) * t39;
t33 = t39 * qJDD(2);
t72 = t41 * qJDD(2);
t71 = t42 * qJDD(1);
t70 = t42 * qJDD(2);
t69 = qJD(1) * qJD(2);
t68 = qJD(2) * qJD(3);
t66 = t39 * t68;
t31 = t40 * t69;
t65 = qJD(3) * t38;
t64 = t41 * t31 + t39 * t67 + t58 * t90;
t17 = qJDD(2) * pkin(5) + t40 * qJDD(1) + t42 * t69;
t26 = -t77 - t82;
t63 = -qJD(2) * t26 - t17;
t46 = pkin(3) * t66 + qJDD(4) + t31 - t75;
t5 = t46 - t71;
t62 = t5 - t75;
t60 = 0.2e1 * t41 * t68;
t59 = -qJ(4) * qJDD(2) - t17;
t55 = -g(1) * (t36 * t41 - t37 * t91) - g(2) * (-t36 * t91 - t37 * t41) + t39 * t94;
t54 = -g(1) * (-t36 * t39 - t37 * t89) - g(2) * (-t36 * t89 + t37 * t39) + g(3) * t90;
t53 = qJD(3) * t61;
t52 = -0.2e1 * qJDD(2) * pkin(2) + pkin(5) * t43 + t31 - t71;
t50 = 0.2e1 * t66 - t72;
t49 = qJD(2) * qJD(4) - t59;
t48 = (-qJD(4) - t15) * qJD(2) + t59;
t47 = -pkin(5) * qJDD(3) + (t26 - t82) * qJD(3);
t45 = (-t58 - t69) * t40;
t21 = t38 * t41;
t20 = t38 * t39;
t10 = -t39 * qJD(4) - t41 * t65;
t9 = t41 * qJD(4) - t39 * t65;
t4 = -t50 * t42 + (-t83 * t41 - t73) * t40;
t3 = (-qJDD(3) * t40 - 0.2e1 * t42 * t68) * t41 + (t83 * t40 - t70) * t39;
t2 = -t39 * t53 + t49 * t41;
t1 = -t49 * t39 - t41 * t53 + t79;
t11 = [t76, 0, -t44 * t40 + t70, -t101, 0, 0, 0, 0, 0, t4, t3, t4, t3, t101 * t84, -g(3) + (-t5 - t100) * t42 + (t78 - t1 * t39 + t2 * t41 + (-t39 * t8 - t41 * t6) * qJD(3)) * t40; 0, qJDD(2), t99, -t76 * t40 + t51, t34 * qJDD(2) + t39 * t60, 0.2e1 * t39 * t72 - 0.2e1 * t85 * t68, t43 * t41 + t73, qJDD(3) * t41 - t43 * t39, 0, t47 * t39 + (-t52 - t93) * t41 + t64, t47 * t41 + (t45 + t52) * t39 + t86, -t20 * qJDD(3) + (-t62 - t93) * t41 + (t10 + (t15 + (-t32 - t92) * qJD(2)) * t39) * qJD(3) + t64, -t21 * qJDD(3) + (t15 * t41 - t9 + (-t32 * t41 + t97) * qJD(2)) * qJD(3) + (t45 + t62) * t39 + t86, (-qJD(3) * t6 + qJDD(2) * t21 + t2) * t41 + (-t8 * qJD(3) + qJDD(2) * t20 - t1) * t39 + (-t10 * t39 + t41 * t9 + (t20 * t41 - t21 * t39) * qJD(3) - t84 * t77) * qJD(2) - t102, t2 * t21 + t8 * t9 - t1 * t20 + t6 * t10 - t5 * t32 + t15 * t39 * t81 - g(3) * (t42 * t32 + t40 * t38) + (-t15 * t40 + t57 * t42) * qJD(1) + t58 * (t32 * t40 - t38 * t42); 0, 0, 0, 0, -t39 * t88, t85 * t44, t33, t72, qJDD(3), t63 * t39 + t55, t63 * t41 + t54, 0.2e1 * t79 + (pkin(3) * t88 + t48) * t39 + t55, t48 * t41 - t44 * t97 + t54, -pkin(3) * t33 + (-t81 + t98) * t41 * qJD(2), t98 * t8 + (t1 + (-g(1) * t36 + g(2) * t37) * t41 + (-t78 + t102) * t39) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t33 + t60, -t84 * t44, t46 - t99 + t100;];
tau_reg = t11;
