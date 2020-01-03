% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% tau_reg [5x16]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:51:04
% EndTime: 2019-12-31 17:51:05
% DurationCPUTime: 0.42s
% Computational Cost: add. (481->118), mult. (761->141), div. (0->0), fcn. (380->8), ass. (0->73)
t42 = sin(qJ(4));
t44 = cos(qJ(4));
t40 = cos(pkin(7));
t26 = -pkin(1) * t40 - pkin(2);
t19 = -pkin(6) + t26;
t14 = qJD(1) * t19 + qJD(3);
t73 = qJ(5) * qJD(1);
t89 = t14 - t73;
t5 = t44 * qJD(2) + t42 * t89;
t92 = qJD(4) * t5;
t39 = sin(pkin(7));
t23 = t39 * pkin(1) + qJ(3);
t81 = t23 * qJDD(1);
t33 = qJ(1) + pkin(7);
t30 = sin(t33);
t31 = cos(t33);
t91 = -g(1) * t30 + g(2) * t31;
t76 = qJ(5) - t19;
t11 = t76 * t42;
t90 = qJD(4) * t11;
t18 = qJD(1) * t23;
t51 = -t18 * qJD(1) + t91;
t82 = t42 * pkin(4);
t55 = t23 + t82;
t10 = qJD(1) * t55 + qJD(5);
t88 = -t10 * qJD(1) + t91;
t87 = t26 * qJDD(1);
t86 = 0.2e1 * qJD(4) * t18 + qJDD(4) * t19;
t9 = t44 * t14;
t4 = -qJD(2) * t42 - t44 * t73 + t9;
t77 = qJD(4) * pkin(4);
t3 = t4 + t77;
t85 = t3 - t4;
t36 = t42 ^ 2;
t37 = t44 ^ 2;
t80 = -t36 - t37;
t79 = t36 - t37;
t46 = qJD(4) ^ 2;
t47 = qJD(1) ^ 2;
t78 = -t46 - t47;
t12 = t76 * t44;
t75 = qJD(4) * t12;
t38 = qJDD(2) - g(3);
t71 = qJDD(4) * t42;
t70 = qJDD(4) * t44;
t69 = t42 * qJDD(1);
t68 = t44 * qJDD(1);
t67 = qJD(1) * qJD(4);
t35 = qJD(3) * qJD(1);
t45 = cos(qJ(1));
t66 = pkin(1) * t45 + pkin(2) * t31 + qJ(3) * t30;
t15 = t35 + t81;
t43 = sin(qJ(1));
t65 = -t43 * pkin(1) + qJ(3) * t31;
t63 = t44 * t67;
t60 = -g(1) * t31 - g(2) * t30;
t59 = g(1) * t43 - g(2) * t45;
t58 = t3 * t44 + t42 * t5;
t54 = qJDD(5) + t15 + (t63 + t69) * pkin(4);
t53 = -qJ(5) * qJDD(1) - qJD(1) * qJD(5);
t52 = qJDD(3) + t87;
t50 = t60 + t81;
t49 = -t19 * t46 + t15 + t35 + t50;
t41 = -qJ(5) - pkin(6);
t17 = -t42 * t46 + t70;
t16 = -t44 * t46 - t71;
t13 = qJDD(1) * t19 + qJDD(3);
t8 = t44 * t13;
t7 = -t42 * qJD(5) - t75;
t6 = -t44 * qJD(5) + t90;
t2 = (qJD(4) * t89 + qJDD(2)) * t44 + (-qJD(2) * qJD(4) + t13 + t53) * t42;
t1 = qJDD(4) * pkin(4) - t42 * qJDD(2) + t53 * t44 + t8 - t92;
t20 = [qJDD(1), t59, g(1) * t45 + g(2) * t43, (t59 + (t39 ^ 2 + t40 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(3) + t91 + 0.2e1 * t87, 0.2e1 * t35 + t50 + t81, t15 * t23 + t18 * qJD(3) + t52 * t26 - g(1) * (-pkin(2) * t30 + t65) - g(2) * t66, qJDD(1) * t37 - 0.2e1 * t42 * t63, -0.2e1 * t42 * t68 + 0.2e1 * t67 * t79, t17, t16, 0, t49 * t42 + t44 * t86, -t42 * t86 + t49 * t44, (-t92 + qJDD(1) * t12 - t1 + (-t6 + t90) * qJD(1)) * t44 + (qJD(4) * t3 + qJDD(1) * t11 - t2 + (-t7 - t75) * qJD(1)) * t42 - t91, -t2 * t11 + t5 * t7 - t1 * t12 + t3 * t6 + t54 * t55 + t10 * (t44 * t77 + qJD(3)) - g(1) * (t31 * t82 + (-pkin(2) + t41) * t30 + t65) - g(2) * (t30 * t82 - t31 * t41 + t66); 0, 0, 0, t38, 0, 0, t38, 0, 0, 0, 0, 0, t16, -t17, 0, -qJD(4) * t58 - t1 * t42 + t2 * t44 - g(3); 0, 0, 0, 0, qJDD(1), -t47, t52 + t51, 0, 0, 0, 0, 0, t42 * t78 + t70, t44 * t78 - t71, t80 * qJDD(1), t1 * t44 + t2 * t42 + (-t3 * t42 + t44 * t5) * qJD(4) + t88; 0, 0, 0, 0, 0, 0, 0, t44 * t47 * t42, -t79 * t47, t68, -t69, qJDD(4), -t38 * t42 + t44 * t51 + t8, t9 * qJD(4) + (-qJD(4) * t14 - t38) * t44 + (-t13 - t51) * t42, -pkin(4) * t68 + (t77 - t85) * t42 * qJD(1), t85 * t5 + (g(3) * t42 + t44 * t88 + t1) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80 * t47, qJD(1) * t58 + t54 + t60;];
tau_reg = t20;
