% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PPPRR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% tau_reg [5x13]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPPRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:39
% EndTime: 2019-12-05 14:59:40
% DurationCPUTime: 0.37s
% Computational Cost: add. (273->87), mult. (656->158), div. (0->0), fcn. (613->10), ass. (0->56)
t42 = qJD(4) ^ 2;
t34 = cos(pkin(9));
t31 = sin(pkin(9));
t36 = cos(pkin(7));
t75 = t36 * t31;
t33 = sin(pkin(7));
t35 = cos(pkin(8));
t76 = t33 * t35;
t16 = t34 * t76 - t75;
t74 = t36 * t34;
t18 = t33 * t31 + t35 * t74;
t40 = cos(qJ(4));
t32 = sin(pkin(8));
t38 = sin(qJ(4));
t78 = t32 * t38;
t19 = t34 * t78 + t35 * t40;
t65 = qJDD(1) * t32;
t22 = t31 * qJDD(2) + t34 * t65;
t25 = -t35 * qJDD(1) + qJDD(3);
t77 = t32 * t40;
t82 = -g(1) * (-t18 * t38 + t36 * t77) - g(2) * (-t16 * t38 + t33 * t77) + g(3) * t19 - t38 * t22 + t40 * t25;
t81 = -t40 * t22 - t38 * t25 + g(1) * (t18 * t40 + t36 * t78) + g(2) * (t16 * t40 + t33 * t78);
t80 = t31 * t32;
t79 = t31 * t40;
t37 = sin(qJ(5));
t29 = t37 ^ 2;
t39 = cos(qJ(5));
t73 = -t39 ^ 2 + t29;
t41 = qJD(5) ^ 2;
t72 = t41 + t42;
t71 = qJD(4) * pkin(4);
t20 = t34 * t77 - t35 * t38;
t70 = qJD(5) * t20;
t67 = t19 * qJD(4);
t66 = qJDD(1) - g(3);
t64 = qJDD(5) * t37;
t63 = qJDD(5) * t40;
t62 = t34 * qJDD(5);
t61 = t39 * qJDD(4);
t60 = t40 * qJDD(4);
t59 = qJD(4) * qJD(5);
t58 = -g(1) * t33 + g(2) * t36;
t57 = t37 * t59;
t56 = 0.2e1 * qJD(4) * t31 * t38;
t51 = -t42 * t38 + t60;
t50 = qJDD(4) * t38 + t42 * t40;
t21 = -t34 * qJDD(2) + t31 * t65;
t49 = g(1) * (-t33 * t34 + t35 * t75) + g(2) * (t31 * t76 + t74) - t21;
t48 = t19 * qJDD(4) + t42 * t20;
t46 = -qJD(5) * t80 + 0.2e1 * t67;
t45 = -qJDD(4) * pkin(6) + qJD(4) * t71 + t81;
t44 = -pkin(6) * qJDD(5) - 0.2e1 * t71 * qJD(5);
t43 = 0.2e1 * qJDD(4) * pkin(4) - pkin(6) * t41 + t82;
t10 = t20 * t39 + t37 * t80;
t9 = -t20 * t37 + t39 * t80;
t1 = [t66, -g(3) + (t32 ^ 2 + t35 ^ 2) * qJDD(1), -t25 * t35 - g(3) + (t21 * t31 + t22 * t34) * t32, 0, -t48, qJD(4) * t67 - t20 * qJDD(4), 0, 0, 0, 0, 0, t9 * qJDD(5) - t48 * t39 + (t46 * t37 - t39 * t70) * qJD(5), -t10 * qJDD(5) + t48 * t37 + (t37 * t70 + t46 * t39) * qJD(5); 0, qJDD(2) + t58, -t21 * t34 + t22 * t31 + t58, 0, -t50 * t31, -t51 * t31, 0, 0, 0, 0, 0, -t39 * t62 + (-t37 * t63 - t50 * t39) * t31 + (t37 * t56 + (t37 * t34 - t39 * t79) * qJD(5)) * qJD(5), t37 * t62 + (t50 * t37 - t39 * t63) * t31 + (t39 * t56 + (t39 * t34 + t37 * t79) * qJD(5)) * qJD(5); 0, 0, qJDD(3) - t66 * t35 + (-g(1) * t36 - g(2) * t33) * t32, 0, t51, -t50, 0, 0, 0, 0, 0, (-0.2e1 * t57 + t61) * t40 + (-t72 * t39 - t64) * t38, (-qJDD(5) * t38 - 0.2e1 * t40 * t59) * t39 + (t72 * t38 - t60) * t37; 0, 0, 0, qJDD(4), t82, g(3) * t20 + t81, t29 * qJDD(4) + 0.2e1 * t39 * t57, 0.2e1 * t37 * t61 - 0.2e1 * t73 * t59, t41 * t39 + t64, qJDD(5) * t39 - t41 * t37, 0, t44 * t37 + t43 * t39, -t43 * t37 + t44 * t39; 0, 0, 0, 0, 0, 0, -t37 * t42 * t39, t73 * t42, t37 * qJDD(4), t61, qJDD(5), -g(3) * t9 + t45 * t37 - t49 * t39, g(3) * t10 + t49 * t37 + t45 * t39;];
tau_reg = t1;
