% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4PRRR5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4PRRR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:46
% EndTime: 2019-12-31 16:33:47
% DurationCPUTime: 0.25s
% Computational Cost: add. (704->74), mult. (999->111), div. (0->0), fcn. (660->8), ass. (0->60)
t62 = sin(qJ(4));
t65 = cos(qJ(4));
t60 = sin(pkin(7));
t61 = cos(pkin(7));
t45 = -t61 * g(1) - t60 * g(2);
t59 = -g(3) + qJDD(1);
t64 = sin(qJ(2));
t67 = cos(qJ(2));
t31 = t67 * t45 + t64 * t59;
t69 = qJD(2) ^ 2;
t29 = -t69 * pkin(2) + t31;
t63 = sin(qJ(3));
t66 = cos(qJ(3));
t30 = -t64 * t45 + t67 * t59;
t70 = qJDD(2) * pkin(2) + t30;
t15 = t66 * t29 + t63 * t70;
t56 = qJD(2) + qJD(3);
t54 = t56 ^ 2;
t55 = qJDD(2) + qJDD(3);
t13 = -t54 * pkin(3) + t55 * pkin(6) + t15;
t72 = -t60 * g(1) + t61 * g(2);
t7 = t62 * t13 - t65 * t72;
t8 = t65 * t13 + t62 * t72;
t3 = t62 * t7 + t65 * t8;
t14 = -t63 * t29 + t66 * t70;
t12 = -t55 * pkin(3) - t54 * pkin(6) - t14;
t80 = -pkin(3) * t12 + pkin(6) * t3;
t46 = t62 * t54 * t65;
t79 = t62 * (qJDD(4) + t46);
t78 = t62 * t55;
t77 = t65 * (qJDD(4) - t46);
t76 = qJD(4) * t56;
t57 = t62 ^ 2;
t49 = t57 * t54;
t68 = qJD(4) ^ 2;
t27 = -t77 - t62 * (-t49 - t68);
t35 = 0.2e1 * t65 * t76 + t78;
t75 = -pkin(3) * t35 + pkin(6) * t27 + t62 * t12;
t58 = t65 ^ 2;
t50 = t58 * t54;
t26 = t65 * (-t50 - t68) - t79;
t48 = t65 * t55;
t36 = -0.2e1 * t62 * t76 + t48;
t74 = pkin(3) * t36 + pkin(6) * t26 - t65 * t12;
t38 = (t57 + t58) * t55;
t41 = t49 + t50;
t73 = pkin(3) * t41 + pkin(6) * t38 + t3;
t39 = -t66 * t54 - t63 * t55;
t71 = t63 * t54 - t66 * t55;
t25 = t79 + t65 * (-t49 + t68);
t24 = t62 * (t50 - t68) + t77;
t21 = t35 * t62;
t20 = t36 * t65;
t19 = t63 * t38 + t66 * t41;
t18 = t65 * t35 + t62 * t36;
t17 = t63 * t27 - t66 * t35;
t16 = t63 * t26 + t66 * t36;
t4 = t66 * t14 + t63 * t15;
t1 = -t66 * t12 + t63 * t3;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, 0, 0, 0, 0, 0, t67 * qJDD(2) - t64 * t69, -t64 * qJDD(2) - t67 * t69, 0, t67 * t30 + t64 * t31, 0, 0, 0, 0, 0, 0, t64 * t39 - t67 * t71, t67 * t39 + t64 * t71, 0, t64 * (-t63 * t14 + t66 * t15) + t67 * t4, 0, 0, 0, 0, 0, 0, t64 * (t66 * t26 - t63 * t36) + t67 * t16, t64 * (t66 * t27 + t63 * t35) + t67 * t17, t64 * (t66 * t38 - t63 * t41) + t67 * t19, t64 * (t63 * t12 + t66 * t3) + t67 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t30, -t31, 0, 0, 0, 0, 0, 0, 0, t55, -pkin(2) * t71 + t14, pkin(2) * t39 - t15, 0, pkin(2) * t4, t21, t18, t25, t20, t24, 0, pkin(2) * t16 + t74, pkin(2) * t17 + t75, pkin(2) * t19 + t73, pkin(2) * t1 + t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t14, -t15, 0, 0, t21, t18, t25, t20, t24, 0, t74, t75, t73, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t49 - t50, t78, t46, t48, qJDD(4), -t7, -t8, 0, 0;];
tauJ_reg = t2;
