% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRPRR9
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRPRR9_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR9_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR9_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:45
% EndTime: 2019-12-31 17:39:46
% DurationCPUTime: 0.36s
% Computational Cost: add. (1124->100), mult. (1697->112), div. (0->0), fcn. (926->8), ass. (0->72)
t86 = pkin(2) + pkin(3);
t53 = (-qJD(2) + qJD(4));
t51 = t53 ^ 2;
t60 = sin(qJ(5));
t63 = cos(qJ(5));
t41 = t60 * t51 * t63;
t35 = qJDD(5) + t41;
t85 = t60 * t35;
t52 = qJDD(2) - qJDD(4);
t84 = t60 * t52;
t36 = qJDD(5) - t41;
t83 = t63 * t36;
t67 = qJD(2) ^ 2;
t54 = qJDD(2) * qJ(3);
t58 = sin(pkin(8));
t59 = cos(pkin(8));
t37 = t58 * g(1) - t59 * g(2);
t38 = -t59 * g(1) - t58 * g(2);
t62 = sin(qJ(2));
t65 = cos(qJ(2));
t82 = -t62 * t37 - t65 * t38;
t79 = (2 * qJD(3) * qJD(2)) - t82;
t76 = t54 + t79;
t16 = -t86 * t67 + t76;
t61 = sin(qJ(4));
t64 = cos(qJ(4));
t77 = t65 * t37 - t62 * t38;
t73 = qJDD(3) - t77;
t70 = -t67 * qJ(3) + t73;
t68 = -t86 * qJDD(2) + t70;
t10 = t64 * t16 + t61 * t68;
t81 = qJDD(2) * pkin(2);
t80 = 2 * qJD(5) * t53;
t57 = g(3) - qJDD(1);
t8 = -t51 * pkin(4) - t52 * pkin(7) + t10;
t5 = -t63 * t57 + t60 * t8;
t6 = t60 * t57 + t63 * t8;
t2 = t60 * t5 + t63 * t6;
t78 = t61 * t16 - t64 * t68;
t7 = t52 * pkin(4) - t51 * pkin(7) + t78;
t75 = -pkin(4) * t7 + pkin(7) * t2;
t43 = t63 * t52;
t74 = t60 * t80 + t43;
t29 = -t64 * t51 + t61 * t52;
t30 = t61 * t51 + t64 * t52;
t26 = t63 * t80 - t84;
t56 = t63 ^ 2;
t46 = t56 * t51;
t66 = qJD(5) ^ 2;
t40 = -t46 - t66;
t24 = t63 * t40 - t85;
t72 = -pkin(4) * t74 + pkin(7) * t24 - t63 * t7;
t55 = t60 ^ 2;
t45 = t55 * t51;
t39 = -t45 - t66;
t25 = -t60 * t39 - t83;
t71 = pkin(4) * t26 - pkin(7) * t25 - t60 * t7;
t28 = (-t55 - t56) * t52;
t31 = t45 + t46;
t69 = pkin(4) * t31 + pkin(7) * t28 + t2;
t23 = t85 + t63 * (-t45 + t66);
t22 = t60 * (t46 - t66) + t83;
t21 = t26 * t60;
t20 = t74 * t63;
t19 = t70 - t81;
t18 = t61 * t28 + t64 * t31;
t17 = t63 * t26 - t60 * t74;
t12 = t61 * t25 - t64 * t26;
t11 = t61 * t24 - t64 * t74;
t3 = t61 * t10 - t64 * t78;
t1 = t61 * t2 - t64 * t7;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, 0, 0, 0, 0, 0, 0, -t63 * t35 - t60 * t40, t60 * t36 - t63 * t39, 0, t63 * t5 - t60 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t77, t82, 0, 0, 0, 0, 0, qJDD(2), 0, 0, -t73 + 0.2e1 * t81, 0, 0.2e1 * t54 + t79, -pkin(2) * t19 + qJ(3) * (-t67 * pkin(2) + t76), 0, 0, 0, 0, 0, t52, qJ(3) * t29 + t86 * t30 + t78, qJ(3) * t30 - t86 * t29 + t10, 0, qJ(3) * (t64 * t10 + t61 * t78) - t86 * t3, -t21, -t17, -t23, t20, -t22, 0, qJ(3) * (t64 * t24 + t61 * t74) - t86 * t11 - t72, qJ(3) * (t64 * t25 + t61 * t26) - t86 * t12 + t71, qJ(3) * (t64 * t28 - t61 * t31) - t86 * t18 - t69, qJ(3) * (t64 * t2 + t61 * t7) - t86 * t1 - t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), 0, -t67, t19, 0, 0, 0, 0, 0, 0, -t30, t29, 0, t3, 0, 0, 0, 0, 0, 0, t11, t12, t18, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t78, -t10, 0, 0, t21, t17, t23, -t20, t22, 0, t72, -t71, t69, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t45 - t46, -t84, t41, -t43, qJDD(5), -t5, -t6, 0, 0;];
tauJ_reg = t4;
