% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RPRR2
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RPRR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:13
% EndTime: 2019-12-31 16:48:14
% DurationCPUTime: 0.32s
% Computational Cost: add. (913->87), mult. (1427->129), div. (0->0), fcn. (824->8), ass. (0->66)
t66 = sin(qJ(4));
t69 = cos(qJ(4));
t73 = qJD(1) ^ 2;
t68 = sin(qJ(1));
t71 = cos(qJ(1));
t79 = t68 * g(1) - t71 * g(2);
t43 = qJDD(1) * pkin(1) + t79;
t76 = t71 * g(1) + t68 * g(2);
t44 = -t73 * pkin(1) - t76;
t64 = sin(pkin(7));
t65 = cos(pkin(7));
t83 = t64 * t43 + t65 * t44;
t21 = -t73 * pkin(2) + t83;
t67 = sin(qJ(3));
t70 = cos(qJ(3));
t77 = t65 * t43 - t64 * t44;
t74 = qJDD(1) * pkin(2) + t77;
t15 = t70 * t21 + t67 * t74;
t60 = qJD(1) + qJD(3);
t58 = t60 ^ 2;
t59 = qJDD(1) + qJDD(3);
t13 = -t58 * pkin(3) + t59 * pkin(6) + t15;
t63 = -g(3) + qJDD(2);
t7 = t66 * t13 - t69 * t63;
t8 = t69 * t13 + t66 * t63;
t3 = t66 * t7 + t69 * t8;
t14 = -t67 * t21 + t70 * t74;
t12 = -t59 * pkin(3) - t58 * pkin(6) - t14;
t87 = -pkin(3) * t12 + pkin(6) * t3;
t49 = t69 * t58 * t66;
t45 = qJDD(4) + t49;
t86 = t66 * t45;
t85 = t66 * t59;
t46 = qJDD(4) - t49;
t84 = t69 * t46;
t82 = qJD(4) * t60;
t61 = t66 ^ 2;
t52 = t61 * t58;
t72 = qJD(4) ^ 2;
t47 = -t52 - t72;
t29 = -t66 * t47 - t84;
t33 = 0.2e1 * t69 * t82 + t85;
t81 = -pkin(3) * t33 + pkin(6) * t29 + t66 * t12;
t62 = t69 ^ 2;
t53 = t62 * t58;
t48 = -t53 - t72;
t28 = t69 * t48 - t86;
t51 = t69 * t59;
t34 = -0.2e1 * t66 * t82 + t51;
t80 = pkin(3) * t34 + pkin(6) * t28 - t69 * t12;
t39 = (t61 + t62) * t59;
t42 = t52 + t53;
t78 = pkin(3) * t42 + pkin(6) * t39 + t3;
t40 = -t70 * t58 - t67 * t59;
t75 = t67 * t58 - t70 * t59;
t27 = t86 + t69 * (-t52 + t72);
t26 = t66 * (t53 - t72) + t84;
t23 = t33 * t66;
t22 = t34 * t69;
t20 = t67 * t39 + t70 * t42;
t18 = t69 * t33 + t66 * t34;
t17 = t67 * t29 - t70 * t33;
t16 = t67 * t28 + t70 * t34;
t4 = t70 * t14 + t67 * t15;
t1 = -t70 * t12 + t67 * t3;
t2 = [0, 0, 0, 0, 0, qJDD(1), t79, t76, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t65 * qJDD(1) - t64 * t73) + t77, pkin(1) * (-t64 * qJDD(1) - t65 * t73) - t83, 0, pkin(1) * (t64 * t83 + t65 * t77), 0, 0, 0, 0, 0, t59, pkin(1) * (t64 * t40 - t65 * t75) - pkin(2) * t75 + t14, pkin(1) * (t65 * t40 + t64 * t75) + pkin(2) * t40 - t15, 0, pkin(1) * (t64 * (-t67 * t14 + t70 * t15) + t65 * t4) + pkin(2) * t4, t23, t18, t27, t22, t26, 0, pkin(1) * (t64 * (t70 * t28 - t67 * t34) + t65 * t16) + pkin(2) * t16 + t80, pkin(1) * (t64 * (t70 * t29 + t67 * t33) + t65 * t17) + pkin(2) * t17 + t81, pkin(1) * (t64 * (t70 * t39 - t67 * t42) + t65 * t20) + pkin(2) * t20 + t78, pkin(1) * (t64 * (t67 * t12 + t70 * t3) + t65 * t1) + pkin(2) * t1 + t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, 0, 0, 0, 0, 0, t69 * t45 + t66 * t48, -t66 * t46 + t69 * t47, 0, t66 * t8 - t69 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t14, -t15, 0, 0, t23, t18, t27, t22, t26, 0, t80, t81, t78, t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t52 - t53, t85, t49, t51, qJDD(4), -t7, -t8, 0, 0;];
tauJ_reg = t2;
