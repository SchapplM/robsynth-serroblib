% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4PRRP3
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
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4PRRP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:56
% EndTime: 2019-12-31 16:26:57
% DurationCPUTime: 0.32s
% Computational Cost: add. (386->87), mult. (847->105), div. (0->0), fcn. (466->6), ass. (0->67)
t73 = (qJD(2) * qJD(4));
t89 = 2 * t73;
t58 = sin(qJ(3));
t60 = cos(qJ(3));
t63 = qJD(2) ^ 2;
t41 = t58 * t63 * t60;
t35 = qJDD(3) + t41;
t88 = pkin(3) * t35;
t49 = t60 * qJDD(2);
t74 = qJD(2) * qJD(3);
t70 = t58 * t74;
t27 = t49 - t70;
t77 = qJD(2) * t58;
t34 = qJD(3) * pkin(3) - qJ(4) * t77;
t52 = -g(3) + qJDD(1);
t55 = sin(pkin(6));
t56 = cos(pkin(6));
t32 = t55 * g(1) - t56 * g(2);
t33 = -t56 * g(1) - t55 * g(2);
t59 = sin(qJ(2));
t61 = cos(qJ(2));
t67 = -t59 * t32 - t61 * t33;
t76 = qJDD(2) * pkin(5);
t9 = -t63 * pkin(2) - t67 + t76;
t6 = t58 * t52 + t60 * t9;
t64 = t27 * qJ(4) - qJD(3) * t34 + t60 * t89 + t6;
t87 = t58 * t9;
t50 = t58 ^ 2;
t86 = t50 * t63;
t51 = t60 ^ 2;
t85 = t51 * t63;
t84 = t58 * t35;
t36 = qJDD(3) - t41;
t83 = t60 * t36;
t28 = t49 - 0.2e1 * t70;
t62 = qJD(3) ^ 2;
t39 = -t62 - t85;
t82 = pkin(5) * (t60 * t39 - t84) + pkin(2) * t28;
t48 = t58 * qJDD(2);
t69 = t60 * t74;
t25 = t48 + 0.2e1 * t69;
t38 = -t62 - t86;
t81 = pkin(5) * (-t58 * t38 - t83) - pkin(2) * t25;
t79 = t50 + t51;
t30 = t79 * t63;
t80 = pkin(2) * t30 + t79 * t76;
t78 = qJ(4) * t58;
t75 = qJ(4) * qJDD(2);
t46 = t60 * t52;
t5 = -t46 + t87;
t72 = t58 * t5 + t60 * t6;
t68 = t61 * t32 - t59 * t33;
t8 = -qJDD(2) * pkin(2) - t63 * pkin(5) - t68;
t26 = t48 + t69;
t65 = -t46 + (t26 - t69) * qJ(4) - t88;
t1 = -0.2e1 * t58 * t73 - t65 - t87;
t3 = -t27 * pkin(3) - qJ(4) * t85 + t34 * t77 + qJDD(4) + t8;
t31 = (t50 - t51) * t63;
t18 = -t58 * t36 + t60 * t38;
t17 = t84 + t60 * (t62 - t86);
t16 = t60 * t35 + t58 * t39;
t15 = t58 * (-t62 + t85) + t83;
t14 = (t26 + t69) * t58;
t13 = (t27 - t70) * t60;
t10 = t60 * t25 + t58 * t28;
t2 = -pkin(3) * t85 + t64;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, 0, 0, 0, 0, 0, t16, t18, 0, -t60 * t5 + t58 * t6, 0, 0, 0, 0, 0, 0, t16, t18, 0, t60 * t1 + t58 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t68, t67, 0, 0, t14, t10, t17, t13, t15, 0, -t60 * t8 + t82, t58 * t8 + t81, t72 + t80, -pkin(2) * t8 + pkin(5) * t72, t14, t10, t17, t13, t15, 0, -t35 * t78 + t60 * (pkin(3) * t28 + qJ(4) * t39 - t3) + t82, t58 * (-qJ(4) * t38 + t3) + t60 * (-pkin(3) * t25 - qJ(4) * t36) + t81, t60 * (t60 * t75 + (t30 - t85) * pkin(3) + t64) + ((t89 + t9 + t75) * t58 + t65) * t58 + t80, -t1 * t78 + t60 * (-pkin(3) * t3 + qJ(4) * t2) - pkin(2) * t3 + pkin(5) * (-t58 * t1 + t60 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t31, t48, t41, t49, qJDD(3), -t5, -t6, 0, 0, -t41, t31, t48, t41, t49, qJDD(3), t1 + t88, (t38 + t85) * pkin(3) - t64, -pkin(3) * t48, pkin(3) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t25, -t30, t3;];
tauJ_reg = t4;
