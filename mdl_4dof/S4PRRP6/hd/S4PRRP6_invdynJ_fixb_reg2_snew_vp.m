% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4PRRP6
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
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4PRRP6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:43
% EndTime: 2019-12-31 16:30:44
% DurationCPUTime: 0.40s
% Computational Cost: add. (421->92), mult. (890->112), div. (0->0), fcn. (489->6), ass. (0->67)
t66 = qJD(3) ^ 2;
t62 = sin(qJ(3));
t55 = t62 ^ 2;
t67 = qJD(2) ^ 2;
t94 = t55 * t67;
t43 = t66 + t94;
t64 = cos(qJ(3));
t47 = t62 * t67 * t64;
t42 = qJDD(3) - t47;
t89 = t64 * t42;
t19 = -t62 * t43 + t89;
t98 = pkin(5) * t19;
t78 = qJD(2) * qJD(3);
t80 = t62 * qJDD(2);
t33 = 0.2e1 * t64 * t78 + t80;
t63 = sin(qJ(2));
t65 = cos(qJ(2));
t97 = t63 * t19 + t65 * t33;
t56 = t64 ^ 2;
t93 = t56 * t67;
t96 = t89 + t62 * (-t66 + t93);
t59 = sin(pkin(6));
t60 = cos(pkin(6));
t40 = -t60 * g(1) - t59 * g(2);
t81 = -g(3) + qJDD(1);
t22 = t65 * t40 + t63 * t81;
t13 = -t67 * pkin(2) + qJDD(2) * pkin(5) + t22;
t39 = -t59 * g(1) + t60 * g(2);
t28 = t64 * t39;
t84 = t62 * qJ(4);
t72 = -t64 * pkin(3) - t84;
t83 = t67 * t72;
t5 = -qJDD(3) * pkin(3) - t66 * qJ(4) + (t13 + t83) * t62 + qJDD(4) - t28;
t95 = 2 * qJD(4);
t41 = qJDD(3) + t47;
t92 = t62 * t41;
t9 = t64 * t13 + t62 * t39;
t45 = -t66 - t93;
t18 = t64 * t45 - t92;
t77 = t62 * t78;
t79 = t64 * qJDD(2);
t34 = -0.2e1 * t77 + t79;
t87 = pkin(2) * t34 + pkin(5) * t18;
t85 = t55 + t56;
t36 = t85 * qJDD(2);
t37 = t85 * t67;
t86 = pkin(2) * t37 + pkin(5) * t36;
t82 = qJD(2) * t62;
t7 = t62 * t13 - t28;
t2 = t62 * t7 + t64 * t9;
t21 = -t63 * t40 + t65 * t81;
t73 = qJDD(3) * qJ(4) + (qJD(3) * t95) + t64 * t83 + t9;
t71 = t64 * t33 + t62 * t34;
t12 = -qJDD(2) * pkin(2) - t67 * pkin(5) - t21;
t70 = pkin(2) - t72;
t69 = t12 - (-t77 + t79) * pkin(3) - qJ(4) * t33;
t68 = t82 * t95 - t69;
t38 = (t55 - t56) * t67;
t17 = t92 + t64 * (t66 - t94);
t16 = t33 * t62;
t15 = t34 * t64;
t11 = t63 * t36 + t65 * t37;
t8 = t63 * t18 + t65 * t34;
t4 = -t66 * pkin(3) + t73;
t3 = (pkin(3) * qJD(3) - (2 * qJD(4))) * t82 + t69;
t1 = t64 * t4 + t62 * t5;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t81, 0, 0, 0, 0, 0, 0, t65 * qJDD(2) - t63 * t67, -t63 * qJDD(2) - t65 * t67, 0, t65 * t21 + t63 * t22, 0, 0, 0, 0, 0, 0, t8, -t97, t11, -t65 * t12 + t63 * t2, 0, 0, 0, 0, 0, 0, t8, t11, t97, t63 * t1 - t65 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t21, -t22, 0, 0, t16, t71, t17, t15, t96, 0, -t64 * t12 + t87, -pkin(2) * t33 + t62 * t12 - t98, t2 + t86, -pkin(2) * t12 + pkin(5) * t2, t16, t17, -t71, 0, -t96, t15, t34 * t84 + t64 * ((t34 - t77) * pkin(3) + t68) + t87, t64 * ((t37 - t66) * pkin(3) + t73) + (qJ(4) * t37 + t5) * t62 + t86, t62 * (-pkin(3) * t77 + t68) + t98 + t70 * t33, pkin(5) * t1 - t70 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t38, t80, t47, t79, qJDD(3), -t7, -t9, 0, 0, -t47, t80, -t38, qJDD(3), -t79, t47, pkin(3) * t41 + qJ(4) * t45 - t5, (-pkin(3) * t62 + qJ(4) * t64) * qJDD(2), qJ(4) * t42 + (t43 - t66) * pkin(3) + t73, -pkin(3) * t5 + qJ(4) * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t80, -t43, t5;];
tauJ_reg = t6;
