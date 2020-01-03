% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
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
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4PRRP5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:29:22
% EndTime: 2019-12-31 16:29:23
% DurationCPUTime: 0.34s
% Computational Cost: add. (426->92), mult. (922->114), div. (0->0), fcn. (505->6), ass. (0->71)
t79 = (qJD(2) * qJD(4));
t94 = 2 * t79;
t67 = sin(qJ(3));
t69 = cos(qJ(3));
t72 = qJD(2) ^ 2;
t50 = t67 * t72 * t69;
t44 = qJDD(3) + t50;
t93 = pkin(3) * t44;
t64 = sin(pkin(6));
t65 = cos(pkin(6));
t42 = -t65 * g(1) - t64 * g(2);
t61 = -g(3) + qJDD(1);
t68 = sin(qJ(2));
t70 = cos(qJ(2));
t25 = t70 * t42 + t68 * t61;
t15 = -t72 * pkin(2) + qJDD(2) * pkin(5) + t25;
t41 = -t64 * g(1) + t65 * g(2);
t10 = t69 * t15 + t67 * t41;
t58 = t69 * qJDD(2);
t80 = qJD(2) * qJD(3);
t77 = t67 * t80;
t35 = t58 - t77;
t82 = qJD(2) * t67;
t43 = qJD(3) * pkin(3) - qJ(4) * t82;
t73 = t35 * qJ(4) - qJD(3) * t43 + t69 * t94 + t10;
t59 = t67 ^ 2;
t92 = t59 * t72;
t60 = t69 ^ 2;
t91 = t60 * t72;
t90 = t67 * t15;
t89 = t67 * t44;
t45 = qJDD(3) - t50;
t88 = t69 * t45;
t71 = qJD(3) ^ 2;
t48 = -t71 - t91;
t22 = t69 * t48 - t89;
t36 = t58 - 0.2e1 * t77;
t87 = pkin(2) * t36 + pkin(5) * t22;
t47 = -t71 - t92;
t23 = -t67 * t47 - t88;
t57 = t67 * qJDD(2);
t76 = t69 * t80;
t33 = t57 + 0.2e1 * t76;
t86 = -pkin(2) * t33 + pkin(5) * t23;
t84 = t59 + t60;
t38 = t84 * qJDD(2);
t39 = t84 * t72;
t85 = pkin(2) * t39 + pkin(5) * t38;
t83 = qJ(4) * t67;
t81 = qJ(4) * qJDD(2);
t31 = t69 * t41;
t7 = -t31 + t90;
t2 = t69 * t10 + t67 * t7;
t24 = -t68 * t42 + t70 * t61;
t14 = -qJDD(2) * pkin(2) - t72 * pkin(5) - t24;
t34 = t57 + t76;
t74 = -t31 + (t34 - t76) * qJ(4) - t93;
t3 = -0.2e1 * t67 * t79 - t74 - t90;
t5 = -t35 * pkin(3) - qJ(4) * t91 + t43 * t82 + qJDD(4) + t14;
t40 = (t59 - t60) * t72;
t21 = t89 + t69 * (t71 - t92);
t20 = t67 * (-t71 + t91) + t88;
t19 = (t34 + t76) * t67;
t18 = (t35 - t77) * t69;
t13 = t68 * t38 + t70 * t39;
t11 = t69 * t33 + t67 * t36;
t9 = t68 * t23 - t70 * t33;
t8 = t68 * t22 + t70 * t36;
t4 = -pkin(3) * t91 + t73;
t1 = -t67 * t3 + t69 * t4;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, 0, 0, 0, 0, 0, t70 * qJDD(2) - t68 * t72, -t68 * qJDD(2) - t70 * t72, 0, t70 * t24 + t68 * t25, 0, 0, 0, 0, 0, 0, t8, t9, t13, -t70 * t14 + t68 * t2, 0, 0, 0, 0, 0, 0, t8, t9, t13, t68 * t1 - t70 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t24, -t25, 0, 0, t19, t11, t21, t18, t20, 0, -t69 * t14 + t87, t67 * t14 + t86, t2 + t85, -pkin(2) * t14 + pkin(5) * t2, t19, t11, t21, t18, t20, 0, -t44 * t83 + t69 * (pkin(3) * t36 + qJ(4) * t48 - t5) + t87, t67 * (-qJ(4) * t47 + t5) + t69 * (-pkin(3) * t33 - qJ(4) * t45) + t86, t69 * (t69 * t81 + (t39 - t91) * pkin(3) + t73) + ((t15 + t94 + t81) * t67 + t74) * t67 + t85, -t3 * t83 + t69 * (-pkin(3) * t5 + qJ(4) * t4) - pkin(2) * t5 + pkin(5) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t40, t57, t50, t58, qJDD(3), -t7, -t10, 0, 0, -t50, t40, t57, t50, t58, qJDD(3), t3 + t93, (t47 + t91) * pkin(3) - t73, -pkin(3) * t57, pkin(3) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t33, -t39, t5;];
tauJ_reg = t6;
