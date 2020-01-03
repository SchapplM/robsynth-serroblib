% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RPRP3
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RPRP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:47
% EndTime: 2019-12-31 16:42:49
% DurationCPUTime: 0.42s
% Computational Cost: add. (566->104), mult. (1190->132), div. (0->0), fcn. (604->6), ass. (0->72)
t88 = (qJD(1) * qJD(4));
t101 = 2 * t88;
t69 = sin(qJ(3));
t71 = cos(qJ(3));
t74 = qJD(1) ^ 2;
t50 = t71 * t74 * t69;
t44 = qJDD(3) + t50;
t100 = pkin(3) * t44;
t59 = t71 * qJDD(1);
t89 = qJD(1) * qJD(3);
t82 = t69 * t89;
t37 = t59 - t82;
t91 = qJD(1) * t69;
t43 = qJD(3) * pkin(3) - qJ(4) * t91;
t70 = sin(qJ(1));
t72 = cos(qJ(1));
t79 = t72 * g(1) + t70 * g(2);
t33 = -t74 * pkin(1) - t79;
t66 = sin(pkin(6));
t67 = cos(pkin(6));
t84 = t70 * g(1) - t72 * g(2);
t77 = qJDD(1) * pkin(1) + t84;
t94 = t67 * t33 + t66 * t77;
t13 = -t74 * pkin(2) + qJDD(1) * pkin(5) + t94;
t63 = -g(3) + qJDD(2);
t8 = t71 * t13 + t69 * t63;
t75 = t37 * qJ(4) - qJD(3) * t43 + t71 * t101 + t8;
t61 = t69 ^ 2;
t99 = t61 * t74;
t62 = t71 ^ 2;
t98 = t62 * t74;
t97 = t69 * t13;
t96 = t69 * t44;
t45 = qJDD(3) - t50;
t95 = t71 * t45;
t93 = t61 + t62;
t92 = qJ(4) * t69;
t58 = t69 * qJDD(1);
t90 = qJ(4) * qJDD(1);
t73 = qJD(3) ^ 2;
t48 = -t73 - t98;
t24 = t71 * t48 - t96;
t38 = t59 - 0.2e1 * t82;
t87 = pkin(5) * t24 + pkin(2) * t38 + pkin(1) * (t66 * t24 + t67 * t38);
t47 = -t73 - t99;
t25 = -t69 * t47 - t95;
t81 = t71 * t89;
t35 = t58 + 0.2e1 * t81;
t86 = pkin(1) * (t66 * t25 - t67 * t35) + pkin(5) * t25 - pkin(2) * t35;
t40 = t93 * qJDD(1);
t41 = t93 * t74;
t85 = pkin(2) * t41 + pkin(5) * t40 + pkin(1) * (t66 * t40 + t67 * t41);
t56 = t71 * t63;
t7 = -t56 + t97;
t2 = t69 * t7 + t71 * t8;
t80 = -t66 * t33 + t67 * t77;
t12 = -qJDD(1) * pkin(2) - t74 * pkin(5) - t80;
t36 = t58 + t81;
t76 = -t56 + (t36 - t81) * qJ(4) - t100;
t3 = -0.2e1 * t69 * t88 - t76 - t97;
t5 = -t37 * pkin(3) - qJ(4) * t98 + t43 * t91 + qJDD(4) + t12;
t42 = (t61 - t62) * t74;
t23 = -t69 * t45 + t71 * t47;
t22 = t96 + t71 * (t73 - t99);
t21 = t71 * t44 + t69 * t48;
t20 = t69 * (-t73 + t98) + t95;
t19 = (t36 + t81) * t69;
t18 = (t37 - t82) * t71;
t14 = t71 * t35 + t69 * t38;
t4 = -pkin(3) * t98 + t75;
t1 = -t69 * t3 + t71 * t4;
t6 = [0, 0, 0, 0, 0, qJDD(1), t84, t79, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t67 * qJDD(1) - t66 * t74) + t80, pkin(1) * (-t66 * qJDD(1) - t67 * t74) - t94, 0, pkin(1) * (t66 * t94 + t67 * t80), t19, t14, t22, t18, t20, 0, -t71 * t12 + t87, t69 * t12 + t86, t2 + t85, -pkin(2) * t12 + pkin(5) * t2 + pkin(1) * (-t67 * t12 + t66 * t2), t19, t14, t22, t18, t20, 0, -t44 * t92 + t71 * (pkin(3) * t38 + qJ(4) * t48 - t5) + t87, t69 * (-qJ(4) * t47 + t5) + t71 * (-pkin(3) * t35 - qJ(4) * t45) + t86, t71 * (t71 * t90 + (t41 - t98) * pkin(3) + t75) + ((t13 + t101 + t90) * t69 + t76) * t69 + t85, -t3 * t92 + t71 * (-pkin(3) * t5 + qJ(4) * t4) - pkin(2) * t5 + pkin(5) * t1 + pkin(1) * (t66 * t1 - t67 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, 0, 0, 0, 0, 0, t21, t23, 0, t69 * t8 - t71 * t7, 0, 0, 0, 0, 0, 0, t21, t23, 0, t71 * t3 + t69 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t42, t58, t50, t59, qJDD(3), -t7, -t8, 0, 0, -t50, t42, t58, t50, t59, qJDD(3), t3 + t100, (t47 + t98) * pkin(3) - t75, -pkin(3) * t58, pkin(3) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t35, -t41, t5;];
tauJ_reg = t6;
