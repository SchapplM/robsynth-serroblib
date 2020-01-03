% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RPRP4
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
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RPRP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:53
% EndTime: 2019-12-31 16:43:55
% DurationCPUTime: 0.44s
% Computational Cost: add. (559->104), mult. (1148->130), div. (0->0), fcn. (583->6), ass. (0->66)
t64 = cos(qJ(3));
t86 = qJD(1) * qJD(3);
t62 = sin(qJ(3));
t88 = t62 * qJDD(1);
t34 = 0.2e1 * t64 * t86 + t88;
t67 = qJD(1) ^ 2;
t66 = qJD(3) ^ 2;
t46 = t64 * t67 * t62;
t41 = qJDD(3) - t46;
t95 = t64 * t41;
t56 = t64 ^ 2;
t99 = t56 * t67;
t103 = t95 + t62 * (-t66 + t99);
t55 = t62 ^ 2;
t100 = t55 * t67;
t42 = t66 + t100;
t102 = -t62 * t42 + t95;
t63 = sin(qJ(1));
t65 = cos(qJ(1));
t76 = t65 * g(1) + t63 * g(2);
t32 = -t67 * pkin(1) - t76;
t59 = sin(pkin(6));
t60 = cos(pkin(6));
t75 = t63 * g(1) - t65 * g(2);
t71 = qJDD(1) * pkin(1) + t75;
t94 = t60 * t32 + t59 * t71;
t12 = -t67 * pkin(2) + qJDD(1) * pkin(5) + t94;
t89 = -g(3) + qJDD(2);
t51 = t64 * t89;
t92 = t62 * qJ(4);
t74 = -t64 * pkin(3) - t92;
t91 = t67 * t74;
t5 = -qJDD(3) * pkin(3) - t66 * qJ(4) + (t12 + t91) * t62 + qJDD(4) - t51;
t101 = 2 * qJD(4);
t40 = qJDD(3) + t46;
t98 = t62 * t40;
t8 = t64 * t12 + t62 * t89;
t93 = t55 + t56;
t90 = qJD(1) * t62;
t87 = t64 * qJDD(1);
t44 = -t66 - t99;
t19 = t64 * t44 - t98;
t82 = t62 * t86;
t35 = -0.2e1 * t82 + t87;
t85 = pkin(5) * t19 + pkin(2) * t35 + pkin(1) * (t59 * t19 + t60 * t35);
t37 = t93 * qJDD(1);
t38 = t93 * t67;
t84 = pkin(1) * (t59 * t37 + t60 * t38) + pkin(5) * t37 + pkin(2) * t38;
t83 = pkin(1) * t59 + pkin(5);
t7 = t62 * t12 - t51;
t2 = t62 * t7 + t64 * t8;
t80 = -t59 * t32 + t60 * t71;
t77 = qJDD(3) * qJ(4) + (qJD(3) * t101) + t64 * t91 + t8;
t73 = t64 * t34 + t62 * t35;
t72 = t62 * t41 + t64 * t42;
t11 = -qJDD(1) * pkin(2) - t67 * pkin(5) - t80;
t70 = pkin(1) * t60 + pkin(2) - t74;
t69 = t11 - (-t82 + t87) * pkin(3) - qJ(4) * t34;
t68 = t90 * t101 - t69;
t39 = (t55 - t56) * t67;
t18 = t98 + t64 * (t66 - t100);
t17 = t64 * t40 + t62 * t44;
t16 = t34 * t62;
t15 = t35 * t64;
t4 = -t66 * pkin(3) + t77;
t1 = [0, 0, 0, 0, 0, qJDD(1), t75, t76, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t60 * qJDD(1) - t59 * t67) + t80, pkin(1) * (-t59 * qJDD(1) - t60 * t67) - t94, 0, pkin(1) * (t59 * t94 + t60 * t80), t16, t73, t18, t15, t103, 0, -t64 * t11 + t85, t62 * t11 - pkin(2) * t34 - pkin(5) * t102 + pkin(1) * (-t102 * t59 - t60 * t34), t2 + t84, -pkin(2) * t11 + pkin(5) * t2 + pkin(1) * (-t60 * t11 + t59 * t2), t16, t18, -t73, 0, -t103, t15, t35 * t92 + t64 * ((t35 - t82) * pkin(3) + t68) + t85, t64 * ((t38 - t66) * pkin(3) + t77) + (qJ(4) * t38 + t5) * t62 + t84, t62 * (-pkin(3) * t82 + t68) + t83 * t102 + t70 * t34, t83 * (t64 * t4 + t62 * t5) - t70 * ((pkin(3) * qJD(3) - (2 * qJD(4))) * t90 + t69); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, 0, 0, 0, 0, 0, 0, t17, -t72, 0, t62 * t8 - t64 * t7, 0, 0, 0, 0, 0, 0, t17, 0, t72, t62 * t4 - t64 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t39, t88, t46, t87, qJDD(3), -t7, -t8, 0, 0, -t46, t88, -t39, qJDD(3), -t87, t46, pkin(3) * t40 + qJ(4) * t44 - t5, (-pkin(3) * t62 + qJ(4) * t64) * qJDD(1), qJ(4) * t41 + (t42 - t66) * pkin(3) + t77, -pkin(3) * t5 + qJ(4) * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, t88, -t42, t5;];
tauJ_reg = t1;
