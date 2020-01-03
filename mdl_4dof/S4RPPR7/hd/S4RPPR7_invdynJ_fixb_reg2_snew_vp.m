% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RPPR7
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RPPR7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:40
% EndTime: 2019-12-31 16:41:43
% DurationCPUTime: 0.62s
% Computational Cost: add. (927->118), mult. (2082->163), div. (0->0), fcn. (1269->6), ass. (0->83)
t60 = sin(pkin(6));
t61 = cos(pkin(6));
t62 = sin(qJ(4));
t64 = cos(qJ(4));
t74 = t60 * t64 + t61 * t62;
t44 = t74 * qJD(1);
t93 = t60 * t62;
t46 = (t61 * t64 - t93) * qJD(1);
t94 = t46 * t44;
t100 = qJDD(4) - t94;
t102 = t100 * t62;
t101 = t100 * t64;
t67 = qJD(1) ^ 2;
t63 = sin(qJ(1));
t65 = cos(qJ(1));
t79 = t63 * g(1) - t65 * g(2);
t75 = qJDD(2) - t79;
t72 = -t67 * qJ(2) + t75;
t88 = -qJ(3) - pkin(1);
t77 = -0.2e1 * qJD(1) * qJD(3) + t88 * qJDD(1) + t72;
t57 = t60 ^ 2;
t58 = t61 ^ 2;
t87 = t57 + t58;
t99 = t87 * t67;
t18 = t74 * qJDD(1);
t41 = t44 ^ 2;
t42 = t46 ^ 2;
t97 = pkin(3) * t67;
t96 = t60 * g(3);
t23 = -t61 * g(3) + t77 * t60;
t83 = t60 * qJDD(1);
t20 = -pkin(5) * t83 - t57 * t97 + t23;
t70 = (-pkin(5) * qJDD(1) - t60 * t97 + t77) * t61;
t8 = t62 * t20 - t64 * (t70 + t96);
t9 = g(3) * t93 + t64 * t20 + t62 * t70;
t2 = t62 * t9 - t64 * t8;
t95 = t61 * t2;
t59 = qJDD(1) * qJ(2);
t76 = t65 * g(1) + t63 * g(2);
t73 = -t59 + t76;
t81 = qJD(2) * qJD(1);
t71 = -qJDD(3) + t73 - 0.2e1 * t81;
t21 = -pkin(3) * t83 + (t87 * pkin(5) - t88) * t67 + t71;
t92 = t62 * t21;
t26 = qJDD(4) + t94;
t91 = t62 * t26;
t90 = t64 * t21;
t89 = t64 * t26;
t86 = qJDD(1) * pkin(1);
t85 = t44 * qJD(4);
t84 = t46 * qJD(4);
t82 = t61 * qJDD(1);
t3 = t62 * t8 + t64 * t9;
t34 = -t88 * t67 + t71;
t78 = -t34 + t59;
t43 = -t62 * t83 + t64 * t82;
t10 = t61 * (t77 * t61 + t96) + t60 * t23;
t66 = qJD(4) ^ 2;
t55 = 0.2e1 * t81;
t49 = t87 * qJDD(1);
t48 = t60 * t99;
t47 = t61 * t99;
t40 = -t72 + t86;
t38 = -t42 - t66;
t37 = -t42 + t66;
t36 = t41 - t66;
t31 = t43 - t85;
t30 = t43 - 0.2e1 * t85;
t29 = -t18 - t84;
t28 = 0.2e1 * t84 + t18;
t24 = -t66 - t41;
t19 = -t41 - t42;
t16 = -t62 * t38 - t89;
t15 = t64 * t38 - t91;
t14 = -t64 * t18 + t62 * t43;
t13 = -t62 * t18 - t64 * t43;
t12 = t64 * t24 - t102;
t11 = t62 * t24 + t101;
t6 = t61 * t15 + t60 * t16;
t5 = t61 * t13 + t60 * t14;
t4 = t61 * t11 + t60 * t12;
t1 = t60 * t3 + t95;
t7 = [0, 0, 0, 0, 0, qJDD(1), t79, t76, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t75 - 0.2e1 * t86, t55 + 0.2e1 * t59 - t76, pkin(1) * t40 + qJ(2) * (-t67 * pkin(1) + t55 - t73), t58 * qJDD(1), -0.2e1 * t60 * t82, 0, t57 * qJDD(1), 0, 0, -t88 * t48 + t78 * t60, -t88 * t47 + t78 * t61, -qJ(2) * t99 - t88 * t49 - t10, -qJ(2) * t34 + t88 * t10, t61 * (t64 * t31 - t62 * t84) - t60 * (t62 * t31 + t64 * t84), t61 * (-t64 * t28 - t62 * t30) - t60 * (-t62 * t28 + t64 * t30), t61 * (-t62 * t37 + t101) - t60 * (t64 * t37 + t102), t61 * (-t62 * t29 + t64 * t85) - t60 * (t64 * t29 + t62 * t85), t61 * (t64 * t36 - t91) - t60 * (t62 * t36 + t89), (t61 * (-t44 * t64 + t46 * t62) - t60 * (-t44 * t62 - t46 * t64)) * qJD(4), t61 * (-pkin(5) * t11 - t92) - t60 * (-pkin(3) * t28 + pkin(5) * t12 + t90) + qJ(2) * t28 + t88 * t4, t61 * (-pkin(5) * t15 - t90) - t60 * (-pkin(3) * t30 + pkin(5) * t16 - t92) + qJ(2) * t30 + t88 * t6, t61 * (-pkin(5) * t13 - t2) - t60 * (-pkin(3) * t19 + pkin(5) * t14 + t3) + qJ(2) * t19 + t88 * t5, -pkin(5) * t95 - t60 * (pkin(3) * t21 + pkin(5) * t3) - qJ(2) * t21 + t88 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t67, -t40, 0, 0, 0, 0, 0, 0, -t48, -t47, -t49, t10, 0, 0, 0, 0, 0, 0, t4, t6, t5, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t82, -t99, -t34, 0, 0, 0, 0, 0, 0, t28, t30, t19, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t42 - t41, t43, -t94, -t18, qJDD(4), -t8, -t9, 0, 0;];
tauJ_reg = t7;
