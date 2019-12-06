% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRPRP6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRPRP6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:41:18
% EndTime: 2019-12-05 15:41:21
% DurationCPUTime: 0.48s
% Computational Cost: add. (656->115), mult. (1253->125), div. (0->0), fcn. (655->6), ass. (0->81)
t107 = pkin(6) + pkin(2);
t70 = sin(qJ(4));
t72 = cos(qJ(4));
t75 = qJD(2) ^ 2;
t49 = t72 * t75 * t70;
t44 = qJDD(4) + t49;
t101 = t70 * t44;
t62 = t72 ^ 2;
t102 = t62 * t75;
t74 = qJD(4) ^ 2;
t48 = t74 + t102;
t19 = t72 * t48 + t101;
t106 = t107 * t19;
t92 = qJD(2) * qJD(4);
t90 = t70 * t92;
t93 = t72 * qJDD(2);
t36 = -0.2e1 * t90 + t93;
t71 = sin(qJ(2));
t73 = cos(qJ(2));
t105 = t73 * t19 + t71 * t36;
t65 = qJDD(2) * pkin(2);
t66 = sin(pkin(7));
t67 = cos(pkin(7));
t43 = -t67 * g(1) - t66 * g(2);
t63 = -g(3) + qJDD(1);
t23 = -t71 * t43 + t73 * t63;
t84 = qJDD(3) - t23;
t14 = -t75 * qJ(3) - t65 + t84;
t12 = -qJDD(2) * pkin(6) + t14;
t42 = -t66 * g(1) + t67 * g(2);
t7 = t70 * t12 + t72 * t42;
t61 = t70 ^ 2;
t103 = t61 * t75;
t45 = qJDD(4) - t49;
t99 = t72 * t45;
t24 = t73 * t43 + t71 * t63;
t98 = t61 + t62;
t97 = t72 * qJ(5);
t82 = t70 * pkin(4) - t97;
t96 = t75 * t82;
t95 = qJD(5) * t72;
t94 = t70 * qJDD(2);
t91 = (qJD(3) * qJD(2));
t89 = t72 * t92;
t6 = -t72 * t12 + t70 * t42;
t58 = qJDD(2) * qJ(3);
t88 = -t75 * pkin(2) + t24 + t58;
t47 = -t74 - t103;
t18 = t70 * t47 + t99;
t33 = 0.2e1 * t89 + t94;
t87 = qJ(3) * t33 - t107 * t18;
t37 = t98 * qJDD(2);
t41 = t98 * t75;
t86 = -qJ(3) * t41 + t107 * t37;
t85 = -t75 * pkin(6) + t88;
t2 = -t72 * t6 + t70 * t7;
t83 = pkin(4) * t72 + qJ(5) * t70;
t81 = t72 * t33 + t70 * t36;
t80 = -t72 * (-t74 + t103) + t101;
t79 = qJ(3) + t82;
t34 = -t89 - t94;
t35 = -t90 + t93;
t78 = -t34 * pkin(4) - t35 * qJ(5) + t85;
t77 = qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t70 * t96 + t7;
t76 = 0.2e1 * qJD(2) * t95 - t78 - (2 * t91);
t5 = -qJDD(4) * pkin(4) - t74 * qJ(5) + t72 * t96 + qJDD(5) + t6;
t55 = 2 * t91;
t40 = (-t61 + t62) * t75;
t39 = -t73 * qJDD(2) + t71 * t75;
t38 = t71 * qJDD(2) + t73 * t75;
t22 = t99 - t70 * (t74 - t102);
t21 = (t35 - t90) * t72;
t17 = (-t34 + t89) * t70;
t15 = t73 * t37 - t71 * t41;
t13 = t55 + t88;
t11 = t55 + t85;
t8 = -t73 * t18 + t71 * t33;
t4 = -t74 * pkin(4) + t77;
t3 = t55 + (t83 * qJD(4) - 0.2e1 * t95) * qJD(2) + t78;
t1 = t70 * t4 - t72 * t5;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, 0, 0, 0, 0, 0, -t39, -t38, 0, t73 * t23 + t71 * t24, 0, 0, 0, 0, 0, 0, 0, t39, t38, t71 * t13 - t73 * t14, 0, 0, 0, 0, 0, 0, t8, t105, t15, t71 * t11 - t73 * t2, 0, 0, 0, 0, 0, 0, t8, t15, -t105, -t73 * t1 + t71 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t23, -t24, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, -0.2e1 * t65 + t84, t55 + 0.2e1 * t58 + t24, -pkin(2) * t14 + qJ(3) * t13, t21, -t81, t22, t17, -t80, 0, t70 * t11 + t87, qJ(3) * t36 + t72 * t11 + t106, -t2 + t86, qJ(3) * t11 - t107 * t2, t21, t22, t81, 0, t80, t17, -t33 * t97 - t70 * (-qJ(5) * t90 + (-t33 - t89) * pkin(4) + t76) + t87, t72 * (qJ(5) * t41 + t5) - t70 * ((t41 - t74) * pkin(4) + t77) + t86, t72 * (-t83 * t92 + t76) - t79 * t36 - t106, -t1 * t107 + t79 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t75, t14, 0, 0, 0, 0, 0, 0, t18, -t19, -t37, t2, 0, 0, 0, 0, 0, 0, t18, -t37, t19, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t40, t93, -t49, -t94, qJDD(4), -t6, -t7, 0, 0, t49, t93, -t40, qJDD(4), t94, -t49, pkin(4) * t45 + qJ(5) * t47 - t5, -t83 * qJDD(2), qJ(5) * t44 + (t48 - t74) * pkin(4) + t77, -pkin(4) * t5 + qJ(5) * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t93, -t48, t5;];
tauJ_reg = t9;
