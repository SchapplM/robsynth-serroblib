% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PPRPR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PPRPR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:01:20
% EndTime: 2019-12-05 15:01:25
% DurationCPUTime: 0.81s
% Computational Cost: add. (1606->145), mult. (3238->235), div. (0->0), fcn. (2405->10), ass. (0->90)
t108 = qJD(3) ^ 2;
t73 = sin(pkin(8));
t75 = cos(pkin(8));
t100 = sin(pkin(7));
t101 = cos(pkin(7));
t87 = -t101 * g(1) - t100 * g(2);
t96 = -g(3) + qJDD(1);
t40 = t73 * t96 + t75 * t87;
t77 = sin(qJ(3));
t79 = cos(qJ(3));
t85 = -t73 * t87 + t75 * t96;
t29 = t79 * t40 + t77 * t85;
t115 = -t108 * pkin(3) + qJDD(3) * qJ(4) + (2 * qJD(3) * qJD(4)) + t29;
t72 = sin(pkin(9));
t74 = cos(pkin(9));
t76 = sin(qJ(5));
t78 = cos(qJ(5));
t49 = (t72 * t76 - t74 * t78) * qJD(3);
t88 = t72 * t78 + t74 * t76;
t51 = t88 * qJD(3);
t38 = t51 * t49;
t109 = qJDD(5) - t38;
t114 = t109 * t76;
t113 = t109 * t78;
t81 = t72 ^ 2;
t83 = t74 ^ 2;
t112 = t81 + t83;
t111 = t108 * t74;
t56 = -t100 * g(1) + t101 * g(2) + qJDD(2);
t53 = t74 * t56;
t110 = t53 + (pkin(4) * t111 - pkin(6) * qJDD(3) - t115) * t72;
t58 = t112 * t108;
t46 = t49 ^ 2;
t47 = t51 ^ 2;
t15 = t115 * t74 + t72 * t56;
t94 = t74 * qJDD(3);
t97 = t83 * t108;
t12 = -pkin(4) * t97 + pkin(6) * t94 + t15;
t5 = -t78 * t110 + t76 * t12;
t6 = t110 * t76 + t78 * t12;
t2 = -t78 * t5 + t76 * t6;
t107 = t72 * t2;
t28 = -t77 * t40 + t79 * t85;
t71 = qJDD(3) * pkin(3);
t24 = -t108 * qJ(4) + qJDD(4) - t28 - t71;
t20 = -pkin(4) * t94 + t24 + (-t108 * t81 - t97) * pkin(6);
t106 = t76 * t20;
t32 = qJDD(5) + t38;
t105 = t76 * t32;
t104 = t78 * t20;
t103 = t78 * t32;
t99 = t49 * qJD(5);
t98 = t51 * qJD(5);
t95 = t72 * qJDD(3);
t93 = t77 * qJDD(3);
t92 = t79 * qJDD(3);
t3 = t76 * t5 + t78 * t6;
t14 = t115 * t72 - t53;
t7 = t72 * t14 + t74 * t15;
t89 = -t24 + t71;
t26 = -t76 * t95 + t78 * t94;
t48 = t88 * qJDD(3);
t80 = qJD(5) ^ 2;
t69 = t83 * qJDD(3);
t68 = t81 * qJDD(3);
t60 = -t79 * t108 - t93;
t59 = -t77 * t108 + t92;
t57 = t69 + t68;
t55 = t112 * t111;
t54 = t72 * t58;
t43 = -t47 - t80;
t42 = -t47 + t80;
t41 = t46 - t80;
t37 = t48 - t99;
t36 = t48 - 0.2e1 * t99;
t35 = t26 - t98;
t34 = -t26 + 0.2e1 * t98;
t30 = -t80 - t46;
t27 = -t46 - t47;
t22 = -t76 * t43 - t103;
t21 = t78 * t43 - t105;
t19 = t78 * t26 + t76 * t48;
t18 = t76 * t26 - t78 * t48;
t17 = t78 * t30 - t114;
t16 = t76 * t30 + t113;
t10 = -t72 * t21 + t74 * t22;
t9 = -t72 * t18 + t74 * t19;
t8 = -t72 * t16 + t74 * t17;
t1 = t74 * t3 - t107;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t96, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t40 + t75 * t85, 0, 0, 0, 0, 0, 0, t75 * t59 + t73 * t60, -t73 * t59 + t75 * t60, 0, t73 * (-t77 * t28 + t79 * t29) + t75 * (t79 * t28 + t77 * t29), 0, 0, 0, 0, 0, 0, t73 * (-t79 * t55 - t74 * t93) + t75 * (-t77 * t55 + t74 * t92), t73 * (t79 * t54 + t72 * t93) + t75 * (t77 * t54 - t72 * t92), t73 * (t79 * t57 - t77 * t58) + t75 * (t77 * t57 + t79 * t58), t73 * (t77 * t24 + t79 * t7) + t75 * (-t79 * t24 + t77 * t7), 0, 0, 0, 0, 0, 0, t73 * (t77 * t34 + t79 * t8) + t75 * (-t79 * t34 + t77 * t8), t73 * (t79 * t10 + t77 * t36) + t75 * (t77 * t10 - t79 * t36), t73 * (t77 * t27 + t79 * t9) + t75 * (-t79 * t27 + t77 * t9), t73 * (t79 * t1 + t77 * t20) + t75 * (t77 * t1 - t79 * t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74 * t14 + t72 * t15, 0, 0, 0, 0, 0, 0, t74 * t16 + t72 * t17, t74 * t21 + t72 * t22, t74 * t18 + t72 * t19, t74 * t2 + t72 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t28, -t29, 0, 0, t68, 0.2e1 * t72 * t94, 0, t69, 0, 0, -qJ(4) * t55 + t89 * t74, qJ(4) * t54 - t89 * t72, pkin(3) * t58 + qJ(4) * t57 + t7, -pkin(3) * t24 + qJ(4) * t7, t72 * (t78 * t37 - t76 * t98) + t74 * (t76 * t37 + t78 * t98), t72 * (-t78 * t34 - t76 * t36) + t74 * (-t76 * t34 + t78 * t36), t72 * (-t76 * t42 + t113) + t74 * (t78 * t42 + t114), t72 * (-t76 * t35 + t78 * t99) + t74 * (t78 * t35 + t76 * t99), t72 * (t78 * t41 - t105) + t74 * (t76 * t41 + t103), (t72 * (-t49 * t78 + t51 * t76) + t74 * (-t49 * t76 - t51 * t78)) * qJD(5), t72 * (-pkin(6) * t16 + t106) + t74 * (-pkin(4) * t34 + pkin(6) * t17 - t104) - pkin(3) * t34 + qJ(4) * t8, t72 * (-pkin(6) * t21 + t104) + t74 * (-pkin(4) * t36 + pkin(6) * t22 + t106) - pkin(3) * t36 + qJ(4) * t10, t72 * (-pkin(6) * t18 - t2) + t74 * (-pkin(4) * t27 + pkin(6) * t19 + t3) - pkin(3) * t27 + qJ(4) * t9, -pkin(6) * t107 + t74 * (-pkin(4) * t20 + pkin(6) * t3) - pkin(3) * t20 + qJ(4) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, t95, -t58, t24, 0, 0, 0, 0, 0, 0, t34, t36, t27, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t47 - t46, t48, -t38, t26, qJDD(5), -t5, -t6, 0, 0;];
tauJ_reg = t4;
