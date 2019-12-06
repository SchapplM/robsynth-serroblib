% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PPRRP3
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PPRRP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:11:16
% EndTime: 2019-12-05 15:11:18
% DurationCPUTime: 0.54s
% Computational Cost: add. (775->113), mult. (1465->150), div. (0->0), fcn. (942->8), ass. (0->77)
t68 = sin(qJ(4));
t60 = t68 ^ 2;
t73 = qJD(3) ^ 2;
t105 = t60 * t73;
t72 = qJD(4) ^ 2;
t49 = t72 + t105;
t70 = cos(qJ(4));
t53 = t68 * t73 * t70;
t48 = qJDD(4) - t53;
t99 = t70 * t48;
t24 = -t68 * t49 + t99;
t87 = qJD(3) * qJD(4);
t89 = t68 * qJDD(3);
t38 = 0.2e1 * t70 * t87 + t89;
t65 = sin(pkin(8));
t66 = cos(pkin(8));
t69 = sin(qJ(3));
t71 = cos(qJ(3));
t112 = (t71 * t24 - t69 * t38) * t65 - t66 * (t68 * t48 + t70 * t49);
t111 = pkin(6) * t24;
t110 = t69 * t24 + t71 * t38;
t61 = t70 ^ 2;
t104 = t61 * t73;
t107 = t99 + t68 * (-t72 + t104);
t92 = sin(pkin(7));
t93 = cos(pkin(7));
t46 = -t93 * g(1) - t92 * g(2);
t62 = -g(3) + qJDD(1);
t28 = t66 * t46 + t65 * t62;
t74 = -t92 * g(1) + t93 * g(2) + qJDD(2);
t15 = t71 * t28 + t69 * t74;
t12 = -t73 * pkin(3) + qJDD(3) * pkin(6) + t15;
t79 = t65 * t46 - t66 * t62;
t27 = t70 * t79;
t94 = t68 * qJ(5);
t81 = -t70 * pkin(4) - t94;
t91 = t73 * t81;
t4 = -qJDD(4) * pkin(4) - t72 * qJ(5) + (t12 + t91) * t68 + qJDD(5) - t27;
t106 = 2 * qJD(5);
t47 = qJDD(4) + t53;
t103 = t68 * t47;
t9 = t70 * t12 + t68 * t79;
t51 = -t72 - t104;
t23 = t70 * t51 - t103;
t86 = t68 * t87;
t88 = t70 * qJDD(3);
t39 = -0.2e1 * t86 + t88;
t97 = pkin(3) * t39 + pkin(6) * t23;
t95 = t60 + t61;
t41 = t95 * qJDD(3);
t44 = t95 * t73;
t96 = pkin(3) * t44 + pkin(6) * t41;
t90 = qJD(3) * t68;
t8 = t68 * t12 - t27;
t2 = t68 * t8 + t70 * t9;
t14 = -t69 * t28 + t71 * t74;
t82 = qJDD(4) * qJ(5) + (qJD(4) * t106) + t70 * t91 + t9;
t80 = t70 * t38 + t68 * t39;
t11 = -qJDD(3) * pkin(3) - t73 * pkin(6) - t14;
t77 = pkin(3) - t81;
t76 = t11 - (-t86 + t88) * pkin(4) - qJ(5) * t38;
t75 = t90 * t106 - t76;
t45 = (t60 - t61) * t73;
t43 = -t69 * qJDD(3) - t71 * t73;
t42 = t71 * qJDD(3) - t69 * t73;
t22 = t103 + t70 * (t72 - t105);
t21 = t38 * t68;
t20 = t39 * t70;
t19 = t66 * t79;
t17 = t69 * t41 + t71 * t44;
t16 = t65 * (t71 * t41 - t69 * t44);
t13 = t69 * t23 + t71 * t39;
t7 = t65 * (t71 * t23 - t69 * t39) + t66 * (-t70 * t47 - t68 * t51);
t5 = (pkin(4) * qJD(4) - (2 * qJD(5))) * t90 + t76;
t3 = -t72 * pkin(4) + t82;
t1 = t70 * t3 + t68 * t4;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65 * t28 - t19, 0, 0, 0, 0, 0, 0, t65 * t43, -t65 * t42, 0, t65 * (-t69 * t14 + t71 * t15) - t19, 0, 0, 0, 0, 0, 0, t7, -t112, t16, t65 * (t69 * t11 + t71 * t2) + t66 * (-t68 * t9 + t70 * t8), 0, 0, 0, 0, 0, 0, t7, t16, t112, t65 * (t71 * t1 + t69 * t5) + t66 * (-t68 * t3 + t70 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, 0, 0, 0, 0, 0, t42, t43, 0, t71 * t14 + t69 * t15, 0, 0, 0, 0, 0, 0, t13, -t110, t17, -t71 * t11 + t69 * t2, 0, 0, 0, 0, 0, 0, t13, t17, t110, t69 * t1 - t71 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t14, -t15, 0, 0, t21, t80, t22, t20, t107, 0, -t70 * t11 + t97, -pkin(3) * t38 + t68 * t11 - t111, t2 + t96, -pkin(3) * t11 + pkin(6) * t2, t21, t22, -t80, 0, -t107, t20, t39 * t94 + t70 * ((t39 - t86) * pkin(4) + t75) + t97, t70 * ((t44 - t72) * pkin(4) + t82) + (qJ(5) * t44 + t4) * t68 + t96, t68 * (-pkin(4) * t86 + t75) + t111 + t77 * t38, pkin(6) * t1 - t77 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t45, t89, t53, t88, qJDD(4), -t8, -t9, 0, 0, -t53, t89, -t45, qJDD(4), -t88, t53, pkin(4) * t47 + qJ(5) * t51 - t4, (-pkin(4) * t68 + qJ(5) * t70) * qJDD(3), qJ(5) * t48 + (t49 - t72) * pkin(4) + t82, -pkin(4) * t4 + qJ(5) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t89, -t49, t4;];
tauJ_reg = t6;
