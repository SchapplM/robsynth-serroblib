% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PPRRP1
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PPRRP1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:07:14
% EndTime: 2019-12-05 15:07:16
% DurationCPUTime: 0.47s
% Computational Cost: add. (804->114), mult. (1551->156), div. (0->0), fcn. (1004->8), ass. (0->79)
t86 = (qJD(3) * qJD(5));
t102 = 2 * t86;
t72 = sin(qJ(4));
t74 = cos(qJ(4));
t77 = qJD(3) ^ 2;
t55 = t72 * t77 * t74;
t49 = qJDD(4) + t55;
t101 = pkin(4) * t49;
t64 = -g(3) + qJDD(1);
t67 = sin(pkin(8));
t69 = cos(pkin(8));
t68 = sin(pkin(7));
t70 = cos(pkin(7));
t82 = -g(1) * t70 - g(2) * t68;
t29 = t64 * t67 + t69 * t82;
t73 = sin(qJ(3));
t75 = cos(qJ(3));
t78 = t69 * t64 - t67 * t82;
t16 = t29 * t75 + t73 * t78;
t14 = -pkin(3) * t77 + qJDD(3) * pkin(6) + t16;
t42 = -g(1) * t68 + g(2) * t70 + qJDD(2);
t10 = t74 * t14 + t72 * t42;
t61 = t74 * qJDD(3);
t87 = qJD(3) * qJD(4);
t84 = t72 * t87;
t39 = t61 - t84;
t90 = qJD(3) * t72;
t48 = qJD(4) * pkin(4) - qJ(5) * t90;
t79 = t39 * qJ(5) - qJD(4) * t48 + t102 * t74 + t10;
t62 = t72 ^ 2;
t100 = t62 * t77;
t63 = t74 ^ 2;
t99 = t63 * t77;
t98 = t72 * t14;
t97 = t72 * t49;
t50 = qJDD(4) - t55;
t96 = t74 * t50;
t76 = qJD(4) ^ 2;
t53 = -t76 - t99;
t26 = t53 * t74 - t97;
t40 = t61 - 0.2e1 * t84;
t95 = pkin(3) * t40 + pkin(6) * t26;
t52 = -t76 - t100;
t27 = -t52 * t72 - t96;
t83 = t74 * t87;
t89 = t72 * qJDD(3);
t37 = 0.2e1 * t83 + t89;
t94 = -pkin(3) * t37 + pkin(6) * t27;
t92 = t62 + t63;
t43 = t92 * qJDD(3);
t46 = t92 * t77;
t93 = pkin(3) * t46 + pkin(6) * t43;
t91 = qJ(5) * t72;
t88 = qJ(5) * qJDD(3);
t31 = t74 * t42;
t9 = -t31 + t98;
t2 = t10 * t74 + t72 * t9;
t15 = -t29 * t73 + t75 * t78;
t13 = -qJDD(3) * pkin(3) - t77 * pkin(6) - t15;
t38 = t83 + t89;
t80 = -t31 + (t38 - t83) * qJ(5) - t101;
t3 = -0.2e1 * t72 * t86 - t80 - t98;
t7 = -t39 * pkin(4) - qJ(5) * t99 + t48 * t90 + qJDD(5) + t13;
t47 = (t62 - t63) * t77;
t45 = -qJDD(3) * t73 - t75 * t77;
t44 = qJDD(3) * t75 - t73 * t77;
t25 = -t50 * t72 + t52 * t74;
t24 = t97 + t74 * (t76 - t100);
t23 = t49 * t74 + t53 * t72;
t22 = t72 * (-t76 + t99) + t96;
t21 = (t38 + t83) * t72;
t20 = (t39 - t84) * t74;
t17 = t37 * t74 + t40 * t72;
t11 = t67 * (t43 * t75 - t46 * t73) + t69 * (t43 * t73 + t46 * t75);
t6 = t67 * (t27 * t75 + t37 * t73) + t69 * (t27 * t73 - t37 * t75);
t5 = t67 * (t26 * t75 - t40 * t73) + t69 * (t26 * t73 + t40 * t75);
t4 = -pkin(4) * t99 + t79;
t1 = -t3 * t72 + t4 * t74;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67 * t29 + t69 * t78, 0, 0, 0, 0, 0, 0, t44 * t69 + t45 * t67, -t44 * t67 + t45 * t69, 0, t67 * (-t15 * t73 + t16 * t75) + t69 * (t15 * t75 + t16 * t73), 0, 0, 0, 0, 0, 0, t5, t6, t11, t67 * (t13 * t73 + t2 * t75) + t69 * (-t13 * t75 + t2 * t73), 0, 0, 0, 0, 0, 0, t5, t6, t11, t67 * (t1 * t75 + t7 * t73) + t69 * (t1 * t73 - t7 * t75); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, t23, t25, 0, t10 * t72 - t74 * t9, 0, 0, 0, 0, 0, 0, t23, t25, 0, t3 * t74 + t4 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t15, -t16, 0, 0, t21, t17, t24, t20, t22, 0, -t13 * t74 + t95, t13 * t72 + t94, t2 + t93, -pkin(3) * t13 + pkin(6) * t2, t21, t17, t24, t20, t22, 0, -t49 * t91 + t74 * (pkin(4) * t40 + qJ(5) * t53 - t7) + t95, t72 * (-qJ(5) * t52 + t7) + t74 * (-pkin(4) * t37 - qJ(5) * t50) + t94, t74 * (t74 * t88 + (t46 - t99) * pkin(4) + t79) + ((t14 + t102 + t88) * t72 + t80) * t72 + t93, -t3 * t91 + t74 * (-pkin(4) * t7 + qJ(5) * t4) - pkin(3) * t7 + pkin(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t47, t89, t55, t61, qJDD(4), -t9, -t10, 0, 0, -t55, t47, t89, t55, t61, qJDD(4), t3 + t101, (t52 + t99) * pkin(4) - t79, -pkin(4) * t89, pkin(4) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, t37, -t46, t7;];
tauJ_reg = t8;
