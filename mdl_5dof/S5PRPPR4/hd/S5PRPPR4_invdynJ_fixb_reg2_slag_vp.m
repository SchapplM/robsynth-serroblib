% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRPPR4
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPPR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:37:00
% EndTime: 2019-12-31 17:37:01
% DurationCPUTime: 0.79s
% Computational Cost: add. (810->175), mult. (1620->221), div. (0->0), fcn. (1171->6), ass. (0->104)
t140 = qJ(3) * qJDD(2) + qJD(2) * qJD(3);
t87 = sin(pkin(8));
t83 = t87 ^ 2;
t88 = cos(pkin(8));
t84 = t88 ^ 2;
t141 = t83 + t84;
t112 = qJD(2) * qJD(4);
t106 = t87 * t112;
t122 = t87 * qJ(4);
t105 = pkin(2) + t122;
t40 = (pkin(3) + pkin(4)) * t88 + t105;
t13 = t40 * qJDD(2) - qJDD(3) + t106;
t85 = pkin(7) + qJ(2);
t80 = sin(t85);
t81 = cos(t85);
t123 = g(1) * t81 + g(2) * t80;
t89 = sin(qJ(5));
t90 = cos(qJ(5));
t42 = t87 * t89 + t88 * t90;
t37 = t42 * qJD(2);
t43 = t87 * t90 - t88 * t89;
t120 = qJD(2) * t88;
t111 = t89 * t120;
t94 = qJD(5) * t111 - t42 * qJDD(2);
t138 = t37 ^ 2;
t121 = qJD(2) * t87;
t110 = t90 * t121;
t39 = t110 - t111;
t137 = t39 ^ 2;
t135 = g(1) * t80;
t69 = g(2) * t81;
t35 = t42 * qJD(5);
t9 = qJD(5) * t110 - t94;
t134 = t35 * t37 - t43 * t9;
t31 = t87 * qJDD(1) + t140 * t88;
t133 = t31 * t87 - g(3);
t27 = t31 * t88;
t132 = t39 * t37;
t117 = qJ(3) * qJD(2);
t45 = t87 * qJD(1) + t88 * t117;
t131 = t45 * t88;
t130 = t80 * t88;
t129 = t81 * t88;
t126 = -pkin(6) + qJ(3);
t125 = qJ(3) * t27 + qJD(3) * t131;
t124 = t81 * pkin(2) + t80 * qJ(3);
t119 = qJD(4) * t87;
t118 = qJDD(2) * pkin(2);
t46 = -t88 * pkin(3) - t105;
t116 = qJDD(2) * t46;
t76 = t87 * qJDD(2);
t115 = t88 * qJDD(2);
t109 = t87 * t115;
t58 = t81 * qJ(3);
t108 = -t80 * pkin(2) + t58;
t107 = -t69 + t135;
t30 = t88 * qJDD(1) - t140 * t87;
t28 = qJDD(4) - t30;
t17 = -pkin(6) * t76 + t28;
t18 = -pkin(6) * t115 + t31;
t104 = t90 * t17 - t89 * t18;
t103 = pkin(3) * t129 + t81 * t122 + t124;
t44 = t88 * qJD(1) - t87 * t117;
t74 = qJDD(3) - t118;
t101 = t89 * t115 - t90 * t76;
t36 = t43 * qJD(5);
t8 = qJD(2) * t35 + t101;
t100 = -t39 * t36 + t42 * t8;
t99 = t89 * t17 + t90 * t18;
t41 = qJD(4) - t44;
t29 = -pkin(6) * t121 + t41;
t32 = -pkin(6) * t120 + t45;
t6 = t90 * t29 - t89 * t32;
t7 = t89 * t29 + t90 * t32;
t48 = t126 * t87;
t49 = t126 * t88;
t15 = t90 * t48 - t89 * t49;
t16 = t89 * t48 + t90 * t49;
t97 = t118 - t74 - t69;
t93 = qJDD(3) + t116;
t19 = t93 - t106;
t96 = -t116 - t19 - t69;
t95 = t140 * t141 - t123 + t27;
t92 = qJD(2) ^ 2;
t91 = qJD(5) ^ 2;
t86 = qJDD(1) - g(3);
t82 = g(3) * t88;
t78 = t84 * qJDD(2);
t75 = t83 * qJDD(2);
t52 = g(1) * t130;
t47 = t141 * t92;
t33 = t46 * qJD(2) + qJD(3);
t25 = t42 * t81;
t24 = t43 * t81;
t23 = t42 * t80;
t22 = t43 * t80;
t21 = t40 * qJD(2) - qJD(3);
t11 = -t36 * qJD(5) - t42 * qJDD(5);
t10 = -t35 * qJD(5) + t43 * qJDD(5);
t5 = t43 * qJD(3) - t16 * qJD(5);
t4 = t42 * qJD(3) + t15 * qJD(5);
t2 = -t7 * qJD(5) + t104;
t1 = t6 * qJD(5) + t99;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t86, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 * t88 + t133, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28 * t88 + t133, 0, 0, 0, 0, 0, 0, t11, -t10, -t100 + t134, t1 * t43 - t2 * t42 - t7 * t35 - t6 * t36 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t107, t123, 0, 0, t75, 0.2e1 * t109, 0, t78, 0, 0, t97 * t88 + t52, (-t97 - t135) * t87, -t30 * t87 + t95, -t74 * pkin(2) - g(1) * t108 - g(2) * t124 + (-t30 * qJ(3) - t44 * qJD(3)) * t87 + t125, t75, 0, -0.2e1 * t109, 0, 0, t78, t52 + (t96 + t106) * t88, t28 * t87 + t95, t83 * t112 + (t96 + t135) * t87, t19 * t46 - g(1) * (-pkin(3) * t130 + t108) - g(2) * t103 + (t28 * qJ(3) + qJ(4) * t135 + t41 * qJD(3) - t33 * qJD(4)) * t87 + t125, -t39 * t35 - t8 * t43, t100 + t134, t10, t37 * t36 + t9 * t42, t11, 0, g(1) * t23 - g(2) * t25 + t5 * qJD(5) + t15 * qJDD(5) + t37 * t119 + t13 * t42 + t21 * t36 + t40 * t9, g(1) * t22 - g(2) * t24 - t4 * qJD(5) - t16 * qJDD(5) + t39 * t119 + t13 * t43 - t21 * t35 - t40 * t8, -t1 * t42 + t15 * t8 - t16 * t9 - t2 * t43 + t6 * t35 - t7 * t36 - t4 * t37 - t5 * t39 + t123, t1 * t16 + t7 * t4 + t2 * t15 + t6 * t5 + t13 * t40 + t21 * t119 - g(1) * (-t81 * pkin(6) + t58) - g(2) * (pkin(4) * t129 + t103) + (-g(1) * (-t88 * pkin(4) + t46) + g(2) * pkin(6)) * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115, t76, -t47, (t44 * t87 - t131) * qJD(2) + t74 - t107, 0, 0, 0, 0, 0, 0, -t115, -t47, -t76, (-t131 + (-qJD(4) - t41) * t87) * qJD(2) + t93 - t107, 0, 0, 0, 0, 0, 0, (-t39 - t110) * qJD(5) + t94, 0.2e1 * t37 * qJD(5) + t101, t137 + t138, -t7 * t37 - t6 * t39 - t107 - t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87 * t92 * t88, t76, -t83 * t92, t82 + (qJD(2) * t33 - t123) * t87 + t28, 0, 0, 0, 0, 0, 0, qJDD(5) * t90 - t37 * t121 - t91 * t89, -qJDD(5) * t89 - t39 * t121 - t91 * t90, t90 * t8 - t89 * t9 + (-t37 * t90 + t39 * t89) * qJD(5), t1 * t89 + t2 * t90 + t82 + (-t6 * t89 + t7 * t90) * qJD(5) + (-qJD(2) * t21 - t123) * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, t137 - t138, -t101, -t132, (t39 - t110) * qJD(5) + t94, qJDD(5), -g(1) * t24 - g(2) * t22 + g(3) * t42 - t21 * t39 + t104, g(1) * t25 + g(2) * t23 + g(3) * t43 + t21 * t37 - t99, 0, 0;];
tau_reg = t3;
