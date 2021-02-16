% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PPRRP4
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% tau_reg [5x16]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:56
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:56:29
% EndTime: 2021-01-15 14:56:32
% DurationCPUTime: 0.49s
% Computational Cost: add. (474->150), mult. (997->183), div. (0->0), fcn. (636->6), ass. (0->94)
t50 = cos(qJ(4));
t49 = sin(qJ(3));
t28 = qJD(3) * pkin(6) + t49 * qJD(2);
t85 = qJ(5) * qJD(3);
t68 = t28 + t85;
t112 = t50 * t68;
t51 = cos(qJ(3));
t53 = qJD(3) ^ 2;
t111 = -qJDD(3) * t49 - t53 * t51;
t48 = sin(qJ(4));
t75 = qJD(3) * qJD(4);
t72 = t48 * t75;
t30 = pkin(4) * t72;
t76 = qJD(2) * qJD(3);
t34 = t49 * t76;
t79 = t51 * qJDD(2);
t67 = t34 - t79;
t37 = t50 * pkin(4) + pkin(3);
t84 = qJDD(3) * t37;
t5 = qJDD(5) + t30 + t67 - t84;
t7 = -t50 * qJD(1) - t48 * t68;
t94 = qJD(4) * pkin(4);
t6 = t7 + t94;
t40 = t48 * qJD(1);
t8 = -t40 + t112;
t62 = t48 * t6 - t50 * t8;
t110 = -qJD(3) * t62 - t5;
t109 = t6 - t7;
t44 = t48 ^ 2;
t108 = pkin(4) * t44;
t92 = sin(pkin(7));
t93 = cos(pkin(7));
t18 = -t49 * t92 - t51 * t93;
t107 = g(1) * t18;
t19 = t49 * t93 - t51 * t92;
t106 = g(1) * t19;
t105 = g(2) * t18;
t104 = g(2) * t19;
t43 = g(3) * t50;
t103 = t19 * t48;
t102 = t50 * t53;
t47 = qJ(5) + pkin(6);
t86 = t51 * qJD(2);
t73 = qJD(4) * t86;
t100 = t50 * t34 + t48 * t73;
t99 = -g(1) * t103 + t50 * t73;
t45 = t50 ^ 2;
t98 = t44 - t45;
t97 = t44 + t45;
t52 = qJD(4) ^ 2;
t96 = t52 + t53;
t95 = qJD(3) * pkin(3);
t90 = qJD(3) * t37;
t15 = qJD(5) - t86 - t90;
t91 = qJD(3) * t15;
t88 = qJDD(4) * pkin(4);
t87 = t28 * qJD(4);
t82 = qJDD(4) * t48;
t38 = t48 * qJDD(3);
t81 = t49 * qJDD(2);
t80 = t50 * qJDD(3);
t78 = t51 * qJDD(3);
t77 = qJD(1) * qJD(4);
t74 = -t50 * t77 + (-qJDD(1) - t87) * t48;
t71 = qJD(4) * t47;
t70 = -t50 * t107 - t74;
t17 = qJDD(3) * pkin(6) + t51 * t76 + t81;
t29 = -t86 - t95;
t69 = -t29 * qJD(3) - t17;
t66 = 0.2e1 * t50 * t75;
t65 = -qJ(5) * qJDD(3) - t17;
t64 = -t105 + t106;
t63 = t104 + t107;
t61 = -g(1) * t92 + g(2) * t93;
t60 = t5 - t84 + t105;
t35 = t48 * t77;
t59 = -g(2) * t103 - t50 * qJDD(1) - t48 * t107 + t35 + t43;
t58 = 0.2e1 * t72 - t80;
t57 = qJD(3) * qJD(5) - t65;
t56 = 0.2e1 * qJDD(3) * pkin(3) - pkin(6) * t52 - t105 - t67;
t55 = (-qJD(5) - t15) * qJD(3) + t65;
t54 = -pkin(6) * qJDD(4) + (t29 - t95) * qJD(4);
t46 = qJDD(1) - g(3);
t24 = t47 * t50;
t23 = t47 * t48;
t22 = qJDD(4) * t50 - t52 * t48;
t21 = t52 * t50 + t82;
t14 = -t48 * qJD(5) - t50 * t71;
t13 = t50 * qJD(5) - t48 * t71;
t4 = -t58 * t51 + (-t50 * t96 - t82) * t49;
t3 = (-qJDD(4) * t49 - 0.2e1 * t51 * t75) * t50 + (t49 * t96 - t78) * t48;
t2 = -qJ(5) * t72 + t50 * t57 + t74;
t1 = t88 + t35 + (-qJD(4) * t68 - qJDD(1)) * t50 - t57 * t48;
t9 = [t46, t46, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t21, -t22, t21, 0, qJD(4) * t62 - t1 * t50 - t2 * t48 - g(3); 0, qJDD(2) + t61, 0, -t53 * t49 + t78, t111, 0, 0, 0, 0, 0, t4, t3, t4, t3, -t111 * t97, t110 * t51 + (t91 - t1 * t48 + t2 * t50 + (-t48 * t8 - t50 * t6) * qJD(4)) * t49 + t61; 0, 0, qJDD(3), t64 + t79, -t63 - t81, t44 * qJDD(3) + t48 * t66, 0.2e1 * t48 * t80 - 0.2e1 * t75 * t98, t21, t22, 0, t54 * t48 + (t56 + t106) * t50 + t100, t54 * t50 + (-t56 - t34) * t48 + t99, -t23 * qJDD(4) + (t14 + (t15 - t90) * t48) * qJD(4) + (-t30 - t60 + t106) * t50 + t100, -t24 * qJDD(4) + (t60 - t34) * t48 + (t15 * t50 - t13 + (-t37 * t50 + t108) * qJD(3)) * qJD(4) + t99, (-qJD(4) * t6 + qJDD(3) * t24 + t2) * t50 + (-t8 * qJD(4) + qJDD(3) * t23 - t1) * t48 + (t13 * t50 - t14 * t48 + (t23 * t50 - t24 * t48) * qJD(4) - t97 * t86) * qJD(3) + t63, t2 * t24 + t8 * t13 - t1 * t23 + t6 * t14 - t5 * t37 + t15 * t48 * t94 - g(1) * (-t18 * t47 - t19 * t37) - g(2) * (t18 * t37 - t19 * t47) + (-t15 * t49 + t51 * t62) * qJD(2); 0, 0, 0, 0, 0, -t48 * t102, t98 * t53, t38, t80, qJDD(4), -t40 * qJD(4) + t48 * t69 + t59, (-g(3) - t87) * t48 + (t69 - t77 - t104) * t50 + t70, 0.2e1 * t88 + (t8 - t112) * qJD(4) + (pkin(4) * t102 + t55) * t48 + t59, -t53 * t108 - g(3) * t48 + (t48 * t85 + t7) * qJD(4) + (t55 - t104) * t50 + t70, -pkin(4) * t38 + (-t94 + t109) * t50 * qJD(3), t109 * t8 + (t43 + t1 + (-t63 - t91) * t48) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t38 + t66, -t97 * t53, -t110 - t64;];
tau_reg = t9;
