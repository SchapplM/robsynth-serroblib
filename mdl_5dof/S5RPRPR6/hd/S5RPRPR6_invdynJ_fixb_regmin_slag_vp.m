% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:51
% EndTime: 2019-12-31 18:17:52
% DurationCPUTime: 0.41s
% Computational Cost: add. (582->106), mult. (916->128), div. (0->0), fcn. (502->12), ass. (0->82)
t49 = qJ(1) + pkin(8);
t45 = qJ(3) + t49;
t40 = sin(t45);
t41 = cos(t45);
t107 = g(1) * t41 + g(2) * t40;
t54 = cos(pkin(8));
t42 = pkin(1) * t54 + pkin(2);
t56 = sin(qJ(3));
t53 = sin(pkin(8));
t103 = pkin(1) * t53;
t86 = qJD(3) * t103;
t59 = cos(qJ(3));
t94 = qJD(3) * t59;
t74 = t42 * t94 - t56 * t86;
t13 = -qJD(4) - t74;
t72 = t103 * t59 + t42 * t56;
t17 = qJ(4) + t72;
t47 = qJDD(1) + qJDD(3);
t48 = qJD(1) + qJD(3);
t109 = -t13 * t48 + t17 * t47 - t107;
t108 = -qJD(1) * t86 + t42 * qJDD(1);
t106 = g(1) * t40 - g(2) * t41;
t29 = t42 * qJD(1);
t100 = t29 * t56;
t87 = qJD(1) * t103;
t12 = t59 * t87 + t100;
t95 = qJ(4) * t48;
t10 = t12 + t95;
t105 = -t10 * t48 - t106;
t46 = t48 ^ 2;
t61 = -pkin(3) - pkin(7);
t43 = t47 * qJ(4);
t44 = t48 * qJD(4);
t85 = qJDD(1) * t103;
t68 = -t108 * t56 - t29 * t94 - t59 * t85;
t4 = -t43 - t44 + t68;
t55 = sin(qJ(5));
t58 = cos(qJ(5));
t93 = qJD(5) * t58;
t104 = t10 * t93 - t4 * t55;
t102 = pkin(3) * t47;
t99 = t58 * t47;
t98 = t41 * pkin(3) + t40 * qJ(4);
t62 = qJD(5) ^ 2;
t97 = -t46 - t62;
t51 = t58 ^ 2;
t96 = t55 ^ 2 - t51;
t11 = -t29 * t59 + t56 * t87;
t92 = qJD(4) + t11;
t52 = qJDD(2) - g(3);
t82 = -t56 * t103 + t42 * t59;
t18 = -pkin(3) - t82;
t15 = -pkin(7) + t18;
t91 = qJDD(5) * t15;
t90 = qJDD(5) * t55;
t89 = qJDD(5) * t58;
t88 = qJDD(5) * t61;
t84 = -pkin(3) * t40 + t41 * qJ(4);
t14 = t72 * qJD(3);
t83 = t17 * t48 + t14;
t81 = -t12 + t95;
t80 = -qJD(3) * t100 + t108 * t59 - t56 * t85;
t57 = sin(qJ(1));
t60 = cos(qJ(1));
t77 = g(1) * t57 - g(2) * t60;
t75 = qJDD(4) - t80;
t73 = t80 + t106;
t5 = t75 - t102;
t71 = -t47 * t61 - t105 - t75;
t70 = t12 * t48 + t73;
t69 = -t14 * t48 + t73;
t67 = t68 + t107;
t66 = -t15 * t62 + t109;
t65 = -t11 * t48 + t67;
t64 = t48 * t92 - t61 * t62 - t107 + t43;
t26 = -t55 * t62 + t89;
t25 = -t58 * t62 - t90;
t16 = -0.2e1 * t48 * t55 * t93 + t47 * t51;
t9 = -pkin(3) * t48 + t92;
t8 = 0.2e1 * qJD(5) * t48 * t96 - 0.2e1 * t55 * t99;
t2 = t4 * t58;
t1 = [qJDD(1), t77, g(1) * t60 + g(2) * t57, (t77 + (t53 ^ 2 + t54 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t47, t47 * t82 + t69, -t47 * t72 - t48 * t74 + t67, qJDD(4) + (-pkin(3) + t18) * t47 - t69, -t4 + t109, -t4 * t17 - t10 * t13 + t5 * t18 + t9 * t14 - g(1) * (-pkin(2) * sin(t49) - t57 * pkin(1) + t84) - g(2) * (pkin(2) * cos(t49) + t60 * pkin(1) + t98), t16, t8, t26, t25, 0, (qJD(5) * t83 + t91) * t58 + t66 * t55 + t104, -t2 + (-t91 + (-t10 - t83) * qJD(5)) * t55 + t66 * t58; 0, 0, 0, t52, 0, 0, 0, 0, 0, t52, 0, 0, 0, 0, 0, t25, -t26; 0, 0, 0, 0, t47, t70, t65, qJDD(4) - t70 - 0.2e1 * t102, 0.2e1 * t43 + 0.2e1 * t44 - t65, -t5 * pkin(3) - g(1) * t84 - g(2) * t98 - t4 * qJ(4) + t10 * t92 - t9 * t12, t16, t8, t26, t25, 0, (qJD(5) * t81 + t88) * t58 + t64 * t55 + t104, -t2 + (-t88 + (-t10 - t81) * qJD(5)) * t55 + t64 * t58; 0, 0, 0, 0, 0, 0, 0, t47, -t46, t5 + t105, 0, 0, 0, 0, 0, t55 * t97 + t89, t58 * t97 - t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t46 * t55, -t96 * t46, t99, -t55 * t47, qJDD(5), -t52 * t55 - t58 * t71, -t52 * t58 + t55 * t71;];
tau_reg = t1;
