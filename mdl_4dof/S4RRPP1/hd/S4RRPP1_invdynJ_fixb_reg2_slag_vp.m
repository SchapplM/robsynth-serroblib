% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRPP1
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
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:51:42
% EndTime: 2018-11-14 13:51:43
% DurationCPUTime: 0.29s
% Computational Cost: add. (395->102), mult. (605->118), div. (0->0), fcn. (331->10), ass. (0->74)
t61 = cos(qJ(2));
t83 = pkin(1) * qJD(2);
t76 = qJD(1) * t83;
t59 = sin(qJ(2));
t82 = qJDD(1) * t59;
t98 = pkin(1) * t82 + t61 * t76;
t56 = qJ(1) + qJ(2);
t50 = sin(t56);
t51 = cos(t56);
t97 = g(1) * t50 - g(2) * t51;
t53 = qJDD(1) + qJDD(2);
t46 = t53 * qJ(4);
t54 = qJD(1) + qJD(2);
t47 = t54 * qJD(4);
t96 = t46 + t47;
t95 = pkin(3) * t53 - qJDD(4);
t94 = pkin(2) * t50;
t57 = sin(pkin(6));
t92 = t57 * pkin(2);
t58 = cos(pkin(6));
t91 = t58 * pkin(2);
t60 = sin(qJ(1));
t90 = t60 * pkin(1);
t89 = t61 * pkin(1);
t88 = t57 * t59;
t87 = t58 * t59;
t86 = t58 * t61;
t44 = pkin(2) + t89;
t20 = pkin(1) * t87 + t44 * t57;
t85 = g(1) * t51 + g(2) * t50;
t84 = pkin(1) * qJD(1);
t81 = pkin(1) * t88;
t45 = qJDD(1) * t89;
t12 = pkin(2) * t53 - t59 * t76 + t45;
t4 = t57 * t12 + t58 * t98;
t80 = -t58 * t12 + t57 * t98;
t49 = pkin(6) + t56;
t37 = sin(t49);
t38 = cos(t49);
t43 = pkin(2) * t51;
t79 = pkin(3) * t38 + qJ(4) * t37 + t43;
t78 = t59 * t84;
t75 = qJD(1) * (-qJD(2) + t54);
t74 = qJD(2) * (-qJD(1) - t54);
t72 = t45 + t97;
t71 = g(1) * t38 + g(2) * t37 - t4;
t62 = cos(qJ(1));
t70 = g(1) * t60 - g(2) * t62;
t22 = pkin(2) * t54 + t61 * t84;
t8 = t22 * t57 + t58 * t78;
t18 = -qJD(2) * t81 + t83 * t86;
t69 = -pkin(3) * t37 + qJ(4) * t38 - t94;
t19 = t44 * t58 - t81;
t68 = g(1) * t37 - g(2) * t38 - t80;
t67 = pkin(1) * (t57 * t61 + t87);
t7 = t22 * t58 - t57 * t78;
t17 = (t86 - t88) * t84;
t66 = t17 * t54 + t71;
t15 = qJD(1) * t67;
t65 = t15 * t54 + t68;
t16 = qJD(2) * t67;
t64 = -t16 * t54 + t68;
t55 = qJDD(3) - g(3);
t52 = t62 * pkin(1);
t39 = -pkin(3) - t91;
t36 = qJ(4) + t92;
t14 = -pkin(3) - t19;
t13 = qJ(4) + t20;
t11 = qJD(4) + t18;
t6 = qJ(4) * t54 + t8;
t5 = -pkin(3) * t54 + qJD(4) - t7;
t2 = t80 - t95;
t1 = t4 + t96;
t3 = [0, 0, 0, 0, 0, qJDD(1), t70, g(1) * t62 + g(2) * t60, 0, 0, 0, 0, 0, 0, 0, t53 (t53 * t61 + t59 * t74) * pkin(1) + t72 ((-qJDD(1) - t53) * t59 + t61 * t74) * pkin(1) + t85, 0 (t70 + (t59 ^ 2 + t61 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t53, t19 * t53 + t64, -t18 * t54 - t20 * t53 + t71, 0, t4 * t20 + t8 * t18 - t80 * t19 - t7 * t16 - g(1) * (-t90 - t94) - g(2) * (t43 + t52) 0, 0, 0, t53, 0, 0, -t14 * t53 + t64 + t95, 0, t11 * t54 + t13 * t53 - t71 + t96, t1 * t13 + t6 * t11 + t2 * t14 + t5 * t16 - g(1) * (t69 - t90) - g(2) * (t52 + t79); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, pkin(1) * t59 * t75 + t72 (t61 * t75 - t82) * pkin(1) + t85, 0, 0, 0, 0, 0, 0, 0, t53, t53 * t91 + t65, -t53 * t92 + t66, 0, t7 * t15 - t8 * t17 + (t4 * t57 - t58 * t80 + t97) * pkin(2), 0, 0, 0, t53, 0, 0, -t39 * t53 + t65 + t95, 0, t36 * t53 + t46 + 0.2e1 * t47 - t66, t1 * t36 + t2 * t39 - t5 * t15 - g(1) * t69 - g(2) * t79 + (qJD(4) - t17) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, 0, -t54 ^ 2, -t54 * t6 - t68 - t95;];
tau_reg  = t3;
