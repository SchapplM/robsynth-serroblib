% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRRP1
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:55
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:54:33
% EndTime: 2018-11-14 13:54:33
% DurationCPUTime: 0.39s
% Computational Cost: add. (604->107), mult. (1023->129), div. (0->0), fcn. (543->10), ass. (0->78)
t55 = cos(qJ(2));
t88 = t55 * pkin(1);
t41 = qJDD(1) * t88;
t48 = qJDD(1) + qJDD(2);
t52 = sin(qJ(2));
t83 = pkin(1) * qJD(1);
t77 = t52 * t83;
t18 = t48 * pkin(2) - qJD(2) * t77 + t41;
t49 = qJD(1) + qJD(2);
t25 = t49 * pkin(2) + t55 * t83;
t51 = sin(qJ(3));
t72 = qJD(3) * t77;
t26 = t51 * t72;
t54 = cos(qJ(3));
t78 = qJDD(1) * t52;
t82 = qJD(2) * t55;
t64 = (qJD(1) * t82 + t78) * pkin(1);
t7 = t51 * t18 + (qJD(3) * t25 + t64) * t54 - t26;
t50 = qJ(1) + qJ(2);
t46 = qJ(3) + t50;
t37 = sin(t46);
t38 = cos(t46);
t97 = g(1) * t38 + g(2) * t37;
t96 = g(1) * t37 - g(2) * t38;
t44 = sin(t50);
t45 = cos(t50);
t95 = g(1) * t44 - g(2) * t45;
t40 = pkin(2) + t88;
t87 = t51 * t52;
t66 = t54 * t55 - t87;
t79 = qJD(3) * t54;
t80 = qJD(3) * t51;
t11 = t40 * t79 + (t66 * qJD(2) - t52 * t80) * pkin(1);
t16 = t51 * t25 + t54 * t77;
t86 = t52 * t54;
t23 = pkin(1) * t86 + t51 * t40;
t94 = t16 * t11 + t7 * t23;
t42 = qJDD(3) + t48;
t93 = pkin(2) * t42;
t92 = pkin(2) * t44;
t53 = sin(qJ(1));
t89 = t53 * pkin(1);
t85 = g(1) * t45 + g(2) * t44;
t36 = pkin(2) * t45;
t56 = cos(qJ(1));
t84 = t56 * pkin(1) + t36;
t76 = t52 * t79;
t74 = qJD(1) * (-qJD(2) + t49);
t73 = qJD(2) * (-qJD(1) - t49);
t22 = -pkin(1) * t87 + t54 * t40;
t71 = t41 + t95;
t21 = t66 * t83;
t70 = -t16 * t21 + (t16 * t79 + t51 * t7) * pkin(2);
t69 = -pkin(3) * t37 - t92;
t68 = g(1) * t53 - g(2) * t56;
t67 = -t51 * t55 - t86;
t15 = t54 * t25 - t51 * t77;
t43 = qJD(3) + t49;
t63 = (-pkin(2) * t43 - t25) * qJD(3) - t64;
t62 = -t7 + t97;
t17 = t54 * t18;
t8 = -t25 * t80 + t17 + (-t51 * t78 + (-t51 * t82 - t76) * qJD(1)) * pkin(1);
t20 = t67 * t83;
t61 = -t20 * t43 + t63 * t51 + t17 + t96;
t60 = t8 + t96;
t12 = -t40 * t80 + (t67 * qJD(2) - t76) * pkin(1);
t59 = t12 * t43 + t60;
t58 = t16 * t43 + t60;
t39 = t54 * pkin(2) + pkin(3);
t32 = t42 * pkin(3);
t31 = pkin(3) * t38;
t19 = pkin(3) + t22;
t13 = t43 * pkin(3) + t15;
t5 = t32 + t8;
t3 = t15 * t43 + t62;
t2 = t21 * t43 + t26 + (-t18 - t93) * t51 + t63 * t54 + t97;
t1 = -t11 * t43 - t23 * t42 + t62;
t4 = [0, 0, 0, 0, 0, qJDD(1), t68, g(1) * t56 + g(2) * t53, 0, 0, 0, 0, 0, 0, 0, t48 (t48 * t55 + t52 * t73) * pkin(1) + t71 ((-qJDD(1) - t48) * t52 + t55 * t73) * pkin(1) + t85, 0 (t68 + (t52 ^ 2 + t55 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t42, t22 * t42 + t59, t1, 0, t8 * t22 + t15 * t12 - g(1) * (-t89 - t92) - g(2) * t84 + t94, 0, 0, 0, 0, 0, t42, t19 * t42 + t32 + t59, t1, 0, t5 * t19 + t13 * t12 - g(1) * (t69 - t89) - g(2) * (t31 + t84) + t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t52 * pkin(1) * t74 + t71 (t55 * t74 - t78) * pkin(1) + t85, 0, 0, 0, 0, 0, 0, 0, t42 (-t72 + t93) * t54 + t61, t2, 0, -t15 * t20 + (-t15 * t80 + t54 * t8 + t95) * pkin(2) + t70, 0, 0, 0, 0, 0, t42, t39 * t42 - t54 * t72 + t32 + t61, t2, 0, t5 * t39 - g(1) * t69 - g(2) * (t31 + t36) + (-pkin(2) * t80 - t20) * t13 + t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t58, t3, 0, 0, 0, 0, 0, 0, 0, t42, 0.2e1 * t32 + t58, t3, 0 (t13 - t15) * t16 + (t5 + t96) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4) - g(3);];
tau_reg  = t4;
