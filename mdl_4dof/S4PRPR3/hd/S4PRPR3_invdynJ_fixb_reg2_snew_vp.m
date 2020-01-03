% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4PRPR3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4PRPR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:55
% EndTime: 2019-12-31 16:20:57
% DurationCPUTime: 0.51s
% Computational Cost: add. (805->107), mult. (1842->164), div. (0->0), fcn. (1271->8), ass. (0->75)
t67 = qJD(2) ^ 2;
t85 = sin(pkin(6));
t86 = cos(pkin(6));
t44 = -t86 * g(1) - t85 * g(2);
t64 = sin(qJ(2));
t72 = t85 * g(1) - t86 * g(2);
t93 = cos(qJ(2));
t70 = -t93 * t44 - t64 * t72;
t98 = -t67 * pkin(2) + qJDD(2) * qJ(3) + (2 * qJD(2) * qJD(3)) - t70;
t61 = sin(pkin(7));
t56 = t61 ^ 2;
t62 = cos(pkin(7));
t57 = t62 ^ 2;
t101 = t56 + t57;
t63 = sin(qJ(4));
t65 = cos(qJ(4));
t39 = (t61 * t63 - t62 * t65) * qJD(2);
t73 = t61 * t65 + t62 * t63;
t41 = t73 * qJD(2);
t31 = t41 * t39;
t96 = qJDD(4) - t31;
t100 = t63 * t96;
t99 = t65 * t96;
t59 = -g(3) + qJDD(1);
t50 = t62 * t59;
t97 = t50 + (pkin(3) * t62 * t67 - pkin(5) * qJDD(2) - t98) * t61;
t58 = t67 * qJ(3);
t60 = qJDD(2) * pkin(2);
t75 = -t64 * t44 + t93 * t72;
t25 = qJDD(3) - t58 - t60 - t75;
t95 = t101 * t58 + t25 - t60;
t36 = t39 ^ 2;
t37 = t41 ^ 2;
t18 = t61 * t59 + t98 * t62;
t81 = t62 * qJDD(2);
t92 = t57 * t67;
t11 = -pkin(3) * t92 + pkin(5) * t81 + t18;
t4 = t63 * t11 - t65 * t97;
t5 = t65 * t11 + t97 * t63;
t1 = -t65 * t4 + t63 * t5;
t94 = t61 * t1;
t15 = -pkin(3) * t81 + t25 + (-t56 * t67 - t92) * pkin(5);
t91 = t63 * t15;
t23 = qJDD(4) + t31;
t90 = t63 * t23;
t89 = t65 * t15;
t88 = t65 * t23;
t84 = t39 * qJD(4);
t83 = t41 * qJD(4);
t82 = t61 * qJDD(2);
t2 = t63 * t4 + t65 * t5;
t17 = t98 * t61 - t50;
t76 = t61 * t17 + t62 * t18;
t16 = -t63 * t82 + t65 * t81;
t38 = t73 * qJDD(2);
t66 = qJD(4) ^ 2;
t55 = t57 * qJDD(2);
t54 = t56 * qJDD(2);
t43 = t101 * t67;
t34 = -t37 - t66;
t33 = -t37 + t66;
t32 = t36 - t66;
t29 = t38 - t84;
t28 = t38 - 0.2e1 * t84;
t27 = t16 - t83;
t26 = -t16 + 0.2e1 * t83;
t21 = -t66 - t36;
t19 = -t36 - t37;
t13 = -t63 * t34 - t88;
t12 = t65 * t34 - t90;
t9 = t65 * t16 + t63 * t38;
t8 = t63 * t16 - t65 * t38;
t7 = t65 * t21 - t100;
t6 = t63 * t21 + t99;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62 * t17 + t61 * t18, 0, 0, 0, 0, 0, 0, t62 * t6 + t61 * t7, t62 * t12 + t61 * t13, t61 * t9 + t62 * t8, t62 * t1 + t61 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t75, t70, 0, 0, t54, 0.2e1 * t61 * t81, 0, t55, 0, 0, -t95 * t62, t95 * t61, pkin(2) * t43 + qJ(3) * (t55 + t54) + t76, -pkin(2) * t25 + qJ(3) * t76, t61 * (t65 * t29 - t63 * t83) + t62 * (t63 * t29 + t65 * t83), t61 * (-t65 * t26 - t63 * t28) + t62 * (-t63 * t26 + t65 * t28), t61 * (-t63 * t33 + t99) + t62 * (t65 * t33 + t100), t61 * (-t63 * t27 + t65 * t84) + t62 * (t65 * t27 + t63 * t84), t61 * (t65 * t32 - t90) + t62 * (t63 * t32 + t88), (t61 * (-t39 * t65 + t41 * t63) + t62 * (-t39 * t63 - t41 * t65)) * qJD(4), t61 * (-pkin(5) * t6 + t91) + t62 * (-pkin(3) * t26 + pkin(5) * t7 - t89) - pkin(2) * t26 + qJ(3) * (-t61 * t6 + t62 * t7), t61 * (-pkin(5) * t12 + t89) + t62 * (-pkin(3) * t28 + pkin(5) * t13 + t91) - pkin(2) * t28 + qJ(3) * (-t61 * t12 + t62 * t13), t61 * (-pkin(5) * t8 - t1) + t62 * (-pkin(3) * t19 + pkin(5) * t9 + t2) - pkin(2) * t19 + qJ(3) * (-t61 * t8 + t62 * t9), -pkin(5) * t94 + t62 * (-pkin(3) * t15 + pkin(5) * t2) - pkin(2) * t15 + qJ(3) * (t62 * t2 - t94); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, t82, -t43, t25, 0, 0, 0, 0, 0, 0, t26, t28, t19, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t37 - t36, t38, -t31, t16, qJDD(4), -t4, -t5, 0, 0;];
tauJ_reg = t3;
