% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PPPRR1
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
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PPPRR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:58:07
% EndTime: 2019-12-05 14:58:09
% DurationCPUTime: 0.31s
% Computational Cost: add. (786->81), mult. (1280->133), div. (0->0), fcn. (1061->10), ass. (0->66)
t62 = sin(qJ(5));
t64 = cos(qJ(5));
t67 = qJD(4) ^ 2;
t44 = t62 * t67 * t64;
t40 = qJDD(5) + t44;
t74 = t62 * t40;
t41 = qJDD(5) - t44;
t73 = t64 * t41;
t72 = t62 * qJDD(4);
t71 = qJD(4) * qJD(5);
t57 = sin(pkin(7));
t60 = cos(pkin(7));
t39 = -t60 * g(1) - t57 * g(2);
t53 = -g(3) + qJDD(1);
t56 = sin(pkin(8));
t59 = cos(pkin(8));
t30 = t59 * t39 + t56 * t53;
t34 = -t57 * g(1) + t60 * g(2) + qJDD(2);
t55 = sin(pkin(9));
t58 = cos(pkin(9));
t19 = -t55 * t30 + t58 * t34;
t20 = t58 * t30 + t55 * t34;
t63 = sin(qJ(4));
t65 = cos(qJ(4));
t13 = t63 * t19 + t65 * t20;
t11 = -t67 * pkin(4) + qJDD(4) * pkin(6) + t13;
t69 = -t56 * t39 + t59 * t53;
t29 = qJDD(3) - t69;
t8 = t62 * t11 - t64 * t29;
t9 = t64 * t11 + t62 * t29;
t4 = t62 * t8 + t64 * t9;
t12 = t65 * t19 - t63 * t20;
t36 = t65 * qJDD(4) - t63 * t67;
t37 = -t63 * qJDD(4) - t65 * t67;
t70 = -t55 * t36 + t58 * t37;
t47 = t64 * qJDD(4);
t33 = -0.2e1 * t62 * t71 + t47;
t68 = t58 * t36 + t55 * t37;
t32 = 0.2e1 * t64 * t71 + t72;
t66 = qJD(5) ^ 2;
t52 = t64 ^ 2;
t51 = t62 ^ 2;
t50 = t52 * t67;
t48 = t51 * t67;
t43 = -t50 - t66;
t42 = -t48 - t66;
t38 = t48 + t50;
t35 = (t51 + t52) * qJDD(4);
t28 = -t62 * t42 - t73;
t27 = t64 * t43 - t74;
t26 = -t62 * t41 + t64 * t42;
t25 = t64 * t40 + t62 * t43;
t23 = t59 * t29;
t22 = t65 * t35 - t63 * t38;
t21 = t63 * t35 + t65 * t38;
t18 = t65 * t28 + t63 * t32;
t17 = t65 * t27 - t63 * t33;
t16 = t63 * t28 - t65 * t32;
t15 = t63 * t27 + t65 * t33;
t10 = -qJDD(4) * pkin(4) - t67 * pkin(6) - t12;
t6 = -t63 * t12 + t65 * t13;
t5 = t65 * t12 + t63 * t13;
t3 = t62 * t9 - t64 * t8;
t2 = t63 * t10 + t65 * t4;
t1 = -t65 * t10 + t63 * t4;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t30 + t59 * t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * (-t55 * t19 + t58 * t20) - t23, 0, 0, 0, 0, 0, 0, t56 * t70, -t56 * t68, 0, t56 * (-t55 * t5 + t58 * t6) - t23, 0, 0, 0, 0, 0, 0, t56 * (-t55 * t15 + t58 * t17) - t59 * t25, t56 * (-t55 * t16 + t58 * t18) - t59 * t26, t56 * (-t55 * t21 + t58 * t22), t56 * (-t55 * t1 + t58 * t2) - t59 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t19 + t55 * t20, 0, 0, 0, 0, 0, 0, t68, t70, 0, t58 * t5 + t55 * t6, 0, 0, 0, 0, 0, 0, t58 * t15 + t55 * t17, t58 * t16 + t55 * t18, t58 * t21 + t55 * t22, t58 * t1 + t55 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, t25, t26, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4), t12, -t13, 0, 0, t32 * t62, t64 * t32 + t62 * t33, t74 + t64 * (-t48 + t66), t33 * t64, t62 * (t50 - t66) + t73, 0, pkin(4) * t33 + pkin(6) * t27 - t64 * t10, -pkin(4) * t32 + pkin(6) * t28 + t62 * t10, pkin(4) * t38 + pkin(6) * t35 + t4, -pkin(4) * t10 + pkin(6) * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t48 - t50, t72, t44, t47, qJDD(5), -t8, -t9, 0, 0;];
tauJ_reg = t7;
