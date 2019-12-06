% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PPRPR3
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
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PPRPR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:22
% EndTime: 2019-12-05 15:05:24
% DurationCPUTime: 0.34s
% Computational Cost: add. (982->92), mult. (1598->148), div. (0->0), fcn. (1219->10), ass. (0->70)
t67 = sin(qJ(5));
t69 = cos(qJ(5));
t72 = qJD(3) ^ 2;
t49 = t67 * t72 * t69;
t45 = qJDD(5) + t49;
t80 = t67 * t45;
t46 = qJDD(5) - t49;
t79 = t69 * t46;
t62 = sin(pkin(7));
t65 = cos(pkin(7));
t44 = -t65 * g(1) - t62 * g(2);
t58 = -g(3) + qJDD(1);
t61 = sin(pkin(8));
t64 = cos(pkin(8));
t33 = t64 * t44 + t61 * t58;
t37 = -t62 * g(1) + t65 * g(2) + qJDD(2);
t68 = sin(qJ(3));
t70 = cos(qJ(3));
t23 = t70 * t33 + t68 * t37;
t21 = -t72 * pkin(3) + t23;
t60 = sin(pkin(9));
t63 = cos(pkin(9));
t22 = -t68 * t33 + t70 * t37;
t73 = qJDD(3) * pkin(3) + t22;
t13 = t63 * t21 + t60 * t73;
t78 = t67 * qJDD(3);
t77 = qJD(3) * qJD(5);
t11 = -t72 * pkin(4) + qJDD(3) * pkin(6) + t13;
t75 = -t61 * t44 + t64 * t58;
t32 = qJDD(4) - t75;
t8 = t67 * t11 - t69 * t32;
t9 = t69 * t11 + t67 * t32;
t4 = t67 * t8 + t69 * t9;
t12 = -t60 * t21 + t63 * t73;
t38 = t63 * qJDD(3) - t60 * t72;
t39 = -t60 * qJDD(3) - t63 * t72;
t76 = -t68 * t38 + t70 * t39;
t53 = t69 * qJDD(3);
t36 = -0.2e1 * t67 * t77 + t53;
t74 = t70 * t38 + t68 * t39;
t35 = 0.2e1 * t69 * t77 + t78;
t71 = qJD(5) ^ 2;
t57 = t69 ^ 2;
t56 = t67 ^ 2;
t55 = t57 * t72;
t54 = t56 * t72;
t48 = -t55 - t71;
t47 = -t54 - t71;
t43 = t54 + t55;
t42 = -t68 * qJDD(3) - t70 * t72;
t41 = t70 * qJDD(3) - t68 * t72;
t40 = (t56 + t57) * qJDD(3);
t31 = -t67 * t47 - t79;
t30 = t69 * t48 - t80;
t29 = -t67 * t46 + t69 * t47;
t28 = t69 * t45 + t67 * t48;
t27 = t64 * t75;
t25 = t63 * t40 - t60 * t43;
t24 = t60 * t40 + t63 * t43;
t20 = t63 * t31 + t60 * t35;
t19 = t63 * t30 - t60 * t36;
t18 = t60 * t31 - t63 * t35;
t17 = t60 * t30 + t63 * t36;
t10 = -qJDD(3) * pkin(4) - t72 * pkin(6) - t12;
t6 = -t60 * t12 + t63 * t13;
t5 = t63 * t12 + t60 * t13;
t3 = t67 * t9 - t69 * t8;
t2 = t60 * t10 + t63 * t4;
t1 = -t63 * t10 + t60 * t4;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61 * t33 + t27, 0, 0, 0, 0, 0, 0, t61 * t42, -t61 * t41, 0, t61 * (-t68 * t22 + t70 * t23) + t27, 0, 0, 0, 0, 0, 0, t61 * t76, -t61 * t74, 0, t61 * (-t68 * t5 + t70 * t6) - t64 * t32, 0, 0, 0, 0, 0, 0, t61 * (-t68 * t17 + t70 * t19) - t64 * t28, t61 * (-t68 * t18 + t70 * t20) - t64 * t29, t61 * (-t68 * t24 + t70 * t25), t61 * (-t68 * t1 + t70 * t2) - t64 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0, 0, 0, 0, 0, t41, t42, 0, t70 * t22 + t68 * t23, 0, 0, 0, 0, 0, 0, t74, t76, 0, t70 * t5 + t68 * t6, 0, 0, 0, 0, 0, 0, t70 * t17 + t68 * t19, t70 * t18 + t68 * t20, t70 * t24 + t68 * t25, t70 * t1 + t68 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t22, -t23, 0, 0, 0, 0, 0, 0, 0, qJDD(3), pkin(3) * t38 + t12, pkin(3) * t39 - t13, 0, pkin(3) * t5, t35 * t67, t69 * t35 + t67 * t36, t80 + t69 * (-t54 + t71), t36 * t69, t67 * (t55 - t71) + t79, 0, pkin(3) * t17 + pkin(4) * t36 + pkin(6) * t30 - t69 * t10, pkin(3) * t18 - pkin(4) * t35 + pkin(6) * t31 + t67 * t10, pkin(3) * t24 + pkin(4) * t43 + pkin(6) * t40 + t4, pkin(3) * t1 - pkin(4) * t10 + pkin(6) * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, t28, t29, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t54 - t55, t78, t49, t53, qJDD(5), -t8, -t9, 0, 0;];
tauJ_reg = t7;
