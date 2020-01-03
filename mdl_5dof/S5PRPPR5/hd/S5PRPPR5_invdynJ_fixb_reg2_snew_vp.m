% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRPPR5
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRPPR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:37
% EndTime: 2019-12-31 17:38:39
% DurationCPUTime: 0.39s
% Computational Cost: add. (1005->101), mult. (1645->133), div. (0->0), fcn. (856->8), ass. (0->68)
t79 = 2 * qJD(2);
t78 = (pkin(2) + pkin(3));
t61 = sin(qJ(5));
t63 = cos(qJ(5));
t66 = qJD(2) ^ 2;
t43 = t61 * t66 * t63;
t39 = qJDD(5) + t43;
t77 = t61 * t39;
t40 = qJDD(5) - t43;
t76 = t63 * t40;
t50 = qJDD(2) * qJ(3);
t57 = sin(pkin(7));
t59 = cos(pkin(7));
t38 = -t59 * g(1) - t57 * g(2);
t54 = -g(3) + qJDD(1);
t62 = sin(qJ(2));
t64 = cos(qJ(2));
t27 = t64 * t38 + t62 * t54;
t73 = (qJD(3) * t79) + t27;
t70 = t50 + t73;
t19 = -(t78 * t66) + t70;
t56 = sin(pkin(8));
t58 = cos(pkin(8));
t55 = qJDD(2) * pkin(2);
t26 = -t62 * t38 + t64 * t54;
t69 = -qJDD(3) + t26;
t23 = -(t66 * qJ(3)) - t55 - t69;
t67 = -qJDD(2) * pkin(3) + t23;
t11 = t58 * t19 + t56 * t67;
t75 = t61 * qJDD(2);
t74 = t63 * qJDD(2);
t72 = qJD(5) * t79;
t71 = t56 * t19 - t58 * t67;
t68 = t57 * g(1) - t59 * g(2) + qJDD(4);
t9 = -(t66 * pkin(4)) - qJDD(2) * pkin(6) + t11;
t6 = t61 * t9 - t63 * t68;
t7 = t61 * t68 + t63 * t9;
t3 = t61 * t6 + t63 * t7;
t30 = t63 * t72 + t75;
t31 = t61 * t72 - t74;
t65 = qJD(5) ^ 2;
t53 = t63 ^ 2;
t52 = t61 ^ 2;
t48 = t53 * t66;
t47 = t52 * t66;
t42 = -t48 - t65;
t41 = -t47 - t65;
t37 = t47 + t48;
t36 = t64 * qJDD(2) - t62 * t66;
t35 = t62 * qJDD(2) + t64 * t66;
t34 = (-t52 - t53) * qJDD(2);
t33 = t58 * qJDD(2) + t56 * t66;
t32 = -t56 * qJDD(2) + t58 * t66;
t25 = -t61 * t41 - t76;
t24 = t63 * t42 - t77;
t22 = t58 * t34 - t56 * t37;
t21 = t56 * t34 + t58 * t37;
t20 = -(t66 * pkin(2)) + t70;
t15 = t58 * t25 - t56 * t30;
t14 = t58 * t24 - t56 * t31;
t13 = t56 * t25 + t58 * t30;
t12 = t56 * t24 + t58 * t31;
t8 = qJDD(2) * pkin(4) - t66 * pkin(6) + t71;
t5 = t58 * t11 + t56 * t71;
t4 = t56 * t11 - t58 * t71;
t2 = t58 * t3 + t56 * t8;
t1 = t56 * t3 - t58 * t8;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, 0, 0, 0, 0, 0, t36, -t35, 0, t64 * t26 + t62 * t27, 0, 0, 0, 0, 0, 0, t36, 0, t35, t62 * t20 - t64 * t23, 0, 0, 0, 0, 0, 0, -t62 * t32 + t64 * t33, t64 * t32 + t62 * t33, 0, -t64 * t4 + t62 * t5, 0, 0, 0, 0, 0, 0, -t64 * t12 + t62 * t14, -t64 * t13 + t62 * t15, -t64 * t21 + t62 * t22, -t64 * t1 + t62 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t26, -t27, 0, 0, 0, 0, 0, qJDD(2), 0, 0, 0.2e1 * t55 + t69, 0, 0.2e1 * t50 + t73, -pkin(2) * t23 + qJ(3) * t20, 0, 0, 0, 0, 0, qJDD(2), -qJ(3) * t32 + t78 * t33 + t71, qJ(3) * t33 + t78 * t32 + t11, 0, qJ(3) * t5 - t78 * t4, t30 * t61, t63 * t30 - t61 * t31, -t77 - t63 * (-t47 + t65), -t31 * t63, -t61 * (t48 - t65) - t76, 0, -pkin(4) * t31 - pkin(6) * t24 + qJ(3) * t14 - t78 * t12 + t63 * t8, -pkin(4) * t30 - pkin(6) * t25 + qJ(3) * t15 - t78 * t13 - t61 * t8, -pkin(4) * t37 - pkin(6) * t34 + qJ(3) * t22 - t78 * t21 - t3, pkin(4) * t8 - pkin(6) * t3 + qJ(3) * t2 - t78 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), 0, -t66, t23, 0, 0, 0, 0, 0, 0, -t33, -t32, 0, t4, 0, 0, 0, 0, 0, 0, t12, t13, t21, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, 0, 0, 0, 0, 0, t63 * t39 + t61 * t42, -t61 * t40 + t63 * t41, 0, -t63 * t6 + t61 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, t47 - t48, -t75, t43, -t74, qJDD(5), -t6, -t7, 0, 0;];
tauJ_reg = t10;
