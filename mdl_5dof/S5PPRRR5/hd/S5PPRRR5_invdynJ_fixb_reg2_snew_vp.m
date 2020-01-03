% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PPRRR5
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PPRRR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:45
% EndTime: 2019-12-31 17:35:47
% DurationCPUTime: 0.35s
% Computational Cost: add. (828->82), mult. (1181->117), div. (0->0), fcn. (822->8), ass. (0->64)
t65 = sin(qJ(5));
t68 = cos(qJ(5));
t63 = sin(pkin(8));
t64 = cos(pkin(8));
t45 = -t63 * g(1) + t64 * g(2) + qJDD(2);
t46 = -t64 * g(1) - t63 * g(2);
t67 = sin(qJ(3));
t70 = cos(qJ(3));
t31 = t67 * t45 + t70 * t46;
t72 = qJD(3) ^ 2;
t21 = -t72 * pkin(3) + t31;
t66 = sin(qJ(4));
t69 = cos(qJ(4));
t30 = t70 * t45 - t67 * t46;
t73 = qJDD(3) * pkin(3) + t30;
t15 = t69 * t21 + t66 * t73;
t59 = qJD(3) + qJD(4);
t57 = t59 ^ 2;
t58 = qJDD(3) + qJDD(4);
t13 = -t57 * pkin(4) + t58 * pkin(7) + t15;
t62 = g(3) - qJDD(1);
t7 = t65 * t13 - t68 * t62;
t8 = t68 * t13 + t65 * t62;
t3 = t65 * t7 + t68 * t8;
t14 = -t66 * t21 + t69 * t73;
t12 = -t58 * pkin(4) - t57 * pkin(7) - t14;
t82 = -pkin(4) * t12 + pkin(7) * t3;
t49 = t65 * t57 * t68;
t43 = qJDD(5) + t49;
t81 = t65 * t43;
t80 = t65 * t58;
t44 = qJDD(5) - t49;
t79 = t68 * t44;
t78 = qJD(5) * t59;
t60 = t65 ^ 2;
t52 = t60 * t57;
t71 = qJD(5) ^ 2;
t47 = -t52 - t71;
t29 = -t65 * t47 - t79;
t35 = 0.2e1 * t68 * t78 + t80;
t77 = -pkin(4) * t35 + pkin(7) * t29 + t65 * t12;
t61 = t68 ^ 2;
t53 = t61 * t57;
t48 = -t53 - t71;
t28 = t68 * t48 - t81;
t51 = t68 * t58;
t36 = -0.2e1 * t65 * t78 + t51;
t76 = pkin(4) * t36 + pkin(7) * t28 - t68 * t12;
t38 = (t60 + t61) * t58;
t42 = t52 + t53;
t75 = pkin(4) * t42 + pkin(7) * t38 + t3;
t39 = -t69 * t57 - t66 * t58;
t74 = t66 * t57 - t69 * t58;
t27 = t81 + t68 * (-t52 + t71);
t26 = t65 * (t53 - t71) + t79;
t23 = t35 * t65;
t22 = t36 * t68;
t19 = t66 * t38 + t69 * t42;
t18 = t68 * t35 + t65 * t36;
t17 = t66 * t29 - t69 * t35;
t16 = t66 * t28 + t69 * t36;
t4 = t69 * t14 + t66 * t15;
t1 = -t69 * t12 + t66 * t3;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, 0, 0, 0, 0, 0, 0, -t68 * t43 - t65 * t48, t65 * t44 - t68 * t47, 0, -t65 * t8 + t68 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, t70 * qJDD(3) - t67 * t72, -t67 * qJDD(3) - t70 * t72, 0, t70 * t30 + t67 * t31, 0, 0, 0, 0, 0, 0, t67 * t39 - t70 * t74, t70 * t39 + t67 * t74, 0, t67 * (-t66 * t14 + t69 * t15) + t70 * t4, 0, 0, 0, 0, 0, 0, t67 * (t69 * t28 - t66 * t36) + t70 * t16, t67 * (t69 * t29 + t66 * t35) + t70 * t17, t67 * (t69 * t38 - t66 * t42) + t70 * t19, t67 * (t66 * t12 + t69 * t3) + t70 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t30, -t31, 0, 0, 0, 0, 0, 0, 0, t58, -pkin(3) * t74 + t14, pkin(3) * t39 - t15, 0, pkin(3) * t4, t23, t18, t27, t22, t26, 0, pkin(3) * t16 + t76, pkin(3) * t17 + t77, pkin(3) * t19 + t75, pkin(3) * t1 + t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t14, -t15, 0, 0, t23, t18, t27, t22, t26, 0, t76, t77, t75, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t52 - t53, t80, t49, t51, qJDD(5), -t7, -t8, 0, 0;];
tauJ_reg = t2;
