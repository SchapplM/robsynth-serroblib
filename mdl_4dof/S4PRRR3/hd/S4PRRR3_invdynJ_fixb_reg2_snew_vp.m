% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4PRRR3
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4PRRR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:39
% EndTime: 2019-12-31 16:31:40
% DurationCPUTime: 0.23s
% Computational Cost: add. (541->61), mult. (841->87), div. (0->0), fcn. (572->8), ass. (0->56)
t55 = sin(pkin(7));
t56 = cos(pkin(7));
t37 = t55 * g(1) - t56 * g(2);
t38 = -t56 * g(1) - t55 * g(2);
t59 = sin(qJ(2));
t62 = cos(qJ(2));
t65 = t62 * t37 - t59 * t38;
t20 = qJDD(2) * pkin(2) + t65;
t64 = -t59 * t37 - t62 * t38;
t21 = -qJD(2) ^ 2 * pkin(2) - t64;
t58 = sin(qJ(3));
t61 = cos(qJ(3));
t13 = t58 * t20 + t61 * t21;
t51 = qJD(2) + qJD(3);
t49 = t51 ^ 2;
t50 = qJDD(2) + qJDD(3);
t11 = -t49 * pkin(3) + t50 * pkin(6) + t13;
t54 = -g(3) + qJDD(1);
t57 = sin(qJ(4));
t60 = cos(qJ(4));
t5 = t57 * t11 - t60 * t54;
t6 = t60 * t11 + t57 * t54;
t2 = t57 * t5 + t60 * t6;
t12 = t61 * t20 - t58 * t21;
t10 = -t50 * pkin(3) - t49 * pkin(6) - t12;
t73 = -pkin(3) * t10 + pkin(6) * t2;
t41 = t57 * t49 * t60;
t35 = qJDD(4) + t41;
t72 = t57 * t35;
t71 = t57 * t50;
t36 = qJDD(4) - t41;
t70 = t60 * t36;
t69 = qJD(4) * t51;
t52 = t57 ^ 2;
t44 = t52 * t49;
t63 = qJD(4) ^ 2;
t39 = -t44 - t63;
t25 = -t57 * t39 - t70;
t29 = 0.2e1 * t60 * t69 + t71;
t68 = -pkin(3) * t29 + pkin(6) * t25 + t57 * t10;
t53 = t60 ^ 2;
t45 = t53 * t49;
t40 = -t45 - t63;
t24 = t60 * t40 - t72;
t43 = t60 * t50;
t30 = -0.2e1 * t57 * t69 + t43;
t67 = pkin(3) * t30 + pkin(6) * t24 - t60 * t10;
t32 = (t52 + t53) * t50;
t33 = t44 + t45;
t66 = pkin(3) * t33 + pkin(6) * t32 + t2;
t23 = t72 + t60 * (-t44 + t63);
t22 = t57 * (t45 - t63) + t70;
t17 = t29 * t57;
t16 = t30 * t60;
t14 = t60 * t29 + t57 * t30;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, 0, 0, 0, 0, 0, t60 * t35 + t57 * t40, -t57 * t36 + t60 * t39, 0, -t60 * t5 + t57 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t65, t64, 0, 0, 0, 0, 0, 0, 0, t50, pkin(2) * (-t58 * t49 + t61 * t50) + t12, pkin(2) * (-t61 * t49 - t58 * t50) - t13, 0, pkin(2) * (t61 * t12 + t58 * t13), t17, t14, t23, t16, t22, 0, pkin(2) * (t58 * t24 + t61 * t30) + t67, pkin(2) * (t58 * t25 - t61 * t29) + t68, pkin(2) * (t58 * t32 + t61 * t33) + t66, pkin(2) * (-t61 * t10 + t58 * t2) + t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t12, -t13, 0, 0, t17, t14, t23, t16, t22, 0, t67, t68, t66, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t44 - t45, t71, t41, t43, qJDD(4), -t5, -t6, 0, 0;];
tauJ_reg = t1;
