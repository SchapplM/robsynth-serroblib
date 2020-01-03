% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PPRR3
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PPRR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:29
% EndTime: 2019-12-31 16:17:30
% DurationCPUTime: 0.24s
% Computational Cost: add. (230->76), mult. (517->109), div. (0->0), fcn. (355->6), ass. (0->60)
t25 = cos(qJ(3));
t27 = qJD(3) ^ 2;
t64 = t27 * t25;
t22 = sin(qJ(4));
t19 = t22 ^ 2;
t24 = cos(qJ(4));
t20 = t24 ^ 2;
t63 = t19 - t20;
t62 = t19 + t20;
t26 = qJD(4) ^ 2;
t61 = t26 + t27;
t60 = qJD(3) * pkin(3);
t23 = sin(qJ(3));
t11 = qJD(3) * pkin(5) + t23 * qJD(2);
t3 = -t24 * qJD(1) - t22 * t11;
t59 = t3 * qJD(4);
t58 = cos(pkin(6));
t57 = sin(pkin(6));
t54 = t25 * qJD(2);
t12 = -t54 - t60;
t56 = qJD(3) * t12;
t55 = qJDD(3) * pkin(3);
t53 = qJDD(4) * t22;
t52 = t23 * qJDD(2);
t51 = t24 * qJDD(3);
t50 = t25 * qJDD(2);
t49 = t25 * qJDD(3);
t48 = qJD(1) * qJD(4);
t47 = qJD(2) * qJD(3);
t46 = qJD(3) * qJD(4);
t45 = t22 * t27 * t24;
t44 = t22 * t46;
t43 = t23 * t47;
t42 = t25 * t47;
t41 = t62 * qJDD(3);
t40 = -qJD(4) * t11 - qJDD(1);
t39 = t24 * t44;
t7 = -t57 * t23 - t58 * t25;
t8 = t58 * t23 - t57 * t25;
t38 = g(1) * t8 - g(2) * t7;
t37 = g(1) * t7 + g(2) * t8;
t4 = -t22 * qJD(1) + t24 * t11;
t36 = t22 * t3 - t24 * t4;
t35 = -g(1) * t57 + g(2) * t58;
t34 = -g(3) - t40;
t5 = t43 - t50 - t55;
t33 = t38 - t5;
t6 = qJDD(3) * pkin(5) + t42 + t52;
t32 = -t37 - t6 - t56;
t31 = -pkin(5) * qJDD(4) + (t12 + t54 - t60) * qJD(4);
t1 = -t22 * qJDD(1) + t24 * t6 + t59;
t16 = t22 * t48;
t2 = -t22 * t6 + t40 * t24 + t16;
t30 = t1 * t24 + (-t22 * t4 - t24 * t3) * qJD(4) - t2 * t22;
t29 = -pkin(5) * t26 + t33 + t43 + t55;
t28 = t30 + t37;
t21 = qJDD(1) - g(3);
t10 = qJDD(4) * t24 - t26 * t22;
t9 = t26 * t24 + t53;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, -t10, t9, 0, t36 * qJD(4) - t1 * t22 - t2 * t24 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) + t35, 0, 0, 0, 0, 0, 0, -t27 * t23 + t49, -qJDD(3) * t23 - t64, 0, (t23 ^ 2 + t25 ^ 2) * qJDD(2) + t35, 0, 0, 0, 0, 0, 0, (-0.2e1 * t44 + t51) * t25 + (-t61 * t24 - t53) * t23, (-qJDD(4) * t23 - 0.2e1 * t25 * t46) * t24 + (t61 * t23 - t49) * t22, t23 * t41 + t62 * t64, (-t36 * qJD(3) - t5) * t25 + (t30 + t56) * t23 + t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t38 + t50, -t37 - t52, 0, 0, t19 * qJDD(3) + 0.2e1 * t39, 0.2e1 * t22 * t51 - 0.2e1 * t63 * t46, t9, t20 * qJDD(3) - 0.2e1 * t39, t10, 0, t31 * t22 + t29 * t24, -t29 * t22 + t31 * t24, pkin(5) * t41 - t62 * t42 + t28, (-t12 * t23 + t36 * t25) * qJD(2) + t33 * pkin(3) + t28 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t63 * t27, t22 * qJDD(3), t45, t51, qJDD(4), t4 * qJD(4) + t32 * t22 - t34 * t24 + t16, t59 + t34 * t22 + (t32 + t48) * t24, 0, 0;];
tau_reg = t13;
