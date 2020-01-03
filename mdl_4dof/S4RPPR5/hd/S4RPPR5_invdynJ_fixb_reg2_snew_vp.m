% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RPPR5
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RPPR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:49
% EndTime: 2019-12-31 16:39:50
% DurationCPUTime: 0.28s
% Computational Cost: add. (600->86), mult. (1059->107), div. (0->0), fcn. (450->6), ass. (0->57)
t65 = 2 * qJD(1);
t64 = (pkin(1) + pkin(2));
t43 = sin(qJ(4));
t45 = cos(qJ(4));
t48 = qJD(1) ^ 2;
t30 = t45 * t48 * t43;
t26 = qJDD(4) + t30;
t63 = t43 * t26;
t27 = qJDD(4) - t30;
t62 = t45 * t27;
t37 = qJDD(1) * qJ(2);
t44 = sin(qJ(1));
t46 = cos(qJ(1));
t54 = t46 * g(1) + t44 * g(2);
t52 = (qJD(2) * t65) - t54;
t51 = t37 + t52;
t18 = -(t64 * t48) + t51;
t40 = sin(pkin(6));
t41 = cos(pkin(6));
t56 = t44 * g(1) - t46 * g(2);
t53 = qJDD(2) - t56;
t50 = -(t48 * qJ(2)) + t53;
t49 = -(t64 * qJDD(1)) + t50;
t9 = t41 * t18 + t40 * t49;
t61 = qJDD(1) * pkin(1);
t60 = g(3) + qJDD(3);
t59 = t43 * qJDD(1);
t58 = t45 * qJDD(1);
t57 = qJD(4) * t65;
t55 = t40 * t18 - t41 * t49;
t7 = -(t48 * pkin(3)) - (qJDD(1) * pkin(5)) + t9;
t4 = t43 * t7 - t45 * t60;
t5 = t43 * t60 + t45 * t7;
t2 = t43 * t4 + t45 * t5;
t20 = t45 * t57 + t59;
t21 = t43 * t57 - t58;
t47 = qJD(4) ^ 2;
t39 = t45 ^ 2;
t38 = t43 ^ 2;
t34 = t39 * t48;
t33 = t38 * t48;
t29 = -t34 - t47;
t28 = -t33 - t47;
t25 = t33 + t34;
t24 = (-t38 - t39) * qJDD(1);
t23 = t41 * qJDD(1) + t40 * t48;
t22 = -t40 * qJDD(1) + t41 * t48;
t19 = -t50 + t61;
t15 = -t43 * t28 - t62;
t14 = t45 * t29 - t63;
t12 = t40 * t24 + t41 * t25;
t11 = t40 * t15 + t41 * t20;
t10 = t40 * t14 + t41 * t21;
t6 = (qJDD(1) * pkin(3)) - (t48 * pkin(5)) + t55;
t3 = t40 * t9 - t41 * t55;
t1 = t40 * t2 - t41 * t6;
t8 = [0, 0, 0, 0, 0, qJDD(1), t56, t54, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -t53 + (2 * t61), 0, 0.2e1 * t37 + t52, pkin(1) * t19 + qJ(2) * (-(t48 * pkin(1)) + t51), 0, 0, 0, 0, 0, qJDD(1), -qJ(2) * t22 + t64 * t23 + t55, qJ(2) * t23 + t64 * t22 + t9, 0, qJ(2) * (t40 * t55 + t41 * t9) - t64 * t3, t20 * t43, t45 * t20 - t43 * t21, -t63 - t45 * (-t33 + t47), -t21 * t45, -t43 * (t34 - t47) - t62, 0, qJ(2) * (t41 * t14 - t40 * t21) + t45 * t6 - pkin(3) * t21 - pkin(5) * t14 - t64 * t10, qJ(2) * (t41 * t15 - t40 * t20) - t43 * t6 - pkin(3) * t20 - pkin(5) * t15 - t64 * t11, qJ(2) * (t41 * t24 - t40 * t25) - pkin(3) * t25 - pkin(5) * t24 - t64 * t12 - t2, qJ(2) * (t41 * t2 + t40 * t6) + pkin(3) * t6 - pkin(5) * t2 - t64 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t48, -t19, 0, 0, 0, 0, 0, 0, -t23, -t22, 0, t3, 0, 0, 0, 0, 0, 0, t10, t11, t12, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, 0, 0, 0, 0, t45 * t26 + t43 * t29, -t43 * t27 + t45 * t28, 0, -t45 * t4 + t43 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t33 - t34, -t59, t30, -t58, qJDD(4), -t4, -t5, 0, 0;];
tauJ_reg = t8;
