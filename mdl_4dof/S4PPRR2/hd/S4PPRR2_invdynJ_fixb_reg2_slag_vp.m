% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PPRR2
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
%   pkin=[a2,a3,a4,d3,d4,theta2]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PPRR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:17:02
% EndTime: 2019-03-08 18:17:02
% DurationCPUTime: 0.26s
% Computational Cost: add. (317->69), mult. (720->84), div. (0->0), fcn. (624->10), ass. (0->47)
t36 = sin(pkin(6));
t37 = cos(pkin(6));
t39 = sin(qJ(3));
t41 = cos(qJ(3));
t59 = -t39 * t36 + t41 * t37;
t16 = t59 * qJD(1);
t33 = pkin(6) + qJ(3);
t31 = qJ(4) + t33;
t27 = sin(t31);
t28 = cos(t31);
t58 = g(1) * t28 + g(2) * t27;
t57 = g(1) * t27 - g(2) * t28;
t21 = t41 * t36 + t39 * t37;
t17 = t21 * qJD(1);
t56 = t21 * qJDD(1);
t19 = t21 * qJD(3);
t44 = t59 * qJDD(1);
t42 = -qJD(1) * t19 + t44;
t10 = qJDD(3) * pkin(3) + t42;
t18 = t59 * qJD(3);
t11 = qJD(1) * t18 + t56;
t38 = sin(qJ(4));
t53 = t38 * t17;
t14 = qJD(4) * t53;
t15 = qJD(3) * pkin(3) + t16;
t40 = cos(qJ(4));
t1 = t38 * t10 - t14 + (qJD(4) * t15 + t11) * t40;
t32 = qJDD(3) + qJDD(4);
t55 = pkin(3) * t32;
t51 = t40 * t17;
t34 = qJD(3) + qJD(4);
t49 = -pkin(3) * t34 - t15;
t48 = t40 * t10 - t38 * t11;
t29 = sin(t33);
t30 = cos(t33);
t45 = g(1) * t29 - g(2) * t30;
t6 = t38 * t15 + t51;
t12 = -t38 * t21 + t40 * t59;
t13 = t40 * t21 + t38 * t59;
t2 = -qJD(4) * t6 + t48;
t35 = qJDD(2) - g(3);
t9 = t40 * t16 - t53;
t8 = -t38 * t16 - t51;
t5 = t40 * t15 - t53;
t4 = -t13 * qJD(4) - t38 * t18 - t40 * t19;
t3 = t12 * qJD(4) + t40 * t18 - t38 * t19;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) + (t36 ^ 2 + t37 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -t19 * qJD(3) + qJDD(3) * t59, -t18 * qJD(3) - t21 * qJDD(3), 0, t11 * t21 - t16 * t19 + t17 * t18 + t42 * t59 - g(2), 0, 0, 0, 0, 0, 0, t12 * t32 + t4 * t34, -t13 * t32 - t3 * t34, 0, t1 * t13 + t2 * t12 + t6 * t3 + t5 * t4 - g(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t44 + t45, g(1) * t30 + g(2) * t29 - t56, 0, 0, 0, 0, 0, 0, 0, t32, t40 * t55 - t8 * t34 + (t49 * t38 - t51) * qJD(4) + t48 + t57, t9 * t34 + t14 + (-t10 - t55) * t38 + (t49 * qJD(4) - t11) * t40 + t58, 0, -t5 * t8 - t6 * t9 + (t1 * t38 + t2 * t40 + (-t38 * t5 + t40 * t6) * qJD(4) + t45) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t6 * t34 + t2 + t57, t5 * t34 - t1 + t58, 0, 0;];
tau_reg  = t7;
