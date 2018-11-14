% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% 
% Output:
% tau_reg [3x(3*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:16
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S3RRR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:16:00
% EndTime: 2018-11-14 10:16:00
% DurationCPUTime: 0.28s
% Computational Cost: add. (290->80), mult. (520->109), div. (0->0), fcn. (278->10), ass. (0->59)
t31 = qJD(1) + qJD(2);
t37 = cos(qJ(2));
t60 = pkin(1) * qJD(1);
t13 = t31 * pkin(2) + t37 * t60;
t33 = sin(qJ(3));
t34 = sin(qJ(2));
t54 = t34 * t60;
t49 = qJD(3) * t54;
t14 = t33 * t49;
t36 = cos(qJ(3));
t55 = qJDD(1) * t34;
t59 = qJD(2) * t37;
t43 = (qJD(1) * t59 + t55) * pkin(1);
t64 = t37 * pkin(1);
t24 = qJDD(1) * t64;
t30 = qJDD(1) + qJDD(2);
t8 = t30 * pkin(2) - qJD(2) * t54 + t24;
t1 = t33 * t8 + (qJD(3) * t13 + t43) * t36 - t14;
t32 = qJ(1) + qJ(2);
t29 = qJ(3) + t32;
t21 = sin(t29);
t22 = cos(t29);
t70 = g(1) * t22 + g(2) * t21;
t69 = g(1) * t21 - g(2) * t22;
t27 = sin(t32);
t28 = cos(t32);
t68 = g(1) * t27 - g(2) * t28;
t25 = qJDD(3) + t30;
t67 = pkin(2) * t25;
t63 = t33 * t34;
t62 = t34 * t36;
t61 = g(1) * t28 + g(2) * t27;
t57 = qJD(3) * t33;
t56 = qJD(3) * t36;
t53 = t34 * t56;
t51 = qJD(1) * (-qJD(2) + t31);
t50 = qJD(2) * (-qJD(1) - t31);
t48 = t24 + t68;
t35 = sin(qJ(1));
t38 = cos(qJ(1));
t47 = g(1) * t35 - g(2) * t38;
t46 = -t33 * t37 - t62;
t45 = t36 * t37 - t63;
t26 = qJD(3) + t31;
t42 = (-pkin(2) * t26 - t13) * qJD(3) - t43;
t41 = -t1 + t70;
t7 = t36 * t8;
t2 = -t13 * t57 + t7 + (-t33 * t55 + (-t33 * t59 - t53) * qJD(1)) * pkin(1);
t40 = t2 + t69;
t23 = pkin(2) + t64;
t12 = pkin(1) * t62 + t33 * t23;
t11 = -pkin(1) * t63 + t36 * t23;
t10 = t45 * t60;
t9 = t46 * t60;
t6 = t33 * t13 + t36 * t54;
t5 = t36 * t13 - t33 * t54;
t4 = -t23 * t57 + (qJD(2) * t46 - t53) * pkin(1);
t3 = t23 * t56 + (qJD(2) * t45 - t34 * t57) * pkin(1);
t15 = [0, 0, 0, 0, 0, qJDD(1), t47, g(1) * t38 + g(2) * t35, 0, 0, 0, 0, 0, 0, 0, t30 (t30 * t37 + t34 * t50) * pkin(1) + t48 ((-qJDD(1) - t30) * t34 + t37 * t50) * pkin(1) + t61, 0 (t47 + (t34 ^ 2 + t37 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t25, t11 * t25 + t4 * t26 + t40, -t12 * t25 - t3 * t26 + t41, 0, t1 * t12 + t6 * t3 + t2 * t11 + t5 * t4 - g(1) * (-t35 * pkin(1) - pkin(2) * t27) - g(2) * (t38 * pkin(1) + pkin(2) * t28); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, pkin(1) * t34 * t51 + t48 (t37 * t51 - t55) * pkin(1) + t61, 0, 0, 0, 0, 0, 0, 0, t25, -t9 * t26 + t7 + (-t49 + t67) * t36 + t42 * t33 + t69, t10 * t26 + t14 + (-t8 - t67) * t33 + t42 * t36 + t70, 0, -t6 * t10 - t5 * t9 + (t1 * t33 + t2 * t36 + (-t33 * t5 + t36 * t6) * qJD(3) + t68) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t6 * t26 + t40, t5 * t26 + t41, 0, 0;];
tau_reg  = t15;
