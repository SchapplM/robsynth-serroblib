% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPRR1
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [4x10]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:51
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:50:35
% EndTime: 2018-11-14 13:50:36
% DurationCPUTime: 0.21s
% Computational Cost: add. (352->78), mult. (620->93), div. (0->0), fcn. (356->12), ass. (0->51)
t32 = qJ(1) + pkin(7) + qJ(3);
t26 = sin(t32);
t27 = cos(t32);
t66 = g(1) * t26 - g(2) * t27;
t36 = cos(pkin(7));
t28 = t36 * pkin(1) + pkin(2);
t41 = cos(qJ(3));
t17 = t41 * t28;
t38 = sin(qJ(3));
t35 = sin(pkin(7));
t63 = pkin(1) * t35;
t65 = -t38 * t63 + t17;
t34 = qJD(1) + qJD(3);
t31 = qJD(4) + t34;
t64 = qJD(4) - t31;
t14 = t28 * qJDD(1);
t56 = qJD(1) * t63;
t51 = t38 * t56;
t16 = t28 * qJD(1);
t59 = qJD(3) * t16;
t3 = t41 * (qJDD(1) * t63 + t59) - qJD(3) * t51 + t38 * t14;
t33 = qJDD(1) + qJDD(3);
t30 = qJDD(4) + t33;
t62 = pkin(3) * t30;
t40 = cos(qJ(4));
t7 = t38 * t16 + t41 * t56;
t60 = t40 * t7;
t29 = qJ(4) + t32;
t24 = sin(t29);
t25 = cos(t29);
t37 = sin(qJ(4));
t57 = qJD(4) * t37 * t7 + g(1) * t25 + g(2) * t24;
t6 = t41 * t16 - t51;
t4 = t34 * pkin(3) + t6;
t55 = -pkin(3) * t31 - t4;
t49 = t41 * t14 - t38 * t59;
t53 = qJD(1) * qJD(3) * t41;
t44 = (-qJDD(1) * t38 - t53) * t63 + t49;
t2 = t33 * pkin(3) + t44;
t54 = -t7 * t31 - t2;
t39 = sin(qJ(1));
t42 = cos(qJ(1));
t50 = g(1) * t39 - g(2) * t42;
t48 = -t37 * t4 - t60;
t11 = t38 * t28 + t41 * t63;
t47 = g(1) * t24 - g(2) * t25 + t40 * t2 - t37 * t3;
t45 = g(1) * t27 + g(2) * t26 - t3;
t10 = pkin(3) + t65;
t9 = t11 * qJD(3);
t8 = t65 * qJD(3);
t1 = [qJDD(1), t50, g(1) * t42 + g(2) * t39 (t50 + (t35 ^ 2 + t36 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t33, t17 * t33 - t9 * t34 + (-t53 + (-qJDD(1) - t33) * t38) * t63 + t49 + t66, -t11 * t33 - t8 * t34 + t45, t30 (-t37 * t8 - t40 * t9) * t31 + (t40 * t10 - t37 * t11) * t30 + ((-t37 * t10 - t40 * t11) * t31 + t48) * qJD(4) + t47 (-(-qJD(4) * t11 - t9) * t31 - t10 * t30 - t2) * t37 + (-(qJD(4) * t10 + t8) * t31 - t11 * t30 - t3 - qJD(4) * t4) * t40 + t57; 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t33, t7 * t34 + t44 + t66, t6 * t34 + t45, t30, t40 * t62 - (-t37 * t6 - t60) * t31 + (t55 * t37 - t60) * qJD(4) + t47 (t54 - t62) * t37 + (t55 * qJD(4) + t6 * t31 - t3) * t40 + t57; 0, 0, 0, 0, 0, 0, 0, t30, t64 * t48 + t47, t54 * t37 + (-t64 * t4 - t3) * t40 + t57;];
tau_reg  = t1;
