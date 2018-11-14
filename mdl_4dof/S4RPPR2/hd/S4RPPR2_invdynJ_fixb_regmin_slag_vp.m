% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPPR2
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
% tau_reg [4x12]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:48
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4RPPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:47:33
% EndTime: 2018-11-14 13:47:33
% DurationCPUTime: 0.25s
% Computational Cost: add. (274->88), mult. (393->107), div. (0->0), fcn. (238->8), ass. (0->57)
t67 = qJD(1) - qJD(4);
t68 = qJD(4) + t67;
t37 = -pkin(1) - pkin(2);
t17 = t37 * qJDD(1) + qJDD(2);
t31 = sin(pkin(6));
t32 = cos(pkin(6));
t60 = qJ(2) * qJDD(1);
t66 = t31 * t17 + t32 * t60;
t34 = sin(qJ(1));
t36 = cos(qJ(1));
t65 = t36 * pkin(1) + t34 * qJ(2);
t64 = g(1) * t34 - g(2) * t36;
t63 = qJ(2) * t31;
t62 = pkin(1) * qJDD(1);
t61 = qJ(2) * qJD(1);
t59 = qJD(1) * qJD(2);
t58 = pkin(6) + qJ(4);
t57 = 0.2e1 * t59;
t56 = t31 * t59;
t55 = t32 * t59;
t54 = -pkin(3) - t63;
t14 = t32 * t37 - t63;
t53 = cos(t58);
t52 = sin(t58);
t51 = qJDD(2) - t62;
t50 = g(1) * t36 + g(2) * t34;
t21 = t37 * qJD(1) + qJD(2);
t16 = t32 * t21;
t6 = t21 * t31 + t32 * t61;
t49 = t31 * (-t31 * t61 + t16) - t32 * t6;
t33 = sin(qJ(4));
t35 = cos(qJ(4));
t4 = t54 * qJD(1) + t16;
t48 = t33 * t6 - t35 * t4;
t47 = -t33 * t4 - t35 * t6;
t11 = -pkin(3) + t14;
t15 = qJ(2) * t32 + t31 * t37;
t46 = t11 * t35 - t15 * t33;
t45 = t11 * t33 + t15 * t35;
t44 = t31 * t35 + t32 * t33;
t43 = t31 * t33 - t32 * t35;
t42 = t67 * t44;
t41 = t43 * t67;
t13 = t32 * t17;
t1 = t54 * qJDD(1) + t13 - t56;
t3 = t55 + t66;
t7 = -t34 * t52 - t36 * t53;
t8 = -t34 * t53 + t36 * t52;
t40 = g(1) * t8 - g(2) * t7 + t35 * t1 - t33 * t3;
t39 = -g(1) * t7 - g(2) * t8 - t33 * t1 - t35 * t3;
t38 = qJD(1) ^ 2;
t29 = qJDD(1) - qJDD(4);
t25 = t36 * qJ(2);
t10 = t31 * t34 + t32 * t36;
t9 = t31 * t36 - t32 * t34;
t2 = t13 + (-t59 - t60) * t31;
t5 = [qJDD(1), t64, t50, -qJDD(2) + 0.2e1 * t62 + t64, -t50 + t57 + 0.2e1 * t60, -t51 * pkin(1) - g(1) * (-pkin(1) * t34 + t25) - g(2) * t65 + (t57 + t60) * qJ(2), 0.2e1 * t56 - g(1) * t9 - g(2) * t10 - t13 + (-t14 + t63) * qJDD(1), -g(1) * t10 + g(2) * t9 + qJDD(1) * t15 + 0.2e1 * t55 + t66, t3 * t15 + t2 * t14 - g(1) * (t37 * t34 + t25) - g(2) * (pkin(2) * t36 + t65) - t49 * qJD(2), t29, -t46 * t29 + qJD(2) * t42 + (t45 * t67 - t47) * qJD(4) - t40, t45 * t29 - qJD(2) * t41 + (t46 * t67 - t48) * qJD(4) - t39; 0, 0, 0, -qJDD(1), -t38, -qJ(2) * t38 + t51 - t64, -qJDD(1) * t32 - t31 * t38, qJDD(1) * t31 - t32 * t38, t49 * qJD(1) + t2 * t32 + t3 * t31 - t64, 0, t43 * t29 - t67 * t42, t44 * t29 + t67 * t41; 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) + g(3), 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, t68 * t47 + t40, t68 * t48 + t39;];
tau_reg  = t5;
