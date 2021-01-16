% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:52:19
% EndTime: 2021-01-15 20:52:22
% DurationCPUTime: 0.42s
% Computational Cost: add. (490->88), mult. (1091->146), div. (0->0), fcn. (862->4), ass. (0->48)
t49 = cos(qJ(2));
t41 = t49 * qJD(2);
t47 = sin(qJ(4));
t48 = sin(qJ(2));
t64 = cos(qJ(4));
t56 = t48 * t64;
t70 = qJD(2) * t56 - t47 * t41;
t63 = t48 * t47;
t65 = pkin(6) - pkin(7);
t69 = t65 * t63;
t68 = -t49 * pkin(2) - t48 * qJ(3);
t67 = 2 * qJD(3);
t66 = -pkin(2) - pkin(3);
t62 = qJ(3) * t41 + t48 * qJD(3);
t42 = qJD(4) * t47;
t61 = t48 * qJD(2);
t60 = -0.2e1 * pkin(1) * qJD(2);
t28 = -pkin(1) + t68;
t59 = pkin(6) * t61;
t58 = pkin(6) * t41;
t43 = qJD(4) * t64;
t22 = t49 * pkin(3) - t28;
t55 = t65 * t64;
t54 = t64 * t66;
t52 = t48 * t55;
t51 = t49 * t64 + t63;
t29 = t65 * t49;
t3 = -qJD(4) * t52 + t29 * t42 + t70 * t65;
t12 = t66 * t61 + t62;
t27 = t64 * qJ(3) + t47 * t66;
t50 = t68 * qJD(2) + t49 * qJD(3);
t4 = t29 * t43 - t55 * t41 + (-qJD(2) + qJD(4)) * t69;
t26 = -t47 * qJ(3) - pkin(4) + t54;
t25 = -t49 * t47 + t56;
t17 = pkin(2) * t61 - t62;
t16 = t47 * qJD(3) + t27 * qJD(4);
t15 = qJ(3) * t42 - t64 * qJD(3) - qJD(4) * t54;
t14 = 0.2e1 * t16;
t13 = 0.2e1 * t15;
t11 = pkin(4) * t51 + t22;
t10 = t51 * qJD(2) - t48 * t42 - t49 * t43;
t9 = -t49 * t42 + t48 * t43 - t70;
t7 = -qJ(5) * t51 + t64 * t29 + t69;
t6 = -t25 * qJ(5) - t47 * t29 + t52;
t5 = t9 * pkin(4) + t12;
t2 = t10 * qJ(5) + t25 * qJD(5) + t4;
t1 = t9 * qJ(5) + qJD(5) * t51 + t3;
t8 = [0, 0, 0, 0.2e1 * t48 * t41, 0.2e1 * (-t48 ^ 2 + t49 ^ 2) * qJD(2), 0, 0, 0, t48 * t60, t49 * t60, -0.2e1 * t17 * t49 + 0.2e1 * t28 * t61, 0, -0.2e1 * t17 * t48 - 0.2e1 * t28 * t41, 0.2e1 * t28 * t17, 0.2e1 * t25 * t10, -0.2e1 * t10 * t51 - 0.2e1 * t25 * t9, 0, 0, 0, 0.2e1 * t12 * t51 + 0.2e1 * t22 * t9, 0.2e1 * t22 * t10 + 0.2e1 * t12 * t25, 0.2e1 * t11 * t9 + 0.2e1 * t5 * t51, 0.2e1 * t11 * t10 + 0.2e1 * t5 * t25, 0.2e1 * t1 * t51 - 0.2e1 * t6 * t10 + 0.2e1 * t2 * t25 - 0.2e1 * t7 * t9, -0.2e1 * t7 * t1 + 0.2e1 * t11 * t5 - 0.2e1 * t6 * t2; 0, 0, 0, 0, 0, t41, -t61, 0, -t58, t59, -t58, t50, -t59, t50 * pkin(6), 0, 0, -t10, t9, 0, t4, -t3, t2, -t1, -t26 * t10 + t15 * t51 + t16 * t25 - t27 * t9, -t1 * t27 - t7 * t15 - t6 * t16 - t2 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, qJ(3) * t67, 0, 0, 0, 0, 0, t14, -t13, t14, -t13, 0, -0.2e1 * t27 * t15 - 0.2e1 * t26 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64 * t10 - t47 * t9 + (t25 * t47 - t51 * t64) * qJD(4), -t2 * t64 - t1 * t47 + (-t47 * t6 + t64 * t7) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t43, t42, t43, 0, -t16 * t64 - t15 * t47 + (-t26 * t47 + t64 * t27) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, -t4, t3, -t2, t1, -t10 * pkin(4), -t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t15, -t16, t15, 0, -t16 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t43, -t42, -t43, 0, -pkin(4) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
