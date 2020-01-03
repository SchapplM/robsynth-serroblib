% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:54
% EndTime: 2019-12-31 17:04:55
% DurationCPUTime: 0.21s
% Computational Cost: add. (280->47), mult. (688->107), div. (0->0), fcn. (611->6), ass. (0->42)
t37 = sin(pkin(7));
t50 = pkin(2) * t37;
t49 = -qJ(3) - pkin(5);
t40 = sin(qJ(2));
t42 = cos(qJ(2));
t44 = qJD(2) * t49;
t23 = t42 * qJD(3) + t40 * t44;
t24 = -t40 * qJD(3) + t42 * t44;
t38 = cos(pkin(7));
t8 = t38 * t23 + t37 * t24;
t32 = t49 * t40;
t33 = t49 * t42;
t14 = t37 * t32 - t38 * t33;
t48 = qJD(2) * t40;
t47 = qJD(2) * t42;
t46 = -0.2e1 * pkin(1) * qJD(2);
t36 = pkin(2) * t48;
t45 = -t42 * pkin(2) - pkin(1);
t7 = -t37 * t23 + t38 * t24;
t13 = t38 * t32 + t37 * t33;
t27 = -t37 * t40 + t38 * t42;
t28 = t37 * t42 + t38 * t40;
t39 = sin(qJ(4));
t41 = cos(qJ(4));
t43 = t41 * t27 - t39 * t28;
t12 = t39 * t27 + t41 * t28;
t35 = t38 * pkin(2) + pkin(3);
t26 = -t37 * t48 + t38 * t47;
t25 = t28 * qJD(2);
t21 = (-t35 * t39 - t41 * t50) * qJD(4);
t20 = (-t35 * t41 + t39 * t50) * qJD(4);
t16 = -t27 * pkin(3) + t45;
t15 = t25 * pkin(3) + t36;
t10 = t27 * pkin(6) + t14;
t9 = -t28 * pkin(6) + t13;
t6 = -t25 * pkin(6) + t8;
t5 = -t26 * pkin(6) + t7;
t4 = t12 * qJD(4) + t41 * t25 + t39 * t26;
t3 = t43 * qJD(4) - t39 * t25 + t41 * t26;
t2 = -t39 * t6 + t41 * t5 + (-t10 * t41 - t39 * t9) * qJD(4);
t1 = -t39 * t5 - t41 * t6 + (t10 * t39 - t41 * t9) * qJD(4);
t11 = [0, 0, 0, 0.2e1 * t40 * t47, 0.2e1 * (-t40 ^ 2 + t42 ^ 2) * qJD(2), 0, 0, 0, t40 * t46, t42 * t46, -0.2e1 * t13 * t26 - 0.2e1 * t14 * t25 + 0.2e1 * t8 * t27 - 0.2e1 * t7 * t28, 0.2e1 * t13 * t7 + 0.2e1 * t14 * t8 + 0.2e1 * t45 * t36, 0.2e1 * t12 * t3, -0.2e1 * t12 * t4 + 0.2e1 * t3 * t43, 0, 0, 0, -0.2e1 * t15 * t43 + 0.2e1 * t16 * t4, 0.2e1 * t15 * t12 + 0.2e1 * t16 * t3; 0, 0, 0, 0, 0, t47, -t48, 0, -pkin(5) * t47, pkin(5) * t48, (-t25 * t37 - t26 * t38) * pkin(2), (t37 * t8 + t38 * t7) * pkin(2), 0, 0, t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t21, 0.2e1 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
