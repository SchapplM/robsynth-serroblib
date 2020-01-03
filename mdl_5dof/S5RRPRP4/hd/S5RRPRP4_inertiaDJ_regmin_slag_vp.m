% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRP4
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
% MMD_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:56
% EndTime: 2019-12-31 19:52:57
% DurationCPUTime: 0.28s
% Computational Cost: add. (182->64), mult. (402->90), div. (0->0), fcn. (220->4), ass. (0->45)
t34 = sin(qJ(4));
t36 = cos(qJ(4));
t58 = -t34 * pkin(4) + t36 * qJ(5);
t39 = 2 * qJD(3);
t57 = 2 * qJD(5);
t18 = qJ(3) - t58;
t35 = sin(qJ(2));
t55 = t35 * pkin(1);
t13 = t18 + t55;
t37 = cos(qJ(2));
t51 = pkin(1) * qJD(2);
t44 = t37 * t51;
t45 = qJD(4) * qJ(5);
t47 = t36 * qJD(4);
t8 = pkin(4) * t47 - t36 * qJD(5) + t34 * t45 + qJD(3);
t2 = t8 + t44;
t56 = t13 * t47 + t2 * t34;
t54 = t18 * t47 + t8 * t34;
t21 = qJD(3) + t44;
t23 = qJ(3) + t55;
t53 = t21 * t34 + t23 * t47;
t46 = qJ(3) * qJD(4);
t52 = qJD(3) * t34 + t36 * t46;
t49 = t34 * qJD(4);
t48 = t34 * qJD(5);
t28 = t35 * t51;
t38 = -pkin(2) - pkin(7);
t43 = t38 * t49;
t42 = t38 * t47;
t41 = -t37 * pkin(1) - pkin(2);
t32 = t34 ^ 2;
t33 = t36 ^ 2;
t9 = (t32 + t33) * t28;
t40 = t58 * qJD(4) + t48;
t30 = qJD(3) * t36;
t22 = -pkin(7) + t41;
t20 = -0.2e1 * t34 * t47;
t16 = t21 * t36;
t14 = 0.2e1 * (t32 - t33) * qJD(4);
t11 = t18 * t49;
t10 = -pkin(4) * t49 + t36 * t45 + t48;
t6 = t13 * t49;
t5 = t22 * t47 + t34 * t28;
t4 = t22 * t49 - t36 * t28;
t1 = [0, 0, 0, 0, -0.2e1 * t28, -0.2e1 * t44, 0.2e1 * t28, 0.2e1 * t21, 0.2e1 * t23 * t21 + 0.2e1 * t41 * t28, t20, t14, 0, 0, 0, 0.2e1 * t53, -0.2e1 * t23 * t49 + 0.2e1 * t16, 0.2e1 * t56, -0.2e1 * t9, -0.2e1 * t2 * t36 + 0.2e1 * t6, 0.2e1 * t13 * t2 + 0.2e1 * t22 * t9; 0, 0, 0, 0, -t28, -t44, t28, t39 + t44, -pkin(2) * t28 + t21 * qJ(3) + t23 * qJD(3), t20, t14, 0, 0, 0, t52 + t53, t16 + t30 + (-qJ(3) - t23) * t49, t54 + t56, -t9, t11 + t6 + (-t2 - t8) * t36, t13 * t8 + t2 * t18 + t38 * t9; 0, 0, 0, 0, 0, 0, 0, t39, qJ(3) * t39, t20, t14, 0, 0, 0, 0.2e1 * t52, -0.2e1 * t34 * t46 + 0.2e1 * t30, 0.2e1 * t54, 0, -0.2e1 * t8 * t36 + 0.2e1 * t11, 0.2e1 * t18 * t8; 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t47, 0, -t4, -t5, -t4, -t10, t5, (pkin(4) * t36 + qJ(5) * t34) * t28 + t40 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t47, 0, -t43, -t42, -t43, -t10, t42, t40 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t47, -t49, 0, t47, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, qJ(5) * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, 0, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
