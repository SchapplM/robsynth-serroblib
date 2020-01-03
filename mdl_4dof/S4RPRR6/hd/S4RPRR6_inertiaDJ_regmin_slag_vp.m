% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:41
% EndTime: 2019-12-31 16:52:42
% DurationCPUTime: 0.19s
% Computational Cost: add. (240->38), mult. (594->76), div. (0->0), fcn. (574->6), ass. (0->36)
t22 = sin(pkin(7));
t23 = cos(pkin(7));
t25 = sin(qJ(3));
t27 = cos(qJ(3));
t29 = t25 * t22 - t27 * t23;
t37 = pkin(5) + qJ(2);
t17 = t37 * t22;
t18 = t37 * t23;
t31 = -t27 * t17 - t25 * t18;
t40 = t29 * qJD(2) - t31 * qJD(3);
t19 = -t23 * pkin(2) - pkin(1);
t39 = 0.2e1 * t19;
t16 = t27 * t22 + t25 * t23;
t14 = t16 * qJD(3);
t38 = pkin(3) * t14;
t36 = pkin(3) * qJD(4);
t24 = sin(qJ(4));
t35 = t24 * t36;
t26 = cos(qJ(4));
t34 = t26 * t36;
t33 = 0.2e1 * (t22 ^ 2 + t23 ^ 2) * qJD(2);
t32 = -t24 * t16 - t26 * t29;
t10 = t26 * t16 - t24 * t29;
t30 = t25 * t17 - t27 * t18;
t28 = -t16 * qJD(2) + t30 * qJD(3);
t13 = t29 * qJD(3);
t11 = pkin(3) * t29 + t19;
t8 = -pkin(6) * t29 - t30;
t7 = -t16 * pkin(6) + t31;
t6 = t13 * pkin(6) + t28;
t5 = -t14 * pkin(6) - t40;
t4 = t10 * qJD(4) - t24 * t13 + t26 * t14;
t3 = t32 * qJD(4) - t26 * t13 - t24 * t14;
t2 = -t24 * t5 + t26 * t6 + (-t24 * t7 - t26 * t8) * qJD(4);
t1 = -t24 * t6 - t26 * t5 + (t24 * t8 - t26 * t7) * qJD(4);
t9 = [0, 0, 0, 0, 0, t33, qJ(2) * t33, -0.2e1 * t16 * t13, 0.2e1 * t13 * t29 - 0.2e1 * t16 * t14, 0, 0, 0, t14 * t39, -t13 * t39, 0.2e1 * t10 * t3, -0.2e1 * t10 * t4 + 0.2e1 * t3 * t32, 0, 0, 0, 0.2e1 * t11 * t4 - 0.2e1 * t32 * t38, 0.2e1 * t10 * t38 + 0.2e1 * t11 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, 0, t28, t40, 0, 0, t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t35, -0.2e1 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
