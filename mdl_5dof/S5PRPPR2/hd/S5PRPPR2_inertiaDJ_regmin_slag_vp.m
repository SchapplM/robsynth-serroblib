% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPPR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:24:41
% EndTime: 2019-12-05 15:24:42
% DurationCPUTime: 0.16s
% Computational Cost: add. (98->31), mult. (296->65), div. (0->0), fcn. (275->8), ass. (0->29)
t16 = sin(pkin(9));
t18 = cos(pkin(9));
t20 = sin(qJ(5));
t22 = cos(qJ(5));
t8 = t20 * t16 - t22 * t18;
t3 = t8 * qJD(5);
t33 = -0.2e1 * t3;
t17 = sin(pkin(8));
t19 = cos(pkin(8));
t21 = sin(qJ(2));
t23 = cos(qJ(2));
t9 = t17 * t23 + t19 * t21;
t1 = t9 * qJD(2);
t7 = t17 * t21 - t19 * t23;
t32 = t7 * t1;
t13 = t17 * pkin(2) + qJ(4);
t31 = pkin(6) + t13;
t28 = t16 ^ 2 + t18 ^ 2;
t27 = -t19 * pkin(2) - pkin(3);
t2 = t7 * qJD(2);
t26 = t28 * t2;
t25 = t28 * qJD(4);
t24 = 0.2e1 * t25;
t10 = t22 * t16 + t20 * t18;
t4 = t10 * qJD(5);
t11 = -t18 * pkin(4) + t27;
t6 = t31 * t18;
t5 = t31 * t16;
t12 = [0, 0, 0, 0, -0.2e1 * t9 * t2 + 0.2e1 * t32, 0, 0, 0, -0.2e1 * t9 * t26 + 0.2e1 * t32, 0, 0, 0, 0, 0, 0, 0; 0, 0, -qJD(2) * t21, -qJD(2) * t23, (-t1 * t19 - t17 * t2) * pkin(2), -t1 * t18, t1 * t16, -t26, t1 * t27 - t13 * t26 + t9 * t25, 0, 0, 0, 0, 0, t1 * t8 + t7 * t4, t1 * t10 - t7 * t3; 0, 0, 0, 0, 0, 0, 0, t24, t13 * t24, t10 * t33, -0.2e1 * t10 * t4 + 0.2e1 * t3 * t8, 0, 0, 0, 0.2e1 * t11 * t4, t11 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t2 + t9 * t3, -t8 * t2 + t9 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, (t20 * t5 - t22 * t6) * qJD(5) - t10 * qJD(4), (t20 * t6 + t22 * t5) * qJD(5) + t8 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
