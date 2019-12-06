% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x13]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPPRR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:37
% EndTime: 2019-12-05 14:59:38
% DurationCPUTime: 0.13s
% Computational Cost: add. (32->24), mult. (131->56), div. (0->0), fcn. (124->8), ass. (0->25)
t5 = sin(pkin(9));
t6 = sin(pkin(8));
t24 = t5 * t6;
t7 = cos(pkin(9));
t23 = t6 * t7;
t12 = cos(qJ(4));
t22 = t12 * t5;
t9 = sin(qJ(5));
t21 = t9 * qJD(5);
t10 = sin(qJ(4));
t20 = t10 * qJD(4);
t11 = cos(qJ(5));
t19 = t11 * qJD(5);
t18 = t12 * qJD(4);
t17 = -0.2e1 * pkin(4) * qJD(5);
t16 = t5 * t20;
t15 = t11 * t20;
t8 = cos(pkin(8));
t4 = -t8 * t10 + t12 * t23;
t3 = t10 * t23 + t8 * t12;
t14 = t10 * t19 + t9 * t18;
t13 = t10 * t21 - t11 * t18;
t2 = t3 * qJD(4);
t1 = t4 * qJD(4);
t25 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t1, t2, 0, 0, 0, 0, 0, -t1 * t11 + t3 * t21, t1 * t9 + t3 * t19; 0, 0, 0, 0, -t5 * t18, t16, 0, 0, 0, 0, 0, t13 * t5, t14 * t5; 0, 0, 0, 0, -t20, -t18, 0, 0, 0, 0, 0, -t12 * t21 - t15, -t12 * t19 + t9 * t20; 0, 0, 0, 0, 0, 0, 0.2e1 * t9 * t19, 0.2e1 * (t11 ^ 2 - t9 ^ 2) * qJD(5), 0, 0, 0, t9 * t17, t11 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 * t9 + (-t11 * t4 - t9 * t24) * qJD(5), t2 * t11 + (-t11 * t24 + t4 * t9) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9 * t16 + (-t11 * t22 + t7 * t9) * qJD(5), t5 * t15 + (t11 * t7 + t9 * t22) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t13; 0, 0, 0, 0, 0, 0, 0, 0, t19, -t21, 0, -pkin(6) * t19, pkin(6) * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t25;
