% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PPPRR1
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPPRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:58:07
% EndTime: 2019-12-05 14:58:08
% DurationCPUTime: 0.13s
% Computational Cost: add. (35->18), mult. (135->43), div. (0->0), fcn. (132->8), ass. (0->20)
t10 = sin(pkin(9));
t12 = cos(pkin(9));
t15 = sin(qJ(4));
t17 = cos(qJ(4));
t7 = t15 * t10 - t17 * t12;
t14 = sin(qJ(5));
t20 = t14 * qJD(5);
t16 = cos(qJ(5));
t19 = t16 * qJD(5);
t18 = -0.2e1 * pkin(4) * qJD(5);
t8 = t17 * t10 + t15 * t12;
t11 = sin(pkin(8));
t4 = t7 * t11;
t6 = t8 * qJD(4);
t13 = cos(pkin(8));
t5 = t7 * qJD(4);
t3 = t8 * t11;
t2 = qJD(4) * t4;
t1 = t11 * t6;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t2, t1, 0, 0, 0, 0, 0, t2 * t16 + t3 * t20, -t2 * t14 + t3 * t19; 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, -t6 * t16 + t7 * t20, t6 * t14 + t7 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0.2e1 * t14 * t19, 0.2e1 * (-t14 ^ 2 + t16 ^ 2) * qJD(5), 0, 0, 0, t14 * t18, t16 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t1 + (t13 * t14 + t16 * t4) * qJD(5), t16 * t1 + (t13 * t16 - t14 * t4) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t5 - t8 * t19, t16 * t5 + t8 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t19; 0, 0, 0, 0, 0, 0, 0, 0, t19, -t20, 0, -pkin(6) * t19, pkin(6) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
