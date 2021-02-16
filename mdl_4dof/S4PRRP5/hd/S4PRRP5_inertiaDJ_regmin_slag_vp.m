% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x15]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 22:36
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRRP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:36:08
% EndTime: 2021-01-14 22:36:09
% DurationCPUTime: 0.13s
% Computational Cost: add. (62->27), mult. (206->61), div. (0->0), fcn. (137->4), ass. (0->29)
t29 = 2 * qJD(3);
t16 = cos(qJ(3));
t28 = t16 * pkin(3);
t27 = qJ(4) + pkin(5);
t14 = sin(qJ(3));
t12 = t14 ^ 2;
t13 = t16 ^ 2;
t26 = t12 + t13;
t15 = sin(qJ(2));
t25 = qJD(2) * t15;
t17 = cos(qJ(2));
t24 = qJD(2) * t17;
t23 = qJD(3) * t14;
t11 = qJD(3) * t16;
t22 = qJD(3) * t17;
t21 = -2 * pkin(2) * qJD(3);
t20 = 0.2e1 * t23;
t19 = pkin(3) * t23;
t4 = -t15 * t11 - t14 * t24;
t1 = -t16 * qJD(4) + t27 * t23;
t2 = -t14 * qJD(4) - t27 * t11;
t7 = t27 * t14;
t8 = t27 * t16;
t18 = -t1 * t16 - t2 * t14 + (-t14 * t8 + t16 * t7) * qJD(3);
t9 = -pkin(2) - t28;
t6 = -t14 * t22 - t16 * t25;
t5 = t14 * t25 - t16 * t22;
t3 = t15 * t23 - t16 * t24;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t26) * t15 * t24; 0, 0, -t25, -t24, 0, 0, 0, 0, 0, t6, t5, t6, t5, t26 * t24, (-t19 + (t14 * t7 + t16 * t8) * qJD(2)) * t17 + (qJD(2) * t9 + t18) * t15; 0, 0, 0, 0, t16 * t20, (-t12 + t13) * t29, 0, 0, 0, t14 * t21, t16 * t21, (t9 - t28) * t20, (pkin(3) * t12 + t16 * t9) * t29, 0.2e1 * t18, -0.2e1 * t8 * t1 + 0.2e1 * t9 * t19 - 0.2e1 * t7 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, t4, t3, 0, t4 * pkin(3); 0, 0, 0, 0, 0, 0, t11, -t23, 0, -pkin(5) * t11, pkin(5) * t23, t2, t1, -pkin(3) * t11, t2 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t11, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
