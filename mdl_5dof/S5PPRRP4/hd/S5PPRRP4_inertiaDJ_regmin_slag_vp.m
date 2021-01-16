% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:56
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:56:28
% EndTime: 2021-01-15 14:56:29
% DurationCPUTime: 0.15s
% Computational Cost: add. (71->29), mult. (232->64), div. (0->0), fcn. (156->4), ass. (0->30)
t30 = 2 * qJD(4);
t18 = cos(qJ(4));
t29 = t18 * pkin(4);
t28 = qJ(5) + pkin(6);
t16 = sin(qJ(4));
t14 = t16 ^ 2;
t15 = t18 ^ 2;
t27 = t14 + t15;
t17 = sin(qJ(3));
t26 = qJD(3) * t17;
t19 = cos(qJ(3));
t25 = qJD(3) * t19;
t12 = qJD(4) * t16;
t13 = qJD(4) * t18;
t24 = qJD(4) * t19;
t23 = -2 * pkin(3) * qJD(4);
t22 = 0.2e1 * t12;
t10 = pkin(4) * t12;
t7 = t28 * t16;
t8 = t28 * t18;
t21 = -t16 * t7 - t18 * t8;
t4 = -t17 * t13 - t16 * t25;
t1 = -t18 * qJD(5) + t28 * t12;
t2 = -t16 * qJD(5) - t28 * t13;
t20 = -t1 * t18 - t2 * t16 + (-t16 * t8 + t18 * t7) * qJD(4);
t9 = -pkin(3) - t29;
t6 = -t16 * t24 - t18 * t26;
t5 = t16 * t26 - t18 * t24;
t3 = t17 * t12 - t18 * t25;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t27) * t17 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * qJD(4) + t16 * t1 - t18 * t2; 0, 0, 0, -t26, -t25, 0, 0, 0, 0, 0, t6, t5, t6, t5, t27 * t25, (-t21 * qJD(3) - t10) * t19 + (qJD(3) * t9 + t20) * t17; 0, 0, 0, 0, 0, t18 * t22, (-t14 + t15) * t30, 0, 0, 0, t16 * t23, t18 * t23, (t9 - t29) * t22, (pkin(4) * t14 + t18 * t9) * t30, 0.2e1 * t20, -0.2e1 * t8 * t1 + 0.2e1 * t9 * t10 - 0.2e1 * t7 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t13, t12, t13, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, t4, t3, 0, t4 * pkin(4); 0, 0, 0, 0, 0, 0, 0, t13, -t12, 0, -pkin(6) * t13, pkin(6) * t12, t2, t1, -pkin(4) * t13, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t13, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
