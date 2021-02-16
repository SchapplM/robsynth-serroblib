% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:48:16
% EndTime: 2021-01-15 14:48:17
% DurationCPUTime: 0.17s
% Computational Cost: add. (111->36), mult. (326->71), div. (0->0), fcn. (267->6), ass. (0->33)
t19 = sin(pkin(8));
t21 = sin(qJ(3));
t30 = cos(pkin(8));
t34 = cos(qJ(3));
t9 = t21 * t19 - t34 * t30;
t36 = 2 * qJD(4);
t22 = cos(qJ(4));
t35 = t22 * pkin(4);
t31 = -qJ(5) - pkin(6);
t12 = t31 * t22;
t33 = t12 * t22;
t20 = sin(qJ(4));
t29 = t20 * qJD(4);
t16 = t22 * qJD(4);
t28 = -2 * pkin(3) * qJD(4);
t27 = 0.2e1 * t29;
t26 = pkin(4) * t29;
t17 = t20 ^ 2;
t18 = t22 ^ 2;
t7 = t9 * qJD(3);
t25 = (-t17 - t18) * t7;
t10 = t34 * t19 + t21 * t30;
t2 = -t10 * t16 + t20 * t7;
t11 = t31 * t20;
t5 = -t22 * qJD(5) - t31 * t29;
t6 = -t20 * qJD(5) + t31 * t16;
t23 = -t6 * t20 - t5 * t22 + (-t11 * t22 + t12 * t20) * qJD(4);
t14 = -pkin(3) - t35;
t8 = t10 * qJD(3);
t4 = -t8 * t22 + t9 * t29;
t3 = t9 * t16 + t8 * t20;
t1 = t10 * t29 + t22 * t7;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t10 * t25 + 0.2e1 * t9 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, t4, t3, t4, t3, t25, t7 * t33 + t8 * t14 + (pkin(4) * qJD(4) * t9 + t11 * t7) * t20 + t23 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20 * t5 + t22 * t6 + (-t11 * t20 - t33) * qJD(4); 0, 0, 0, 0, 0, t22 * t27, (-t17 + t18) * t36, 0, 0, 0, t20 * t28, t22 * t28, (t14 - t35) * t27, (pkin(4) * t17 + t14 * t22) * t36, 0.2e1 * t23, 0.2e1 * t11 * t6 + 0.2e1 * t12 * t5 + 0.2e1 * t14 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, t2, t1, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t16, -t29, -t16, 0, -t26; 0, 0, 0, 0, 0, 0, 0, t16, -t29, 0, -pkin(6) * t16, pkin(6) * t29, t6, t5, -pkin(4) * t16, t6 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t16, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t13;
