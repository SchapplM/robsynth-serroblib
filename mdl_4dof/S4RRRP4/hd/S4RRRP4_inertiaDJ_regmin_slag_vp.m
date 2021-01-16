% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:30:23
% EndTime: 2021-01-15 14:30:24
% DurationCPUTime: 0.23s
% Computational Cost: add. (321->53), mult. (796->106), div. (0->0), fcn. (646->4), ass. (0->41)
t49 = qJD(2) + qJD(3);
t48 = pkin(5) + pkin(6);
t47 = cos(qJ(3));
t29 = sin(qJ(3));
t30 = sin(qJ(2));
t46 = t29 * t30;
t31 = cos(qJ(2));
t45 = t29 * t31;
t44 = qJD(3) * t29;
t43 = t30 * qJD(2);
t42 = t31 * qJD(2);
t41 = -0.2e1 * pkin(1) * qJD(2);
t40 = pkin(2) * t43;
t39 = pkin(2) * t44;
t38 = t47 * pkin(2);
t28 = -t31 * pkin(2) - pkin(1);
t37 = t47 * t31;
t36 = qJD(2) * t48;
t35 = t47 * qJD(3);
t34 = pkin(2) * t35;
t33 = qJD(2) * t37;
t20 = t48 * t30;
t21 = t48 * t31;
t32 = t29 * t20 - t47 * t21;
t18 = t47 * t30 + t45;
t19 = t30 * t36;
t3 = t47 * t19 + t20 * t35 + t21 * t44 + t36 * t45;
t4 = t32 * qJD(3) + t29 * t19 - t48 * t33;
t27 = t38 + pkin(3);
t25 = -0.2e1 * t34;
t24 = -0.2e1 * t39;
t17 = -t37 + t46;
t11 = t17 * pkin(3) + t28;
t10 = t49 * t18;
t9 = -t31 * t35 + t49 * t46 - t33;
t7 = t10 * pkin(3) + t40;
t6 = -t17 * qJ(4) - t32;
t5 = -t18 * qJ(4) - t47 * t20 - t29 * t21;
t2 = t9 * qJ(4) - t18 * qJD(4) + t4;
t1 = t10 * qJ(4) + t17 * qJD(4) + t3;
t8 = [0, 0, 0, 0.2e1 * t30 * t42, 0.2e1 * (-t30 ^ 2 + t31 ^ 2) * qJD(2), 0, 0, 0, t30 * t41, t31 * t41, -0.2e1 * t18 * t9, -0.2e1 * t18 * t10 + 0.2e1 * t9 * t17, 0, 0, 0, 0.2e1 * t28 * t10 + 0.2e1 * t17 * t40, 0.2e1 * t18 * t40 - 0.2e1 * t28 * t9, 0.2e1 * t11 * t10 + 0.2e1 * t7 * t17, -0.2e1 * t11 * t9 + 0.2e1 * t7 * t18, 0.2e1 * t1 * t17 - 0.2e1 * t6 * t10 - 0.2e1 * t2 * t18 + 0.2e1 * t5 * t9, -0.2e1 * t6 * t1 + 0.2e1 * t11 * t7 + 0.2e1 * t5 * t2; 0, 0, 0, 0, 0, t42, -t43, 0, -pkin(5) * t42, pkin(5) * t43, 0, 0, -t9, -t10, 0, t4, t3, t2, t1, t27 * t9 + (-t10 * t29 + (-t47 * t17 + t18 * t29) * qJD(3)) * pkin(2), t2 * t27 + (-t1 * t29 + (-t29 * t5 + t47 * t6) * qJD(3)) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t25, t24, t25, 0, 0.2e1 * (t38 - t27) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, 0, t4, t3, t2, t1, t9 * pkin(3), t2 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t34, -t39, -t34, 0, -pkin(3) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
