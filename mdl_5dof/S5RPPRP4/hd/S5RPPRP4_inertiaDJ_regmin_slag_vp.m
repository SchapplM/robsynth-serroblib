% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 17:13:14
% EndTime: 2021-01-15 17:13:16
% DurationCPUTime: 0.18s
% Computational Cost: add. (162->40), mult. (280->73), div. (0->0), fcn. (207->4), ass. (0->29)
t20 = cos(qJ(4));
t19 = sin(qJ(4));
t29 = t19 * qJD(4);
t18 = cos(pkin(7));
t31 = qJD(2) * t18;
t17 = sin(pkin(7));
t21 = -pkin(1) - pkin(2);
t32 = t18 * qJ(2) + t17 * t21;
t7 = -pkin(6) + t32;
t24 = -t20 * t31 + t7 * t29;
t1 = -qJ(5) * t29 + t20 * qJD(5) + t24;
t28 = t20 * qJD(4);
t33 = qJ(5) - t7;
t2 = (qJD(5) - t31) * t19 + t33 * t28;
t3 = t33 * t19;
t4 = t33 * t20;
t35 = (-t19 * t4 + t20 * t3) * qJD(4) + t1 * t20 + t2 * t19;
t34 = 0.2e1 * qJD(2);
t30 = t17 * qJD(2);
t27 = pkin(4) * t29;
t26 = t17 * t28;
t25 = -t17 * qJ(2) + t18 * t21;
t6 = pkin(3) - t25;
t12 = t18 * t28;
t10 = t18 * t29;
t9 = t17 * t29;
t8 = -t27 + t30;
t5 = t20 * pkin(4) + t6;
t11 = [0, 0, 0, 0, t34, qJ(2) * t34, (-t25 * t17 + t32 * t18) * t34, 0.2e1 * t19 * t28, 0.2e1 * (-t19 ^ 2 + t20 ^ 2) * qJD(4), 0, 0, 0, 0.2e1 * t20 * t30 - 0.2e1 * t6 * t29, -0.2e1 * t19 * t30 - 0.2e1 * t6 * t28, 0.2e1 * t8 * t20 - 0.2e1 * t5 * t29, -0.2e1 * t8 * t19 - 0.2e1 * t5 * t28, 0.2e1 * t35, 0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t5 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t12, t10, t12, 0, -t35 * t17 - t8 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t19 + t2 * t20 + (-t19 * t3 - t20 * t4) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t29, 0, -t19 * t31 - t7 * t28, t24, t2, t1, pkin(4) * t28, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t9, -t26, t9, 0, -pkin(4) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t28, -t29, -t28, 0, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t28, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
