% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:18
% EndTime: 2019-07-18 13:30:19
% DurationCPUTime: 0.18s
% Computational Cost: add. (112->36), mult. (378->61), div. (0->0), fcn. (236->6), ass. (0->32)
t21 = sin(qJ(4));
t22 = sin(qJ(3));
t25 = cos(qJ(3));
t24 = cos(qJ(4));
t35 = qJD(4) * t24;
t38 = t22 * t24;
t40 = ((t21 * t25 + t38) * qJD(3) + t22 * t35) * pkin(2);
t18 = t25 * pkin(2) + pkin(3);
t32 = t21 * t22 * pkin(2);
t11 = -t24 * t18 + t32;
t23 = cos(qJ(5));
t19 = qJD(5) * t23;
t20 = sin(qJ(5));
t36 = qJD(4) * t21;
t5 = t18 * t36 + t40;
t1 = t11 * t19 + t5 * t20;
t37 = pkin(2) * qJD(3);
t34 = qJD(5) * t20;
t33 = qJD(5) * t24;
t31 = t22 * t37;
t30 = t25 * t37;
t29 = pkin(3) * t36;
t28 = pkin(3) * t35;
t2 = t11 * t34 - t5 * t23;
t9 = -pkin(3) * t23 * t33 + t20 * t29;
t4 = -t18 * t35 - t24 * t30 + (qJD(3) + qJD(4)) * t32;
t10 = (-t20 * t33 - t23 * t36) * pkin(3);
t17 = t21 * pkin(3) + pkin(6);
t16 = 0.2e1 * t20 * t19;
t12 = 0.2e1 * (-t20 ^ 2 + t23 ^ 2) * qJD(5);
t8 = pkin(2) * t38 + t21 * t18 + pkin(6);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -0.2e1 * t31, -0.2e1 * t30, 0, -0.2e1 * t5, 0.2e1 * t4, t16, t12, 0, 0, 0, 0.2e1 * t2, 0.2e1 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t31, -t30, 0, (-pkin(3) - t18) * t36 - t40, t4 - t28, t16, t12, 0, 0, 0, t10 + t2, t9 + t1; 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t29, -0.2e1 * t28, t16, t12, 0, 0, 0, 0.2e1 * t10, 0.2e1 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t5, t4, t16, t12, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t28, t16, t12, 0, 0, 0, t10, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t12, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t34, 0, -t8 * t19 + t20 * t4, t23 * t4 + t8 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t34, 0, -t17 * t19 - t20 * t28, t17 * t34 - t23 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t34, 0, -pkin(6) * t19, pkin(6) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
