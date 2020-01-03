% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:49
% EndTime: 2019-12-31 18:17:49
% DurationCPUTime: 0.18s
% Computational Cost: add. (104->33), mult. (247->51), div. (0->0), fcn. (154->6), ass. (0->26)
t13 = cos(pkin(8)) * pkin(1) + pkin(2);
t19 = sin(qJ(3));
t21 = cos(qJ(3));
t32 = pkin(1) * sin(pkin(8));
t34 = -t21 * t13 + t19 * t32;
t23 = 2 * qJD(4);
t18 = sin(qJ(5));
t20 = cos(qJ(5));
t28 = qJD(5) * t20;
t5 = t34 * qJD(3);
t4 = -qJD(4) + t5;
t24 = t19 * t13 + t21 * t32;
t8 = qJ(4) + t24;
t33 = -t4 * t18 + t8 * t28;
t26 = qJ(4) * qJD(5);
t30 = qJD(4) * t18 + t20 * t26;
t29 = qJD(5) * t18;
t27 = qJD(5) * (-pkin(3) - pkin(7));
t25 = -pkin(3) + t34;
t16 = qJD(4) * t20;
t11 = -0.2e1 * t18 * t28;
t9 = 0.2e1 * (t18 ^ 2 - t20 ^ 2) * qJD(5);
t7 = -pkin(7) + t25;
t6 = t24 * qJD(3);
t2 = t4 * t20;
t1 = [0, 0, 0, 0, 0, -0.2e1 * t6, 0.2e1 * t5, 0.2e1 * t6, -0.2e1 * t4, 0.2e1 * t25 * t6 - 0.2e1 * t8 * t4, t11, t9, 0, 0, 0, 0.2e1 * t33, -0.2e1 * t8 * t29 - 0.2e1 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t6, t5, t6, t23 - t5, -t6 * pkin(3) - t4 * qJ(4) + t8 * qJD(4), t11, t9, 0, 0, 0, t30 + t33, t16 - t2 + (-qJ(4) - t8) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t23, qJ(4) * t23, t11, t9, 0, 0, 0, 0.2e1 * t30, -0.2e1 * t18 * t26 + 0.2e1 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t28, 0, t20 * t6 - t7 * t29, -t18 * t6 - t7 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t28, 0, -t18 * t27, -t20 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
