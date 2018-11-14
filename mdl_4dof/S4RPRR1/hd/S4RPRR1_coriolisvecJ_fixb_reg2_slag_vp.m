% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:51
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:50:35
% EndTime: 2018-11-14 13:50:35
% DurationCPUTime: 0.19s
% Computational Cost: add. (379->50), mult. (964->71), div. (0->0), fcn. (546->6), ass. (0->36)
t24 = cos(pkin(7)) * pkin(1) + pkin(2);
t22 = t24 * qJD(1);
t30 = sin(qJ(3));
t32 = cos(qJ(3));
t46 = sin(pkin(7)) * pkin(1);
t42 = qJD(1) * t46;
t14 = t32 * t22 - t30 * t42;
t36 = t32 * t24 - t30 * t46;
t26 = qJD(1) + qJD(3);
t10 = t26 * pkin(3) + t14;
t12 = t14 * qJD(3);
t31 = cos(qJ(4));
t15 = t30 * t22 + t32 * t42;
t13 = t15 * qJD(3);
t29 = sin(qJ(4));
t45 = t29 * t15;
t39 = -qJD(4) * t45 - t29 * t13;
t1 = (qJD(4) * t10 + t12) * t31 + t39;
t44 = t31 * t15;
t25 = qJD(4) + t26;
t41 = -pkin(3) * t25 - t10;
t40 = -t29 * t12 - t31 * t13;
t6 = t29 * t10 + t44;
t18 = pkin(3) + t36;
t19 = t30 * t24 + t32 * t46;
t35 = t31 * t18 - t29 * t19;
t34 = t29 * t18 + t31 * t19;
t2 = -qJD(4) * t6 + t40;
t17 = t19 * qJD(3);
t16 = t36 * qJD(3);
t8 = t31 * t14 - t45;
t7 = -t29 * t14 - t44;
t5 = t31 * t10 - t45;
t4 = -t34 * qJD(4) - t29 * t16 - t31 * t17;
t3 = t35 * qJD(4) + t31 * t16 - t29 * t17;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17 * t26 - t13, -t16 * t26 - t12, 0, t12 * t19 - t13 * t36 - t14 * t17 + t15 * t16, 0, 0, 0, 0, 0, 0, t4 * t25 + t2, -t3 * t25 - t1, 0, t1 * t34 + t2 * t35 + t6 * t3 + t5 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 * t26 - t13, t14 * t26 - t12, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * t25 + (t41 * t29 - t44) * qJD(4) + t40, t8 * t25 + (t41 * qJD(4) - t12) * t31 - t39, 0, -t5 * t7 - t6 * t8 + (t1 * t29 + t2 * t31 + (-t29 * t5 + t31 * t6) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * t25 + t2, t5 * t25 - t1, 0, 0;];
tauc_reg  = t9;
