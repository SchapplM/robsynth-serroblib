% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRRP1
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
% taug_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:55
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRP1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:54:31
% EndTime: 2018-11-14 13:54:32
% DurationCPUTime: 0.19s
% Computational Cost: add. (314->48), mult. (774->73), div. (0->0), fcn. (396->4), ass. (0->44)
t27 = qJD(1) + qJD(2);
t31 = cos(qJ(2));
t49 = pkin(1) * qJD(1);
t45 = t31 * t49;
t21 = t27 * pkin(2) + t45;
t30 = cos(qJ(3));
t43 = qJD(2) * t45;
t28 = sin(qJ(3));
t29 = sin(qJ(2));
t46 = t29 * t49;
t44 = t28 * t46;
t50 = (qJD(2) + qJD(3)) * t44;
t10 = (qJD(3) * t21 + t43) * t30 - t50;
t25 = t31 * pkin(1) + pkin(2);
t52 = t28 * t29;
t37 = t30 * t31 - t52;
t47 = qJD(3) * t30;
t48 = qJD(3) * t28;
t12 = t25 * t47 + (t37 * qJD(2) - t29 * t48) * pkin(1);
t17 = t28 * t21 + t30 * t46;
t51 = t29 * t30;
t53 = t10 * (pkin(1) * t51 + t28 * t25) + t17 * t12;
t42 = -pkin(1) * t52 + t30 * t25;
t19 = t37 * t49;
t41 = -t17 * t19 + (t10 * t28 + t17 * t47) * pkin(2);
t40 = (-qJD(2) + t27) * t49;
t39 = pkin(1) * qJD(2) * (-qJD(1) - t27);
t38 = -t28 * t31 - t51;
t26 = qJD(3) + t27;
t36 = (-pkin(2) * t26 - t21) * qJD(3);
t16 = t30 * t21 - t44;
t33 = (t38 * qJD(2) - t29 * t47) * pkin(1);
t32 = qJD(1) * t33;
t11 = -t21 * t48 + t32;
t18 = t38 * t49;
t14 = t26 * pkin(3) + t16;
t13 = -t25 * t48 + t33;
t6 = t17 * t26 + t11;
t5 = t16 * t26 - t10;
t4 = t19 * t26 + (t36 - t43) * t30 + t50;
t3 = -t18 * t26 + t28 * t36 + t32;
t2 = t13 * t26 + t11;
t1 = -t12 * t26 - t10;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 * t39, t31 * t39, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t11 * t42 + t16 * t13 + t53, 0, 0, 0, 0, 0, 0, t2, t1, 0, t11 * (pkin(3) + t42) + t14 * t13 + t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 * t40, t31 * t40, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, -t16 * t18 + (t11 * t30 - t16 * t48) * pkin(2) + t41, 0, 0, 0, 0, 0, 0, t3, t4, 0, t11 * (t30 * pkin(2) + pkin(3)) + (-pkin(2) * t48 - t18) * t14 + t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, t11 * pkin(3) + (t14 - t16) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tauc_reg  = t7;
