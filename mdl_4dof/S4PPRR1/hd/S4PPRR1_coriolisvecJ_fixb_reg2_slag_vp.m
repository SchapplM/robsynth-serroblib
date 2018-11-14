% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:40
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4PPRR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:40:19
% EndTime: 2018-11-14 13:40:19
% DurationCPUTime: 0.12s
% Computational Cost: add. (100->30), mult. (250->48), div. (0->0), fcn. (176->4), ass. (0->28)
t12 = qJD(3) + qJD(4);
t15 = cos(qJ(4));
t16 = cos(qJ(3));
t26 = t16 * qJD(2);
t24 = qJD(3) * t26;
t11 = qJD(3) * pkin(3) + t26;
t27 = qJD(4) * t11;
t13 = sin(qJ(4));
t14 = sin(qJ(3));
t28 = qJD(2) * t14;
t25 = t13 * t28;
t30 = t12 * t25;
t1 = (t24 + t27) * t15 - t30;
t29 = t15 * t14;
t23 = t13 * t16 + t29;
t22 = -t13 * t14 + t15 * t16;
t21 = (-pkin(3) * t12 - t11) * qJD(4);
t20 = t23 * qJD(3);
t18 = (-qJD(4) * t29 - t20) * qJD(2);
t2 = -t13 * t27 + t18;
t17 = qJD(3) ^ 2;
t8 = t22 * qJD(2);
t7 = t23 * qJD(2);
t6 = t13 * t11 + t15 * t28;
t5 = t15 * t11 - t25;
t4 = -t23 * qJD(4) - t20;
t3 = t12 * t22;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17 * t14, -t17 * t16, 0, 0, 0, 0, 0, 0, 0, 0, t4 * t12, -t3 * t12, 0, t1 * t23 + t2 * t22 + t6 * t3 + t5 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * t12 + t13 * t21 + t18, t8 * t12 + (t21 - t24) * t15 + t30, 0, t5 * t7 - t6 * t8 + (t1 * t13 + t2 * t15 + (-t13 * t5 + t15 * t6) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * t12 + t2, t5 * t12 - t1, 0, 0;];
tauc_reg  = t9;
