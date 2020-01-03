% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
% 
% Output:
% tauc_reg [5x15]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRPR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:30
% EndTime: 2019-12-31 17:33:30
% DurationCPUTime: 0.12s
% Computational Cost: add. (74->31), mult. (190->50), div. (0->0), fcn. (91->4), ass. (0->29)
t8 = sin(qJ(3));
t25 = t8 * qJD(2);
t23 = qJD(3) * qJ(4);
t4 = t23 + t25;
t19 = -t4 + t25;
t33 = qJD(5) * (-t19 + t23);
t16 = t19 * qJD(3);
t7 = sin(qJ(5));
t9 = cos(qJ(5));
t32 = t7 ^ 2 - t9 ^ 2;
t12 = qJD(5) ^ 2;
t31 = t12 * t7;
t30 = t12 * t9;
t13 = qJD(3) ^ 2;
t29 = t13 * t8;
t10 = cos(qJ(3));
t28 = t13 * t10;
t27 = t12 + t13;
t26 = qJD(3) * pkin(3);
t24 = t10 * qJD(2);
t22 = qJD(3) * qJD(5);
t21 = 0.2e1 * t22;
t20 = t10 * t27;
t18 = -0.2e1 * t7 * t22;
t17 = qJD(4) - t24;
t2 = (qJD(4) + t24) * qJD(3);
t14 = t17 * qJD(3) - (-pkin(3) - pkin(6)) * t12 + t2;
t3 = t17 - t26;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t31; 0, 0, 0, -t29, -t28, t29, t28, t2 * t8 + (-t19 * t10 + t3 * t8) * qJD(3), 0, 0, 0, 0, 0, t8 * t9 * t21 + t7 * t20, t8 * t18 + t9 * t20; 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * qJD(4), t2 * qJ(4) + t4 * qJD(4) + (-t10 * t4 + (-t3 - t26) * t8) * qJD(2), t9 * t18, t32 * t21, -t31, -t30, 0, t14 * t7 + t9 * t33, t14 * t9 - t7 * t33; 0, 0, 0, 0, 0, 0, -t13, t16, 0, 0, 0, 0, 0, -t27 * t7, -t27 * t9; 0, 0, 0, 0, 0, 0, 0, 0, t9 * t13 * t7, -t32 * t13, 0, 0, 0, t9 * t16, -t7 * t16;];
tauc_reg = t1;
