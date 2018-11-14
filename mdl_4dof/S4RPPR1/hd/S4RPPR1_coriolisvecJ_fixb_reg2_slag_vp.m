% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:47
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPPR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:46:46
% EndTime: 2018-11-14 13:46:46
% DurationCPUTime: 0.11s
% Computational Cost: add. (165->26), mult. (290->39), div. (0->0), fcn. (129->4), ass. (0->20)
t16 = sin(qJ(4));
t17 = cos(qJ(4));
t10 = -cos(pkin(6)) * pkin(1) - pkin(2) - pkin(3);
t8 = t10 * qJD(1) + qJD(3);
t11 = sin(pkin(6)) * pkin(1) + qJ(3);
t9 = qJD(1) * t11;
t3 = -t16 * t9 + t17 * t8;
t24 = qJD(3) * qJD(1);
t4 = t16 * t8 + t17 * t9;
t2 = -qJD(4) * t4 - t16 * t24;
t23 = qJD(1) - qJD(4);
t28 = -t4 * t23 + t2;
t18 = -qJD(4) * t3 - t17 * t24;
t27 = t23 * t3 - t18;
t21 = t23 ^ 2;
t20 = t17 * t10 - t16 * t11;
t19 = t16 * t10 + t17 * t11;
t6 = -t16 * qJD(3) - t19 * qJD(4);
t5 = t17 * qJD(3) + t20 * qJD(4);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t24, 0.2e1 * t9 * qJD(3), 0, 0, 0, 0, 0, 0, -t23 * t6 - t2, t23 * t5 - t18, 0, -t18 * t19 + t2 * t20 + t3 * t6 + t4 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1) ^ 2, -t9 * qJD(1), 0, 0, 0, 0, 0, 0, -t16 * t21, -t17 * t21, 0, t27 * t16 + t28 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, 0, 0;];
tauc_reg  = t1;
