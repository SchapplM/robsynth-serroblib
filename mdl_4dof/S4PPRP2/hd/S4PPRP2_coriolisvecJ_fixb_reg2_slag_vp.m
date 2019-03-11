% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta2]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PPRP2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:13:21
% EndTime: 2019-03-08 18:13:21
% DurationCPUTime: 0.09s
% Computational Cost: add. (71->26), mult. (213->33), div. (0->0), fcn. (164->4), ass. (0->24)
t13 = cos(pkin(5));
t15 = cos(qJ(3));
t21 = t15 * t13;
t12 = sin(pkin(5));
t14 = sin(qJ(3));
t22 = t14 * t12;
t9 = -t21 + t22;
t5 = t9 * qJD(1);
t24 = qJD(4) + t5;
t10 = t15 * t12 + t14 * t13;
t6 = t10 * qJD(1);
t8 = t10 * qJD(3);
t4 = qJD(1) * t8;
t23 = t4 * t9;
t7 = t9 * qJD(3);
t20 = t7 * qJD(3);
t19 = t8 * qJD(3);
t18 = qJD(1) * t22;
t17 = -t5 + t18;
t11 = qJD(3) * qJD(1) * t21;
t3 = qJD(3) * qJ(4) + t6;
t2 = -qJD(3) * pkin(3) + t24;
t1 = t11 + (qJD(4) - t18) * qJD(3);
t16 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, t20, 0 (-qJD(3) * t18 + t11) * t10 - t6 * t7 + t23 + t5 * t8, 0, 0, 0, 0, 0, 0, -t19, 0, -t20, t1 * t10 + t2 * t8 - t3 * t7 + t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 * qJD(3) - t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 + (0.2e1 * qJD(4) - t17) * qJD(3), -t4 * pkin(3) + t1 * qJ(4) - t2 * t6 + t24 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) ^ 2 (-t3 + t6) * qJD(3);];
tauc_reg  = t16;
