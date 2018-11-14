% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% taug_reg [4x8]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:58
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4PPRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:57:27
% EndTime: 2018-11-14 13:57:28
% DurationCPUTime: 0.10s
% Computational Cost: add. (55->24), mult. (159->33), div. (0->0), fcn. (116->4), ass. (0->26)
t13 = cos(pkin(5));
t15 = cos(qJ(3));
t12 = sin(pkin(5));
t14 = sin(qJ(3));
t26 = t14 * t12;
t16 = -t13 * t15 + t26;
t5 = t16 * qJD(1);
t27 = qJD(4) + t5;
t22 = qJD(1) * qJD(3);
t20 = t15 * t22;
t25 = t14 * t13;
t4 = t12 * t20 + t22 * t25;
t7 = t16 * qJD(3);
t24 = t7 * qJD(3);
t17 = t12 * t15 + t25;
t8 = t17 * qJD(3);
t23 = t8 * qJD(3);
t21 = qJD(1) * t26;
t6 = t17 * qJD(1);
t19 = qJD(3) * t6 - t4;
t18 = -t5 + t21;
t11 = t13 * t20;
t3 = qJD(3) * qJ(4) + t6;
t2 = -qJD(3) * pkin(3) + t27;
t1 = t11 + (qJD(4) - t21) * qJD(3);
t9 = [0, 0, 0, -t23, t24, -t23, -t24, t1 * t17 + t4 * t16 + t2 * t8 - t3 * t7; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t19, t18 * qJD(3) - t11, t19, t11 + (0.2e1 * qJD(4) - t18) * qJD(3), -pkin(3) * t4 + qJ(4) * t1 - t2 * t6 + t27 * t3; 0, 0, 0, 0, 0, 0, -qJD(3) ^ 2, -t3 * qJD(3) + t4;];
tauc_reg  = t9;
