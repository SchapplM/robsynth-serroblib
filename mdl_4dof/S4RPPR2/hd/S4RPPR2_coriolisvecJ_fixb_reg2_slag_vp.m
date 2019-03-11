% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPPR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:28:32
% EndTime: 2019-03-08 18:28:33
% DurationCPUTime: 0.19s
% Computational Cost: add. (277->42), mult. (496->68), div. (0->0), fcn. (266->4), ass. (0->37)
t47 = qJD(1) - qJD(4);
t26 = -pkin(1) - pkin(2);
t19 = t26 * qJD(1) + qJD(2);
t22 = sin(pkin(6));
t23 = cos(pkin(6));
t40 = qJD(1) * qJ(2);
t10 = t22 * t19 + t23 * t40;
t24 = sin(qJ(4));
t25 = cos(qJ(4));
t17 = t23 * t19;
t41 = t22 * qJ(2);
t8 = t17 + (-pkin(3) - t41) * qJD(1);
t6 = t25 * t10 + t24 * t8;
t46 = t6 * qJD(4);
t39 = qJD(1) * qJD(2);
t36 = t23 * t39;
t37 = t22 * t39;
t45 = t25 * t8;
t1 = -(qJD(4) * t10 + t37) * t24 + qJD(4) * t45 + t25 * t36;
t30 = -t24 * t22 + t25 * t23;
t43 = t47 * t30;
t31 = t25 * t22 + t24 * t23;
t42 = t47 * t31;
t38 = 0.2e1 * t39;
t35 = t23 * t26 - t41;
t34 = -t10 * t23 + (-t22 * t40 + t17) * t22;
t15 = -pkin(3) + t35;
t16 = t23 * qJ(2) + t22 * t26;
t33 = t25 * t15 - t24 * t16;
t32 = t24 * t15 + t25 * t16;
t29 = t31 * qJD(2);
t2 = -qJD(1) * t29 - t46;
t27 = qJD(1) ^ 2;
t5 = -t24 * t10 + t45;
t4 = -t32 * qJD(4) - t29;
t3 = t30 * qJD(2) + t33 * qJD(4);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, qJ(2) * t38, 0, 0, 0, 0, 0, 0, 0.2e1 * t37, 0.2e1 * t36, 0 ((t23 * t16 - t22 * t35) * qJD(1) - t34) * qJD(2), 0, 0, 0, 0, 0, 0, t31 * t39 - t4 * t47 + t46, t3 * t47 + t1, 0, t1 * t32 + t2 * t33 + t6 * t3 + t5 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t27 * qJ(2), 0, 0, 0, 0, 0, 0, -t22 * t27, -t23 * t27, 0, t34 * qJD(1), 0, 0, 0, 0, 0, 0, -t42 * t47, -t43 * t47, 0, t1 * t31 + t2 * t30 + t42 * t5 - t43 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47 * t6 + t2, -t47 * t5 - t1, 0, 0;];
tauc_reg  = t7;
