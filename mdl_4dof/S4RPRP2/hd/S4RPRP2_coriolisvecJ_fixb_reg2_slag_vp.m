% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRP2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:30:53
% EndTime: 2019-03-08 18:30:54
% DurationCPUTime: 0.16s
% Computational Cost: add. (230->40), mult. (390->53), div. (0->0), fcn. (151->2), ass. (0->33)
t25 = sin(qJ(3));
t26 = cos(qJ(3));
t27 = -pkin(1) - pkin(2);
t32 = -t25 * qJ(2) + t26 * t27;
t21 = t27 * qJD(1) + qJD(2);
t40 = t25 * qJD(2);
t41 = qJD(3) * t26;
t42 = qJD(3) * t25;
t10 = -(qJ(2) * t41 + t40) * qJD(1) - t21 * t42;
t13 = t26 * qJD(2) + t32 * qJD(3);
t39 = qJD(1) * qJ(2);
t16 = t25 * t21 + t26 * t39;
t30 = t26 * qJ(2) + t25 * t27;
t33 = t25 * t39;
t38 = qJD(1) * qJD(2);
t9 = -qJD(3) * t33 + t21 * t41 + t26 * t38;
t45 = t16 * t13 + t9 * t30;
t44 = t16 * t26;
t37 = qJD(1) - qJD(3);
t36 = t10 * t26 + t16 * t41 + t9 * t25;
t34 = 0.2e1 * t38;
t31 = t37 ^ 2;
t15 = t26 * t21 - t33;
t28 = qJD(1) ^ 2;
t18 = t26 * t31;
t17 = t25 * t31;
t14 = -t30 * qJD(3) - t40;
t11 = -pkin(3) * t37 + t15;
t4 = -t16 * t37 + t10;
t3 = -t15 * t37 - t9;
t2 = -t14 * t37 - t10;
t1 = t13 * t37 + t9;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, qJ(2) * t34, 0, 0, 0, 0, 0, 0, t2, t1, 0, t10 * t32 + t15 * t14 + t45, 0, 0, 0, 0, 0, 0, t2, t1, 0, t10 * (-pkin(3) + t32) + t11 * t14 + t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t28 * qJ(2), 0, 0, 0, 0, 0, 0, -t17, -t18, 0, -t15 * t42 + (t15 * t25 - t44) * qJD(1) + t36, 0, 0, 0, 0, 0, 0, -t17, -t18, 0, -t11 * t42 + (t11 * t25 - t44) * qJD(1) + t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, t10 * pkin(3) + (t11 - t15) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tauc_reg  = t5;
