% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRRR12_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:08:19
% EndTime: 2024-09-28 18:08:19
% DurationCPUTime: 0.00s
% Computational Cost: add. (125->19), mult. (226->55), div. (0->0), fcn. (166->12), ass. (0->31)
t19 = qJD(2) + qJD(3);
t17 = qJD(4) + t19;
t22 = cos(pkin(6));
t15 = t22 * t17 + qJD(5);
t20 = sin(pkin(6));
t40 = t17 * t20;
t31 = cos(qJ(2));
t39 = qJD(1) * sin(pkin(5));
t14 = qJD(2) * pkin(2) + t31 * t39;
t26 = sin(qJ(3));
t30 = cos(qJ(3));
t27 = sin(qJ(2));
t36 = t27 * t39;
t33 = t30 * t14 - t26 * t36;
t10 = t19 * pkin(3) + t33;
t12 = t26 * t14 + t30 * t36;
t25 = sin(qJ(4));
t29 = cos(qJ(4));
t42 = t25 * t10 + t29 * t12;
t43 = t15 * (pkin(10) * t40 + t42);
t16 = t17 ^ 2;
t41 = t16 * t20 ^ 2;
t38 = qJD(1) * cos(pkin(5));
t37 = t15 * t40;
t35 = t29 * t10 - t25 * t12;
t34 = qJD(2) * t39;
t7 = t17 * pkin(4) + t35;
t32 = (t20 * t38 + t22 * t7) * t15 - (-t20 * t7 + t22 * t38) * t40;
t28 = cos(qJ(5));
t24 = sin(qJ(5));
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t31 * t34, -t27 * t34, t19 ^ 2 / 0.2e1, t33 * t19, -t12 * t19, t16 / 0.2e1, t35 * t17, -t42 * t17, t24 ^ 2 * t41 / 0.2e1, t24 * t28 * t41, t24 * t37, t28 * t37, t15 ^ 2 / 0.2e1, -t24 * t43 + t32 * t28, -t32 * t24 - t28 * t43;];
T_reg = t1;
