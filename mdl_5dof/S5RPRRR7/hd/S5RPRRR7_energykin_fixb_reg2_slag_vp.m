% Calculate inertial parameters regressor of fixed base kinetic energy for
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x(5*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR7_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_energykin_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_energykin_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_energykin_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:04:08
% EndTime: 2019-12-31 19:04:08
% DurationCPUTime: 0.16s
% Computational Cost: add. (280->44), mult. (648->115), div. (0->0), fcn. (380->8), ass. (0->36)
t40 = qJD(1) ^ 2;
t33 = t40 / 0.2e1;
t48 = pkin(1) * t40;
t47 = cos(qJ(4));
t46 = cos(qJ(5));
t34 = sin(pkin(9));
t25 = (pkin(1) * t34 + pkin(6)) * qJD(1);
t38 = sin(qJ(3));
t39 = cos(qJ(3));
t18 = t38 * qJD(2) + t39 * t25;
t15 = qJD(3) * pkin(7) + t18;
t35 = cos(pkin(9));
t42 = -pkin(1) * t35 - pkin(2);
t16 = (-pkin(3) * t39 - pkin(7) * t38 + t42) * qJD(1);
t37 = sin(qJ(4));
t6 = t47 * t15 + t37 * t16;
t45 = qJD(1) * t38;
t44 = t39 * qJD(1);
t43 = qJD(1) * qJD(3);
t5 = -t37 * t15 + t47 * t16;
t17 = t39 * qJD(2) - t38 * t25;
t29 = -qJD(4) + t44;
t14 = -qJD(3) * pkin(3) - t17;
t36 = sin(qJ(5));
t27 = -qJD(5) + t29;
t26 = t42 * qJD(1);
t24 = t37 * qJD(3) + t47 * t45;
t22 = -t47 * qJD(3) + t37 * t45;
t10 = -t36 * t22 + t46 * t24;
t8 = t46 * t22 + t36 * t24;
t7 = t22 * pkin(4) + t14;
t4 = -t22 * pkin(8) + t6;
t3 = -t29 * pkin(4) - t24 * pkin(8) + t5;
t2 = t36 * t3 + t46 * t4;
t1 = t46 * t3 - t36 * t4;
t9 = [0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t35 * t48, -t34 * t48, 0, qJD(2) ^ 2 / 0.2e1 + (t34 ^ 2 / 0.2e1 + t35 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t40, t38 ^ 2 * t33, t38 * t40 * t39, t38 * t43, t39 ^ 2 * t33, t39 * t43, qJD(3) ^ 2 / 0.2e1, t17 * qJD(3) - t26 * t44, -t18 * qJD(3) + t26 * t45, (-t17 * t38 + t18 * t39) * qJD(1), t18 ^ 2 / 0.2e1 + t17 ^ 2 / 0.2e1 + t26 ^ 2 / 0.2e1, t24 ^ 2 / 0.2e1, -t24 * t22, -t24 * t29, t22 ^ 2 / 0.2e1, t22 * t29, t29 ^ 2 / 0.2e1, t14 * t22 - t5 * t29, t14 * t24 + t6 * t29, -t6 * t22 - t5 * t24, t6 ^ 2 / 0.2e1 + t5 ^ 2 / 0.2e1 + t14 ^ 2 / 0.2e1, t10 ^ 2 / 0.2e1, -t10 * t8, -t10 * t27, t8 ^ 2 / 0.2e1, t8 * t27, t27 ^ 2 / 0.2e1, -t1 * t27 + t7 * t8, t7 * t10 + t2 * t27, -t1 * t10 - t2 * t8, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t7 ^ 2 / 0.2e1;];
T_reg = t9;
