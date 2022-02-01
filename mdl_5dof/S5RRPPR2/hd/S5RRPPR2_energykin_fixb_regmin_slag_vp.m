% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:50
% EndTime: 2022-01-20 10:05:50
% DurationCPUTime: 0.07s
% Computational Cost: add. (121->26), mult. (198->61), div. (0->0), fcn. (103->8), ass. (0->27)
t26 = qJD(1) + qJD(2);
t24 = t26 ^ 2;
t27 = sin(pkin(9));
t43 = t24 * t27 ^ 2;
t42 = t26 * t27;
t29 = cos(pkin(9));
t41 = t29 * t26;
t40 = pkin(1) * qJD(1);
t36 = cos(qJ(2)) * t40;
t18 = t26 * pkin(2) + t36;
t28 = sin(pkin(8));
t30 = cos(pkin(8));
t37 = sin(qJ(2)) * t40;
t16 = t28 * t18 + t30 * t37;
t31 = sin(qJ(5));
t39 = t31 * t42;
t33 = cos(qJ(5));
t38 = t33 * t42;
t15 = t30 * t18 - t28 * t37;
t35 = qJD(4) - t15;
t19 = -qJD(5) + t41;
t14 = t26 * qJ(4) + t16;
t13 = -t26 * pkin(3) + t35;
t12 = t27 * qJD(3) + t29 * t14;
t10 = -t29 * qJD(3) + t27 * t14;
t9 = (-pkin(4) * t29 - pkin(7) * t27 - pkin(3)) * t26 + t35;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t24 / 0.2e1, t26 * t36, -t26 * t37, t16 ^ 2 / 0.2e1 + t15 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, -t13 * t41, (t10 * t27 + t12 * t29) * t26, t12 ^ 2 / 0.2e1 + t10 ^ 2 / 0.2e1 + t13 ^ 2 / 0.2e1, t33 ^ 2 * t43 / 0.2e1, -t33 * t31 * t43, -t19 * t38, t19 * t39, t19 ^ 2 / 0.2e1, -(-t31 * t12 + t33 * t9) * t19 + t10 * t39, (t33 * t12 + t31 * t9) * t19 + t10 * t38;];
T_reg = t1;
