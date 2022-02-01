% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:31
% EndTime: 2022-01-20 10:48:31
% DurationCPUTime: 0.07s
% Computational Cost: add. (120->28), mult. (198->67), div. (0->0), fcn. (107->8), ass. (0->29)
t34 = qJD(1) + qJD(2);
t32 = t34 ^ 2;
t50 = t32 / 0.2e1;
t38 = sin(qJ(4));
t49 = t34 * t38;
t41 = cos(qJ(4));
t48 = t34 * t41;
t46 = pkin(1) * qJD(1);
t43 = cos(qJ(2)) * t46;
t26 = t34 * pkin(2) + t43;
t35 = sin(pkin(9));
t36 = cos(pkin(9));
t44 = sin(qJ(2)) * t46;
t22 = t35 * t26 + t36 * t44;
t20 = t34 * pkin(7) + t22;
t47 = t38 * qJD(3) + t41 * t20;
t45 = qJD(4) * t34;
t21 = t36 * t26 - t35 * t44;
t40 = cos(qJ(5));
t37 = sin(qJ(5));
t33 = qJD(4) + qJD(5);
t31 = t41 * qJD(3);
t24 = (t37 * t41 + t38 * t40) * t34;
t23 = t37 * t49 - t40 * t48;
t19 = -t34 * pkin(3) - t21;
t17 = (-pkin(4) * t41 - pkin(3)) * t34 - t21;
t16 = pkin(8) * t48 + t47;
t15 = qJD(4) * pkin(4) + t31 + (-pkin(8) * t34 - t20) * t38;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t50, t34 * t43, -t34 * t44, t22 ^ 2 / 0.2e1 + t21 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t38 ^ 2 * t50, t38 * t32 * t41, t38 * t45, t41 * t45, qJD(4) ^ 2 / 0.2e1, -t19 * t48 + (-t38 * t20 + t31) * qJD(4), -t47 * qJD(4) + t19 * t49, t24 ^ 2 / 0.2e1, -t24 * t23, t24 * t33, -t23 * t33, t33 ^ 2 / 0.2e1, t17 * t23 + (t40 * t15 - t37 * t16) * t33, t17 * t24 - (t37 * t15 + t40 * t16) * t33;];
T_reg = t1;
