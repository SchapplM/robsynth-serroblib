% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RPPR4
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
% T_reg [1x14]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RPPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:53
% EndTime: 2019-12-31 16:38:53
% DurationCPUTime: 0.03s
% Computational Cost: add. (25->16), mult. (73->38), div. (0->0), fcn. (21->4), ass. (0->14)
t52 = qJD(1) ^ 2;
t57 = t52 / 0.2e1;
t48 = sin(pkin(6));
t45 = (-pkin(1) * t48 - qJ(3)) * qJD(1);
t56 = t45 * qJD(1);
t55 = qJD(1) * qJD(4);
t49 = cos(pkin(6));
t54 = -pkin(1) * t49 - pkin(2);
t51 = cos(qJ(4));
t50 = sin(qJ(4));
t47 = qJD(2) ^ 2 / 0.2e1;
t44 = t54 * qJD(1) + qJD(3);
t43 = qJD(3) + (-pkin(5) + t54) * qJD(1);
t1 = [t57, 0, 0, t47 + (t48 ^ 2 / 0.2e1 + t49 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t52, t44 * qJD(1), -t56, t47 + t45 ^ 2 / 0.2e1 + t44 ^ 2 / 0.2e1, t51 ^ 2 * t57, -t51 * t52 * t50, t51 * t55, -t50 * t55, qJD(4) ^ 2 / 0.2e1, -t50 * t56 + (-t50 * qJD(2) + t51 * t43) * qJD(4), -t51 * t56 - (t51 * qJD(2) + t50 * t43) * qJD(4);];
T_reg = t1;
