% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% T_reg [1x14]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:59
% EndTime: 2019-12-31 16:21:59
% DurationCPUTime: 0.02s
% Computational Cost: add. (14->11), mult. (46->28), div. (0->0), fcn. (11->2), ass. (0->10)
t38 = qJD(2) ^ 2;
t41 = t38 / 0.2e1;
t40 = t38 * qJ(3);
t39 = qJD(2) * qJD(4);
t37 = cos(qJ(4));
t36 = sin(qJ(4));
t35 = qJD(1) ^ 2 / 0.2e1;
t34 = -qJD(2) * pkin(2) + qJD(3);
t33 = qJD(3) + (-pkin(2) - pkin(5)) * qJD(2);
t1 = [t35, t41, 0, 0, t34 * qJD(2), t40, t35 + qJ(3) ^ 2 * t41 + t34 ^ 2 / 0.2e1, t37 ^ 2 * t41, -t37 * t38 * t36, t37 * t39, -t36 * t39, qJD(4) ^ 2 / 0.2e1, t36 * t40 + (-t36 * qJD(1) + t37 * t33) * qJD(4), t37 * t40 - (t37 * qJD(1) + t36 * t33) * qJD(4);];
T_reg = t1;
