% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPPR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:39:14
% EndTime: 2019-12-31 19:39:14
% DurationCPUTime: 0.12s
% Computational Cost: add. (165->40), mult. (386->87), div. (0->0), fcn. (214->6), ass. (0->32)
t150 = qJD(1) ^ 2;
t158 = t150 / 0.2e1;
t149 = cos(qJ(2));
t157 = t149 * t150;
t147 = sin(qJ(2));
t156 = qJD(1) * t147;
t154 = pkin(6) * t156 + qJD(3);
t126 = -qJ(4) * t156 + (-pkin(2) - pkin(3)) * qJD(2) + t154;
t155 = qJD(1) * t149;
t134 = pkin(6) * t155 + qJD(2) * qJ(3);
t131 = -qJ(4) * t155 + t134;
t143 = sin(pkin(8));
t144 = cos(pkin(8));
t119 = t143 * t126 + t144 * t131;
t153 = qJD(1) * qJD(2);
t132 = -qJD(1) * pkin(1) - pkin(2) * t155 - qJ(3) * t156;
t152 = t147 * t153;
t151 = t149 * t153;
t118 = t144 * t126 - t143 * t131;
t125 = pkin(3) * t155 + qJD(4) - t132;
t148 = cos(qJ(5));
t146 = sin(qJ(5));
t140 = qJD(2) - qJD(5);
t133 = -qJD(2) * pkin(2) + t154;
t130 = (-t143 * t149 + t144 * t147) * qJD(1);
t129 = (t143 * t147 + t144 * t149) * qJD(1);
t122 = -t146 * t129 + t148 * t130;
t121 = t148 * t129 + t146 * t130;
t120 = t129 * pkin(4) + t125;
t117 = -t129 * pkin(7) + t119;
t116 = -qJD(2) * pkin(4) - t130 * pkin(7) + t118;
t1 = [t158, 0, 0, t147 ^ 2 * t158, t147 * t157, t152, t151, qJD(2) ^ 2 / 0.2e1, pkin(1) * t157 - pkin(6) * t152, -t150 * pkin(1) * t147 - pkin(6) * t151, -t133 * qJD(2) - t132 * t155, (t133 * t147 + t134 * t149) * qJD(1), t134 * qJD(2) - t132 * t156, t134 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1 + t133 ^ 2 / 0.2e1, -t118 * qJD(2) + t125 * t129, t119 * qJD(2) + t125 * t130, -t118 * t130 - t119 * t129, t119 ^ 2 / 0.2e1 + t118 ^ 2 / 0.2e1 + t125 ^ 2 / 0.2e1, t122 ^ 2 / 0.2e1, -t122 * t121, -t122 * t140, t121 * t140, t140 ^ 2 / 0.2e1, t120 * t121 - (t148 * t116 - t146 * t117) * t140, t120 * t122 + (t146 * t116 + t148 * t117) * t140;];
T_reg = t1;
