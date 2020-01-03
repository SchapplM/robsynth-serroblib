% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR11_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:35:13
% EndTime: 2019-12-31 21:35:13
% DurationCPUTime: 0.08s
% Computational Cost: add. (178->39), mult. (395->83), div. (0->0), fcn. (241->6), ass. (0->35)
t159 = pkin(3) + pkin(4);
t147 = qJD(1) ^ 2;
t158 = t147 / 0.2e1;
t146 = cos(qJ(2));
t157 = t146 * t147;
t143 = sin(qJ(2));
t128 = (-pkin(2) * t146 - pkin(7) * t143 - pkin(1)) * qJD(1);
t154 = t146 * qJD(1);
t134 = pkin(6) * t154 + qJD(2) * pkin(7);
t142 = sin(qJ(3));
t145 = cos(qJ(3));
t156 = t142 * t128 + t145 * t134;
t155 = qJD(1) * t143;
t153 = qJD(1) * qJD(2);
t137 = -qJD(3) + t154;
t122 = -t137 * qJ(4) + t156;
t152 = t143 * t153;
t151 = t146 * t153;
t133 = -qJD(2) * pkin(2) + pkin(6) * t155;
t150 = t145 * t128 - t142 * t134;
t149 = qJD(4) - t150;
t130 = t142 * qJD(2) + t145 * t155;
t148 = t130 * qJ(4) - t133;
t144 = cos(qJ(5));
t141 = sin(qJ(5));
t136 = qJD(5) + t137;
t129 = -t145 * qJD(2) + t142 * t155;
t125 = t141 * t129 + t144 * t130;
t124 = -t144 * t129 + t141 * t130;
t123 = t129 * pkin(3) - t148;
t121 = t137 * pkin(3) + t149;
t120 = -t159 * t129 + t148;
t119 = t129 * pkin(8) + t122;
t118 = -t130 * pkin(8) + t159 * t137 + t149;
t1 = [t158, 0, 0, t143 ^ 2 * t158, t143 * t157, t152, t151, qJD(2) ^ 2 / 0.2e1, pkin(1) * t157 - pkin(6) * t152, -t147 * pkin(1) * t143 - pkin(6) * t151, t130 ^ 2 / 0.2e1, -t130 * t129, -t130 * t137, t129 * t137, t137 ^ 2 / 0.2e1, t133 * t129 - t150 * t137, t133 * t130 + t156 * t137, t121 * t137 + t123 * t129, t121 * t130 - t122 * t129, -t122 * t137 - t123 * t130, t122 ^ 2 / 0.2e1 + t123 ^ 2 / 0.2e1 + t121 ^ 2 / 0.2e1, t125 ^ 2 / 0.2e1, -t125 * t124, t125 * t136, -t124 * t136, t136 ^ 2 / 0.2e1, (t144 * t118 - t141 * t119) * t136 + t120 * t124, -(t141 * t118 + t144 * t119) * t136 + t120 * t125;];
T_reg = t1;
