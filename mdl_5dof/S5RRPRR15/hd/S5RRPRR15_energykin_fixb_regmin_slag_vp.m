% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR15_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:43:07
% EndTime: 2019-12-31 20:43:08
% DurationCPUTime: 0.08s
% Computational Cost: add. (143->39), mult. (327->85), div. (0->0), fcn. (186->6), ass. (0->35)
t162 = -pkin(2) - pkin(7);
t150 = qJD(1) ^ 2;
t161 = t150 / 0.2e1;
t149 = cos(qJ(2));
t160 = t149 * t150;
t146 = sin(qJ(2));
t152 = -qJ(3) * t146 - pkin(1);
t129 = (t162 * t149 + t152) * qJD(1);
t157 = t146 * qJD(1);
t156 = pkin(6) * t157 + qJD(3);
t130 = pkin(3) * t157 + t162 * qJD(2) + t156;
t145 = sin(qJ(4));
t148 = cos(qJ(4));
t159 = t148 * t129 + t145 * t130;
t158 = qJD(1) * t149;
t137 = -pkin(6) * t158 - qJD(2) * qJ(3);
t155 = qJD(1) * qJD(2);
t132 = pkin(3) * t158 - t137;
t154 = t146 * t155;
t153 = t149 * t155;
t151 = -t145 * t129 + t148 * t130;
t139 = qJD(4) + t157;
t147 = cos(qJ(5));
t144 = sin(qJ(5));
t138 = qJD(5) + t139;
t136 = -qJD(2) * pkin(2) + t156;
t135 = t148 * qJD(2) - t145 * t158;
t134 = t145 * qJD(2) + t148 * t158;
t133 = (-pkin(2) * t149 + t152) * qJD(1);
t125 = t134 * pkin(4) + t132;
t124 = -t144 * t134 + t147 * t135;
t123 = t147 * t134 + t144 * t135;
t122 = -t134 * pkin(8) + t159;
t121 = t139 * pkin(4) - t135 * pkin(8) + t151;
t1 = [t161, 0, 0, t146 ^ 2 * t161, t146 * t160, t154, t153, qJD(2) ^ 2 / 0.2e1, pkin(1) * t160 - pkin(6) * t154, -t150 * pkin(1) * t146 - pkin(6) * t153, (t136 * t146 - t137 * t149) * qJD(1), t136 * qJD(2) + t133 * t158, -t137 * qJD(2) - t133 * t157, t133 ^ 2 / 0.2e1 + t137 ^ 2 / 0.2e1 + t136 ^ 2 / 0.2e1, t135 ^ 2 / 0.2e1, -t135 * t134, t135 * t139, -t134 * t139, t139 ^ 2 / 0.2e1, t132 * t134 + t151 * t139, t132 * t135 - t159 * t139, t124 ^ 2 / 0.2e1, -t124 * t123, t124 * t138, -t123 * t138, t138 ^ 2 / 0.2e1, (t147 * t121 - t144 * t122) * t138 + t125 * t123, -(t144 * t121 + t147 * t122) * t138 + t125 * t124;];
T_reg = t1;
