% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:14:16
% EndTime: 2019-12-31 21:14:16
% DurationCPUTime: 0.13s
% Computational Cost: add. (250->39), mult. (623->86), div. (0->0), fcn. (440->8), ass. (0->39)
t163 = -pkin(7) - pkin(6);
t152 = qJD(1) ^ 2;
t162 = t152 / 0.2e1;
t161 = cos(qJ(3));
t151 = cos(qJ(2));
t160 = t151 * t152;
t148 = sin(qJ(3));
t149 = sin(qJ(2));
t137 = (t148 * t151 + t161 * t149) * qJD(1);
t144 = qJD(2) + qJD(3);
t158 = qJD(1) * t149;
t139 = qJD(2) * pkin(2) + t163 * t158;
t157 = qJD(1) * t151;
t140 = t163 * t157;
t153 = t161 * t139 + t148 * t140;
t124 = t144 * pkin(3) - t137 * qJ(4) + t153;
t136 = t148 * t158 - t161 * t157;
t159 = t148 * t139 - t161 * t140;
t126 = -t136 * qJ(4) + t159;
t145 = sin(pkin(9));
t146 = cos(pkin(9));
t121 = t145 * t124 + t146 * t126;
t156 = qJD(1) * qJD(2);
t155 = t149 * t156;
t154 = t151 * t156;
t130 = -t146 * t136 - t145 * t137;
t141 = (-pkin(2) * t151 - pkin(1)) * qJD(1);
t120 = t146 * t124 - t145 * t126;
t132 = t136 * pkin(3) + qJD(4) + t141;
t150 = cos(qJ(5));
t147 = sin(qJ(5));
t131 = -t145 * t136 + t146 * t137;
t129 = qJD(5) - t130;
t128 = t150 * t131 + t147 * t144;
t127 = t147 * t131 - t150 * t144;
t122 = -t130 * pkin(4) - t131 * pkin(8) + t132;
t119 = t144 * pkin(8) + t121;
t118 = -t144 * pkin(4) - t120;
t1 = [t162, 0, 0, t149 ^ 2 * t162, t149 * t160, t155, t154, qJD(2) ^ 2 / 0.2e1, pkin(1) * t160 - pkin(6) * t155, -t152 * pkin(1) * t149 - pkin(6) * t154, t137 ^ 2 / 0.2e1, -t137 * t136, t137 * t144, -t136 * t144, t144 ^ 2 / 0.2e1, t141 * t136 + t153 * t144, t141 * t137 - t159 * t144, -t120 * t131 + t121 * t130, t121 ^ 2 / 0.2e1 + t120 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1, t128 ^ 2 / 0.2e1, -t128 * t127, t128 * t129, -t127 * t129, t129 ^ 2 / 0.2e1, (-t147 * t119 + t150 * t122) * t129 + t118 * t127, -(t150 * t119 + t147 * t122) * t129 + t118 * t128;];
T_reg = t1;
