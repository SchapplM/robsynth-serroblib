% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPR7
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
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:17:15
% EndTime: 2019-12-31 21:17:16
% DurationCPUTime: 0.13s
% Computational Cost: add. (291->41), mult. (653->90), div. (0->0), fcn. (449->8), ass. (0->39)
t166 = -pkin(7) - pkin(6);
t155 = qJD(1) ^ 2;
t165 = t155 / 0.2e1;
t164 = cos(pkin(9));
t154 = cos(qJ(2));
t163 = t154 * t155;
t150 = sin(qJ(3));
t153 = cos(qJ(3));
t160 = qJD(1) * t154;
t151 = sin(qJ(2));
t161 = qJD(1) * t151;
t138 = t150 * t161 - t153 * t160;
t139 = (t150 * t154 + t151 * t153) * qJD(1);
t144 = (-pkin(2) * t154 - pkin(1)) * qJD(1);
t129 = t138 * pkin(3) - t139 * qJ(4) + t144;
t147 = qJD(2) + qJD(3);
t142 = qJD(2) * pkin(2) + t166 * t161;
t143 = t166 * t160;
t162 = t150 * t142 - t153 * t143;
t132 = t147 * qJ(4) + t162;
t148 = sin(pkin(9));
t123 = t148 * t129 + t164 * t132;
t159 = qJD(1) * qJD(2);
t158 = t151 * t159;
t157 = t154 * t159;
t122 = t164 * t129 - t148 * t132;
t156 = t153 * t142 + t150 * t143;
t131 = -t147 * pkin(3) + qJD(4) - t156;
t152 = cos(qJ(5));
t149 = sin(qJ(5));
t137 = qJD(5) + t138;
t135 = t164 * t139 + t148 * t147;
t134 = t148 * t139 - t164 * t147;
t126 = -t149 * t134 + t152 * t135;
t125 = t152 * t134 + t149 * t135;
t124 = t134 * pkin(4) + t131;
t121 = -t134 * pkin(8) + t123;
t120 = t138 * pkin(4) - t135 * pkin(8) + t122;
t1 = [t165, 0, 0, t151 ^ 2 * t165, t151 * t163, t158, t157, qJD(2) ^ 2 / 0.2e1, pkin(1) * t163 - pkin(6) * t158, -t155 * pkin(1) * t151 - pkin(6) * t157, t139 ^ 2 / 0.2e1, -t139 * t138, t139 * t147, -t138 * t147, t147 ^ 2 / 0.2e1, t144 * t138 + t156 * t147, t144 * t139 - t162 * t147, t122 * t138 + t131 * t134, -t123 * t138 + t131 * t135, -t122 * t135 - t123 * t134, t123 ^ 2 / 0.2e1 + t122 ^ 2 / 0.2e1 + t131 ^ 2 / 0.2e1, t126 ^ 2 / 0.2e1, -t126 * t125, t126 * t137, -t125 * t137, t137 ^ 2 / 0.2e1, (t152 * t120 - t149 * t121) * t137 + t124 * t125, -(t149 * t120 + t152 * t121) * t137 + t124 * t126;];
T_reg = t1;
