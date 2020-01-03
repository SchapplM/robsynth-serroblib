% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:43:01
% EndTime: 2020-01-03 11:43:01
% DurationCPUTime: 0.11s
% Computational Cost: add. (190->41), mult. (508->93), div. (0->0), fcn. (336->8), ass. (0->36)
t151 = sin(pkin(8));
t148 = t151 ^ 2;
t158 = qJD(1) ^ 2;
t166 = t148 * t158;
t153 = cos(pkin(8));
t138 = qJD(2) + (-pkin(2) * t153 - pkin(6) * t151 - pkin(1)) * qJD(1);
t157 = cos(qJ(3));
t137 = t157 * t138;
t163 = t153 * qJD(1);
t144 = -qJD(3) + t163;
t155 = sin(qJ(3));
t129 = -t144 * pkin(3) + t137 + (-qJ(2) * t153 * t155 - qJ(4) * t151 * t157) * qJD(1);
t164 = qJD(1) * t151;
t160 = t155 * t164;
t161 = qJ(2) * t163;
t165 = t155 * t138 + t157 * t161;
t132 = -qJ(4) * t160 + t165;
t150 = sin(pkin(9));
t152 = cos(pkin(9));
t124 = t150 * t129 + t152 * t132;
t162 = t155 * t166;
t139 = pkin(3) * t160 + qJ(2) * t164 + qJD(4);
t123 = t152 * t129 - t150 * t132;
t156 = cos(qJ(5));
t154 = sin(qJ(5));
t149 = t153 ^ 2;
t147 = -qJD(1) * pkin(1) + qJD(2);
t141 = -qJD(5) + t144;
t135 = (-t150 * t155 + t152 * t157) * t164;
t134 = (-t150 * t157 - t152 * t155) * t164;
t131 = -t134 * pkin(4) + t139;
t126 = t154 * t134 + t156 * t135;
t125 = -t156 * t134 + t154 * t135;
t122 = t134 * pkin(7) + t124;
t121 = -t144 * pkin(4) - t135 * pkin(7) + t123;
t1 = [t158 / 0.2e1, 0, 0, -t147 * t163, t147 * t164, (t148 + t149) * t158 * qJ(2), t147 ^ 2 / 0.2e1 + (t149 / 0.2e1 + t148 / 0.2e1) * qJ(2) ^ 2 * t158, t157 ^ 2 * t166 / 0.2e1, -t157 * t162, -t157 * t144 * t164, t144 * t160, t144 ^ 2 / 0.2e1, qJ(2) * t162 - (-t155 * t161 + t137) * t144, qJ(2) * t157 * t166 + t165 * t144, -t123 * t135 + t124 * t134, t124 ^ 2 / 0.2e1 + t123 ^ 2 / 0.2e1 + t139 ^ 2 / 0.2e1, t126 ^ 2 / 0.2e1, -t126 * t125, -t126 * t141, t125 * t141, t141 ^ 2 / 0.2e1, -(t156 * t121 - t154 * t122) * t141 + t131 * t125, (t154 * t121 + t156 * t122) * t141 + t131 * t126;];
T_reg = t1;
