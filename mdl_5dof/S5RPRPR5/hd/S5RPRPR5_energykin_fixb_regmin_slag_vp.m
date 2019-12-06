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
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:57:11
% EndTime: 2019-12-05 17:57:11
% DurationCPUTime: 0.12s
% Computational Cost: add. (190->41), mult. (508->92), div. (0->0), fcn. (336->8), ass. (0->37)
t159 = sin(qJ(3));
t171 = qJ(2) * t159;
t155 = sin(pkin(8));
t152 = t155 ^ 2;
t162 = qJD(1) ^ 2;
t170 = t152 * t162;
t157 = cos(pkin(8));
t142 = qJD(2) + (-pkin(2) * t157 - pkin(6) * t155 - pkin(1)) * qJD(1);
t161 = cos(qJ(3));
t141 = t161 * t142;
t167 = t157 * qJD(1);
t148 = -qJD(3) + t167;
t133 = -t148 * pkin(3) + t141 + (-qJ(4) * t155 * t161 - t157 * t171) * qJD(1);
t168 = qJD(1) * t155;
t164 = t159 * t168;
t165 = qJ(2) * t167;
t169 = t159 * t142 + t161 * t165;
t136 = -qJ(4) * t164 + t169;
t154 = sin(pkin(9));
t156 = cos(pkin(9));
t128 = t154 * t133 + t156 * t136;
t166 = t161 * t170;
t143 = pkin(3) * t164 + qJ(2) * t168 + qJD(4);
t127 = t156 * t133 - t154 * t136;
t160 = cos(qJ(5));
t158 = sin(qJ(5));
t153 = t157 ^ 2;
t151 = -qJD(1) * pkin(1) + qJD(2);
t145 = -qJD(5) + t148;
t139 = (-t154 * t159 + t156 * t161) * t168;
t138 = (-t154 * t161 - t156 * t159) * t168;
t135 = -t138 * pkin(4) + t143;
t130 = t158 * t138 + t160 * t139;
t129 = -t160 * t138 + t158 * t139;
t126 = t138 * pkin(7) + t128;
t125 = -t148 * pkin(4) - t139 * pkin(7) + t127;
t1 = [t162 / 0.2e1, 0, 0, -t151 * t167, t151 * t168, (t152 + t153) * t162 * qJ(2), t151 ^ 2 / 0.2e1 + (t153 / 0.2e1 + t152 / 0.2e1) * qJ(2) ^ 2 * t162, t161 ^ 2 * t170 / 0.2e1, -t159 * t166, -t161 * t148 * t168, t148 * t164, t148 ^ 2 / 0.2e1, t170 * t171 - (-t159 * t165 + t141) * t148, qJ(2) * t166 + t169 * t148, -t127 * t139 + t128 * t138, t128 ^ 2 / 0.2e1 + t127 ^ 2 / 0.2e1 + t143 ^ 2 / 0.2e1, t130 ^ 2 / 0.2e1, -t130 * t129, -t130 * t145, t129 * t145, t145 ^ 2 / 0.2e1, -(t160 * t125 - t158 * t126) * t145 + t135 * t129, (t158 * t125 + t160 * t126) * t145 + t135 * t130;];
T_reg = t1;
