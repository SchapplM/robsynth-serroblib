% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:54
% EndTime: 2020-01-03 11:22:54
% DurationCPUTime: 0.10s
% Computational Cost: add. (208->43), mult. (553->95), div. (0->0), fcn. (353->8), ass. (0->35)
t158 = sin(pkin(7));
t160 = cos(pkin(8));
t170 = t158 * t160;
t161 = cos(pkin(7));
t146 = qJD(2) + (-pkin(2) * t161 - qJ(3) * t158 - pkin(1)) * qJD(1);
t157 = sin(pkin(8));
t168 = qJD(1) * t161;
t167 = qJ(2) * t168;
t139 = t157 * t146 + t160 * t167;
t135 = -qJ(4) * t168 + t139;
t169 = qJD(1) * t158;
t150 = qJ(2) * t169 + qJD(3);
t141 = (pkin(3) * t157 - qJ(4) * t160) * t169 + t150;
t156 = sin(pkin(9));
t159 = cos(pkin(9));
t132 = t159 * t135 + t156 * t141;
t166 = t157 * t169;
t138 = t160 * t146 - t157 * t167;
t131 = -t156 * t135 + t159 * t141;
t134 = pkin(3) * t168 + qJD(4) - t138;
t144 = (t156 * t170 + t159 * t161) * qJD(1);
t164 = qJD(1) ^ 2;
t163 = cos(qJ(5));
t162 = sin(qJ(5));
t155 = t161 ^ 2;
t154 = t158 ^ 2;
t153 = -qJD(1) * pkin(1) + qJD(2);
t145 = (-t156 * t161 + t159 * t170) * qJD(1);
t142 = qJD(5) + t144;
t137 = t163 * t145 + t162 * t166;
t136 = t162 * t145 - t163 * t166;
t130 = pkin(6) * t166 + t132;
t129 = -pkin(4) * t166 - t131;
t128 = t144 * pkin(4) - t145 * pkin(6) + t134;
t1 = [t164 / 0.2e1, 0, 0, -t153 * t168, t153 * t169, (t154 + t155) * t164 * qJ(2), t153 ^ 2 / 0.2e1 + (t155 / 0.2e1 + t154 / 0.2e1) * qJ(2) ^ 2 * t164, (t150 * t157 * t158 - t138 * t161) * qJD(1), (t139 * t161 + t150 * t170) * qJD(1), (-t138 * t160 - t139 * t157) * t169, t139 ^ 2 / 0.2e1 + t138 ^ 2 / 0.2e1 + t150 ^ 2 / 0.2e1, t131 * t166 + t134 * t144, -t132 * t166 + t134 * t145, -t131 * t145 - t132 * t144, t132 ^ 2 / 0.2e1 + t131 ^ 2 / 0.2e1 + t134 ^ 2 / 0.2e1, t137 ^ 2 / 0.2e1, -t137 * t136, t137 * t142, -t136 * t142, t142 ^ 2 / 0.2e1, (t163 * t128 - t162 * t130) * t142 + t129 * t136, -(t162 * t128 + t163 * t130) * t142 + t129 * t137;];
T_reg = t1;
