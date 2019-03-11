% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRPR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:53:46
% EndTime: 2019-03-08 19:53:47
% DurationCPUTime: 0.08s
% Computational Cost: add. (123->42), mult. (257->83), div. (0->0), fcn. (144->8), ass. (0->37)
t162 = qJD(2) ^ 2;
t178 = t162 / 0.2e1;
t177 = qJD(1) ^ 2 / 0.2e1;
t161 = cos(qJ(2));
t175 = qJD(1) * sin(pkin(6));
t164 = -t161 * t175 + qJD(3);
t143 = (-pkin(2) - pkin(8)) * qJD(2) + t164;
t157 = sin(qJ(4));
t160 = cos(qJ(4));
t155 = cos(pkin(6));
t174 = qJD(1) * t155;
t176 = t157 * t143 + t160 * t174;
t173 = qJD(2) * t157;
t158 = sin(qJ(2));
t169 = t158 * t175;
t147 = qJD(2) * qJ(3) + t169;
t172 = t147 * qJD(2);
t171 = t160 * qJD(2);
t170 = qJD(2) * qJD(4);
t168 = qJD(2) * t175;
t167 = -qJ(5) * t160 + qJ(3);
t149 = t157 * t174;
t166 = t160 * t143 - t149;
t165 = pkin(4) * t173 + t169;
t139 = -qJD(4) * qJ(5) - t176;
t159 = cos(qJ(6));
t156 = sin(qJ(6));
t152 = qJD(6) + t171;
t146 = t159 * qJD(4) + t156 * t173;
t145 = t156 * qJD(4) - t159 * t173;
t144 = -qJD(2) * pkin(2) + t164;
t141 = t167 * qJD(2) + t165;
t140 = (pkin(9) * t157 + t167) * qJD(2) + t165;
t138 = -qJD(4) * pkin(4) + qJD(5) - t166;
t137 = -pkin(5) * t173 - t139;
t136 = qJD(5) + t149 + (pkin(5) * qJD(2) - t143) * t160 + (-pkin(4) - pkin(9)) * qJD(4);
t1 = [t177, t178, t161 * t168, -t158 * t168, t144 * qJD(2), t172, t155 ^ 2 * t177 + t147 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1, t160 ^ 2 * t178, -t160 * t162 * t157, t160 * t170, -t157 * t170, qJD(4) ^ 2 / 0.2e1, t166 * qJD(4) + t157 * t172, -t176 * qJD(4) + t147 * t171 (t138 * t160 + t139 * t157) * qJD(2), t138 * qJD(4) - t141 * t173, -t139 * qJD(4) - t141 * t171, t141 ^ 2 / 0.2e1 + t139 ^ 2 / 0.2e1 + t138 ^ 2 / 0.2e1, t146 ^ 2 / 0.2e1, -t146 * t145, t146 * t152, -t145 * t152, t152 ^ 2 / 0.2e1 (t159 * t136 - t156 * t140) * t152 + t137 * t145 -(t156 * t136 + t159 * t140) * t152 + t137 * t146;];
T_reg  = t1;
