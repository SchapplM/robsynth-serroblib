% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:10:56
% EndTime: 2019-03-08 21:10:56
% DurationCPUTime: 0.08s
% Computational Cost: add. (148->40), mult. (335->90), div. (0->0), fcn. (194->8), ass. (0->34)
t175 = -pkin(3) - pkin(4);
t163 = qJD(2) ^ 2;
t174 = t163 / 0.2e1;
t159 = sin(qJ(2));
t172 = qJD(1) * sin(pkin(6));
t145 = qJD(2) * pkin(8) + t159 * t172;
t158 = sin(qJ(3));
t161 = cos(qJ(3));
t171 = qJD(1) * cos(pkin(6));
t173 = t161 * t145 + t158 * t171;
t162 = cos(qJ(2));
t146 = -qJD(2) * pkin(2) - t162 * t172;
t170 = qJD(2) * t161;
t169 = t158 * qJD(2);
t168 = qJD(2) * qJD(3);
t139 = qJD(3) * qJ(4) + t173;
t167 = qJD(2) * t172;
t140 = -pkin(3) * t170 - qJ(4) * t169 + t146;
t166 = -t158 * t145 + t161 * t171;
t165 = qJD(4) - t166;
t137 = pkin(4) * t170 + qJD(5) - t140;
t136 = qJ(5) * t170 - t139;
t164 = -qJ(5) * t169 + t165;
t160 = cos(qJ(6));
t157 = sin(qJ(6));
t149 = qJD(6) + t169;
t144 = -t157 * qJD(3) - t160 * t170;
t143 = -t160 * qJD(3) + t157 * t170;
t138 = -qJD(3) * pkin(3) + t165;
t135 = qJD(3) * pkin(5) - t136;
t134 = t175 * qJD(3) + t164;
t133 = (-pkin(9) + t175) * qJD(3) + t164;
t132 = (pkin(5) * t158 + pkin(9) * t161) * qJD(2) + t137;
t1 = [qJD(1) ^ 2 / 0.2e1, t174, t162 * t167, -t159 * t167, t158 ^ 2 * t174, t158 * t163 * t161, t158 * t168, t161 * t168, qJD(3) ^ 2 / 0.2e1, t166 * qJD(3) - t146 * t170, -t173 * qJD(3) + t146 * t169, -t138 * qJD(3) - t140 * t170 (t138 * t158 + t139 * t161) * qJD(2), t139 * qJD(3) - t140 * t169, t139 ^ 2 / 0.2e1 + t140 ^ 2 / 0.2e1 + t138 ^ 2 / 0.2e1, -t136 * qJD(3) + t137 * t169, t134 * qJD(3) - t137 * t170 (-t134 * t158 + t136 * t161) * qJD(2), t134 ^ 2 / 0.2e1 + t136 ^ 2 / 0.2e1 + t137 ^ 2 / 0.2e1, t144 ^ 2 / 0.2e1, t144 * t143, t144 * t149, t143 * t149, t149 ^ 2 / 0.2e1 (t160 * t132 - t157 * t133) * t149 - t135 * t143 -(t157 * t132 + t160 * t133) * t149 + t135 * t144;];
T_reg  = t1;
