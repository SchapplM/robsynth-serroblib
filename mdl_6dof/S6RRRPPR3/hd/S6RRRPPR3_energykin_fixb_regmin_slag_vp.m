% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:30:06
% EndTime: 2019-03-09 15:30:07
% DurationCPUTime: 0.09s
% Computational Cost: add. (288->50), mult. (627->98), div. (0->0), fcn. (393->6), ass. (0->40)
t177 = -pkin(4) - pkin(9);
t176 = -pkin(8) - pkin(7);
t163 = qJD(1) ^ 2;
t175 = t163 / 0.2e1;
t162 = cos(qJ(2));
t174 = t162 * t163;
t159 = sin(qJ(2));
t172 = qJD(1) * t159;
t148 = qJD(2) * pkin(2) + t176 * t172;
t171 = qJD(1) * t162;
t149 = t176 * t171;
t158 = sin(qJ(3));
t161 = cos(qJ(3));
t173 = t158 * t148 - t161 * t149;
t150 = -qJD(1) * pkin(1) - pkin(2) * t171;
t170 = qJD(1) * qJD(2);
t155 = qJD(2) + qJD(3);
t137 = t155 * qJ(4) + t173;
t169 = t159 * t170;
t168 = t162 * t170;
t167 = t161 * t148 + t158 * t149;
t144 = t158 * t172 - t161 * t171;
t145 = (t158 * t162 + t159 * t161) * qJD(1);
t135 = t144 * pkin(3) - t145 * qJ(4) + t150;
t166 = qJD(4) - t167;
t165 = qJD(5) - t135;
t134 = -t144 * qJ(5) - t137;
t164 = -t145 * qJ(5) + t166;
t160 = cos(qJ(6));
t157 = sin(qJ(6));
t143 = qJD(6) + t145;
t139 = t160 * t144 - t157 * t155;
t138 = t157 * t144 + t160 * t155;
t136 = -t155 * pkin(3) + t166;
t133 = t155 * pkin(5) - t134;
t132 = -t144 * pkin(4) + t165;
t131 = (-pkin(3) - pkin(4)) * t155 + t164;
t130 = (-pkin(3) + t177) * t155 + t164;
t129 = t145 * pkin(5) + t177 * t144 + t165;
t1 = [t175, 0, 0, t159 ^ 2 * t175, t159 * t174, t169, t168, qJD(2) ^ 2 / 0.2e1, pkin(1) * t174 - pkin(7) * t169, -t163 * pkin(1) * t159 - pkin(7) * t168, t145 ^ 2 / 0.2e1, -t145 * t144, t145 * t155, -t144 * t155, t155 ^ 2 / 0.2e1, t150 * t144 + t167 * t155, t150 * t145 - t173 * t155, t135 * t144 - t136 * t155, t136 * t145 - t137 * t144, -t135 * t145 + t137 * t155, t137 ^ 2 / 0.2e1 + t135 ^ 2 / 0.2e1 + t136 ^ 2 / 0.2e1, t132 * t145 - t134 * t155, t131 * t155 + t132 * t144, -t131 * t145 - t134 * t144, t131 ^ 2 / 0.2e1 + t134 ^ 2 / 0.2e1 + t132 ^ 2 / 0.2e1, t139 ^ 2 / 0.2e1, -t139 * t138, t139 * t143, -t138 * t143, t143 ^ 2 / 0.2e1 (t160 * t129 - t157 * t130) * t143 + t133 * t138 -(t157 * t129 + t160 * t130) * t143 + t133 * t139;];
T_reg  = t1;
