% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% T_reg [1x34]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:21:12
% EndTime: 2019-03-09 07:21:12
% DurationCPUTime: 0.13s
% Computational Cost: add. (310->46), mult. (599->97), div. (0->0), fcn. (405->8), ass. (0->39)
t162 = qJD(1) ^ 2;
t172 = t162 / 0.2e1;
t171 = cos(qJ(5));
t170 = t162 * qJ(2);
t153 = qJD(3) + qJD(4);
t161 = cos(qJ(3));
t149 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t165 = -pkin(8) * qJD(1) + t149;
t142 = qJD(3) * pkin(3) + t165 * t161;
t158 = sin(qJ(3));
t144 = t165 * t158;
t157 = sin(qJ(4));
t160 = cos(qJ(4));
t168 = t157 * t142 + t160 * t144;
t132 = t153 * pkin(9) + t168;
t146 = (t157 * t161 + t158 * t160) * qJD(1);
t147 = (-t157 * t158 + t160 * t161) * qJD(1);
t148 = (pkin(3) * t158 + qJ(2)) * qJD(1);
t135 = t146 * pkin(4) - t147 * pkin(9) + t148;
t156 = sin(qJ(5));
t169 = t171 * t132 + t156 * t135;
t167 = qJD(3) * t149;
t166 = qJD(1) * qJD(3);
t164 = -t156 * t132 + t171 * t135;
t163 = t160 * t142 - t157 * t144;
t131 = -t153 * pkin(4) - t163;
t145 = qJD(5) + t146;
t159 = cos(qJ(6));
t155 = sin(qJ(6));
t152 = -qJD(1) * pkin(1) + qJD(2);
t143 = qJD(6) + t145;
t138 = t171 * t147 + t156 * t153;
t137 = t156 * t147 - t171 * t153;
t129 = -t155 * t137 + t159 * t138;
t128 = t159 * t137 + t155 * t138;
t127 = t137 * pkin(5) + t131;
t126 = -t137 * pkin(10) + t169;
t125 = t145 * pkin(5) - t138 * pkin(10) + t164;
t1 = [t172, 0, 0, t152 * qJD(1), t170, qJ(2) ^ 2 * t172 + t152 ^ 2 / 0.2e1, t161 ^ 2 * t172, -t161 * t162 * t158, t161 * t166, -t158 * t166, qJD(3) ^ 2 / 0.2e1, t158 * t170 + t161 * t167, -t158 * t167 + t161 * t170, t147 ^ 2 / 0.2e1, -t147 * t146, t147 * t153, -t146 * t153, t153 ^ 2 / 0.2e1, t148 * t146 + t163 * t153, t148 * t147 - t168 * t153, t138 ^ 2 / 0.2e1, -t138 * t137, t138 * t145, -t137 * t145, t145 ^ 2 / 0.2e1, t131 * t137 + t164 * t145, t131 * t138 - t169 * t145, t129 ^ 2 / 0.2e1, -t129 * t128, t129 * t143, -t128 * t143, t143 ^ 2 / 0.2e1 (t159 * t125 - t155 * t126) * t143 + t127 * t128 -(t155 * t125 + t159 * t126) * t143 + t127 * t129;];
T_reg  = t1;
