% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRR9
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
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:25:18
% EndTime: 2019-03-09 07:25:18
% DurationCPUTime: 0.12s
% Computational Cost: add. (308->47), mult. (619->98), div. (0->0), fcn. (423->8), ass. (0->39)
t163 = qJD(1) ^ 2;
t173 = t163 / 0.2e1;
t172 = cos(qJ(5));
t171 = t163 * qJ(2);
t158 = sin(qJ(4));
t161 = cos(qJ(4));
t162 = cos(qJ(3));
t168 = qJD(1) * t162;
t149 = t158 * qJD(3) + t161 * t168;
t159 = sin(qJ(3));
t153 = t159 * qJD(1) + qJD(4);
t145 = (pkin(3) * t159 - pkin(8) * t162 + qJ(2)) * qJD(1);
t152 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t146 = qJD(3) * pkin(8) + t159 * t152;
t164 = t161 * t145 - t158 * t146;
t134 = t153 * pkin(4) - t149 * pkin(9) + t164;
t148 = -t161 * qJD(3) + t158 * t168;
t169 = t158 * t145 + t161 * t146;
t136 = -t148 * pkin(9) + t169;
t157 = sin(qJ(5));
t170 = t157 * t134 + t172 * t136;
t167 = qJD(3) * t152;
t166 = qJD(1) * qJD(3);
t165 = t172 * t134 - t157 * t136;
t147 = -qJD(3) * pkin(3) - t162 * t152;
t151 = qJD(5) + t153;
t140 = t148 * pkin(4) + t147;
t160 = cos(qJ(6));
t156 = sin(qJ(6));
t154 = -qJD(1) * pkin(1) + qJD(2);
t150 = qJD(6) + t151;
t139 = -t157 * t148 + t172 * t149;
t138 = t172 * t148 + t157 * t149;
t131 = t138 * pkin(5) + t140;
t130 = -t156 * t138 + t160 * t139;
t129 = t160 * t138 + t156 * t139;
t128 = -t138 * pkin(10) + t170;
t127 = t151 * pkin(5) - t139 * pkin(10) + t165;
t1 = [t173, 0, 0, t154 * qJD(1), t171, qJ(2) ^ 2 * t173 + t154 ^ 2 / 0.2e1, t162 ^ 2 * t173, -t162 * t163 * t159, t162 * t166, -t159 * t166, qJD(3) ^ 2 / 0.2e1, t159 * t171 + t162 * t167, -t159 * t167 + t162 * t171, t149 ^ 2 / 0.2e1, -t149 * t148, t149 * t153, -t148 * t153, t153 ^ 2 / 0.2e1, t147 * t148 + t164 * t153, t147 * t149 - t169 * t153, t139 ^ 2 / 0.2e1, -t139 * t138, t139 * t151, -t138 * t151, t151 ^ 2 / 0.2e1, t140 * t138 + t165 * t151, t140 * t139 - t170 * t151, t130 ^ 2 / 0.2e1, -t130 * t129, t130 * t150, -t129 * t150, t150 ^ 2 / 0.2e1 (t160 * t127 - t156 * t128) * t150 + t131 * t129 -(t156 * t127 + t160 * t128) * t150 + t131 * t130;];
T_reg  = t1;
