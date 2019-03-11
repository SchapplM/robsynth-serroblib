% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:19:44
% EndTime: 2019-03-09 03:19:45
% DurationCPUTime: 0.14s
% Computational Cost: add. (255->48), mult. (637->93), div. (0->0), fcn. (419->6), ass. (0->38)
t173 = pkin(3) + pkin(8);
t172 = cos(qJ(5));
t171 = pkin(7) + qJ(2);
t159 = sin(qJ(3));
t160 = cos(qJ(3));
t157 = cos(pkin(9));
t167 = qJD(1) * t157;
t156 = sin(pkin(9));
t168 = qJD(1) * t156;
t146 = t159 * t168 - t160 * t167;
t147 = (t156 * t160 + t157 * t159) * qJD(1);
t150 = qJD(2) + (-pkin(2) * t157 - pkin(1)) * qJD(1);
t163 = -t147 * qJ(4) + t150;
t131 = t173 * t146 + t163;
t148 = t171 * t168;
t149 = t171 * t167;
t165 = -t160 * t148 - t159 * t149;
t164 = qJD(4) - t165;
t134 = t147 * pkin(4) - t173 * qJD(3) + t164;
t158 = sin(qJ(5));
t170 = t172 * t131 + t158 * t134;
t169 = -t159 * t148 + t160 * t149;
t138 = -qJD(3) * qJ(4) - t169;
t166 = -t158 * t131 + t172 * t134;
t135 = -t146 * pkin(4) - t138;
t161 = qJD(1) ^ 2;
t154 = t157 ^ 2;
t153 = t156 ^ 2;
t152 = -qJD(1) * pkin(1) + qJD(2);
t142 = qJD(5) + t147;
t140 = t172 * qJD(3) + t158 * t146;
t139 = t158 * qJD(3) - t172 * t146;
t137 = -qJD(3) * pkin(3) + t164;
t136 = t146 * pkin(3) + t163;
t129 = t139 * pkin(5) + qJD(6) + t135;
t128 = -t139 * qJ(6) + t170;
t127 = t142 * pkin(5) - t140 * qJ(6) + t166;
t1 = [t161 / 0.2e1, 0, 0, -t152 * t167, t152 * t168 (t153 + t154) * t161 * qJ(2), t152 ^ 2 / 0.2e1 + (t154 / 0.2e1 + t153 / 0.2e1) * qJ(2) ^ 2 * t161, t147 ^ 2 / 0.2e1, -t147 * t146, t147 * qJD(3), -t146 * qJD(3), qJD(3) ^ 2 / 0.2e1, qJD(3) * t165 + t150 * t146, -t169 * qJD(3) + t150 * t147, t137 * t147 + t138 * t146, t137 * qJD(3) - t136 * t146, -t138 * qJD(3) - t136 * t147, t136 ^ 2 / 0.2e1 + t138 ^ 2 / 0.2e1 + t137 ^ 2 / 0.2e1, t140 ^ 2 / 0.2e1, -t140 * t139, t140 * t142, -t139 * t142, t142 ^ 2 / 0.2e1, t135 * t139 + t142 * t166, t135 * t140 - t170 * t142, -t127 * t140 - t128 * t139, t128 ^ 2 / 0.2e1 + t127 ^ 2 / 0.2e1 + t129 ^ 2 / 0.2e1;];
T_reg  = t1;
