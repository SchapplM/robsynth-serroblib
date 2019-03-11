% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:14:00
% EndTime: 2019-03-08 23:14:00
% DurationCPUTime: 0.12s
% Computational Cost: add. (240->46), mult. (546->97), div. (0->0), fcn. (388->10), ass. (0->42)
t196 = pkin(4) + pkin(10);
t182 = qJD(2) ^ 2;
t195 = t182 / 0.2e1;
t177 = sin(qJ(2));
t192 = qJD(1) * sin(pkin(6));
t166 = qJD(2) * pkin(8) + t177 * t192;
t180 = cos(qJ(3));
t191 = qJD(1) * cos(pkin(6));
t169 = t180 * t191;
t176 = sin(qJ(3));
t156 = qJD(3) * pkin(3) + t169 + (-pkin(9) * qJD(2) - t166) * t176;
t189 = qJD(2) * t180;
t193 = t180 * t166 + t176 * t191;
t157 = pkin(9) * t189 + t193;
t175 = sin(qJ(4));
t179 = cos(qJ(4));
t194 = t175 * t156 + t179 * t157;
t190 = qJD(2) * t176;
t188 = qJD(2) * qJD(3);
t181 = cos(qJ(2));
t187 = t181 * t192;
t186 = qJD(2) * t192;
t185 = t179 * t156 - t175 * t157;
t171 = qJD(3) + qJD(4);
t150 = -t171 * qJ(5) - t194;
t184 = qJD(5) - t185;
t164 = (t175 * t180 + t176 * t179) * qJD(2);
t161 = -t187 + (-pkin(3) * t180 - pkin(2)) * qJD(2);
t183 = -t164 * qJ(5) + t161;
t178 = cos(qJ(6));
t174 = sin(qJ(6));
t167 = -qJD(2) * pkin(2) - t187;
t163 = t175 * t190 - t179 * t189;
t162 = qJD(6) + t164;
t159 = t174 * t163 + t178 * t171;
t158 = -t178 * t163 + t174 * t171;
t152 = t163 * pkin(4) + t183;
t151 = t196 * t163 + t183;
t149 = -t171 * pkin(4) + t184;
t148 = -t163 * pkin(5) - t150;
t147 = t164 * pkin(5) - t196 * t171 + t184;
t1 = [qJD(1) ^ 2 / 0.2e1, t195, t181 * t186, -t177 * t186, t176 ^ 2 * t195, t176 * t182 * t180, t176 * t188, t180 * t188, qJD(3) ^ 2 / 0.2e1 (-t176 * t166 + t169) * qJD(3) - t167 * t189, -t193 * qJD(3) + t167 * t190, t164 ^ 2 / 0.2e1, -t164 * t163, t164 * t171, -t163 * t171, t171 ^ 2 / 0.2e1, t161 * t163 + t185 * t171, t161 * t164 - t194 * t171, t149 * t164 + t150 * t163, t149 * t171 - t152 * t163, -t150 * t171 - t152 * t164, t152 ^ 2 / 0.2e1 + t150 ^ 2 / 0.2e1 + t149 ^ 2 / 0.2e1, t159 ^ 2 / 0.2e1, -t159 * t158, t159 * t162, -t158 * t162, t162 ^ 2 / 0.2e1 (t178 * t147 - t174 * t151) * t162 + t148 * t158 -(t174 * t147 + t178 * t151) * t162 + t148 * t159;];
T_reg  = t1;
