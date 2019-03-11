% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRPR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:49:33
% EndTime: 2019-03-08 19:49:33
% DurationCPUTime: 0.08s
% Computational Cost: add. (197->42), mult. (411->88), div. (0->0), fcn. (268->10), ass. (0->40)
t178 = qJD(2) ^ 2;
t193 = t178 / 0.2e1;
t192 = qJD(1) ^ 2 / 0.2e1;
t191 = cos(pkin(11));
t177 = cos(qJ(2));
t189 = qJD(1) * sin(pkin(6));
t180 = -t177 * t189 + qJD(3);
t159 = (-pkin(2) - pkin(8)) * qJD(2) + t180;
t173 = sin(qJ(4));
t176 = cos(qJ(4));
t171 = cos(pkin(6));
t188 = qJD(1) * t171;
t190 = t173 * t159 + t176 * t188;
t153 = qJD(4) * qJ(5) + t190;
t174 = sin(qJ(2));
t183 = t174 * t189;
t156 = t183 + (pkin(4) * t173 - qJ(5) * t176 + qJ(3)) * qJD(2);
t169 = sin(pkin(11));
t147 = t191 * t153 + t169 * t156;
t187 = qJD(2) * t176;
t163 = qJD(2) * qJ(3) + t183;
t186 = t163 * qJD(2);
t185 = t173 * qJD(2);
t184 = qJD(2) * qJD(4);
t182 = qJD(2) * t189;
t146 = -t169 * t153 + t191 * t156;
t181 = t176 * t159 - t173 * t188;
t152 = -qJD(4) * pkin(4) + qJD(5) - t181;
t175 = cos(qJ(6));
t172 = sin(qJ(6));
t167 = qJD(6) + t185;
t162 = -qJD(2) * pkin(2) + t180;
t161 = t169 * qJD(4) + t191 * t187;
t160 = -t191 * qJD(4) + t169 * t187;
t151 = -t172 * t160 + t175 * t161;
t150 = t175 * t160 + t172 * t161;
t148 = t160 * pkin(5) + t152;
t145 = -t160 * pkin(9) + t147;
t144 = pkin(5) * t185 - t161 * pkin(9) + t146;
t1 = [t192, t193, t177 * t182, -t174 * t182, t162 * qJD(2), t186, t171 ^ 2 * t192 + t163 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1, t176 ^ 2 * t193, -t176 * t178 * t173, t176 * t184, -t173 * t184, qJD(4) ^ 2 / 0.2e1, t181 * qJD(4) + t163 * t185, -t190 * qJD(4) + t176 * t186, t146 * t185 + t152 * t160, -t147 * t185 + t152 * t161, -t146 * t161 - t147 * t160, t147 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1 + t152 ^ 2 / 0.2e1, t151 ^ 2 / 0.2e1, -t151 * t150, t151 * t167, -t150 * t167, t167 ^ 2 / 0.2e1 (t175 * t144 - t172 * t145) * t167 + t148 * t150 -(t172 * t144 + t175 * t145) * t167 + t148 * t151;];
T_reg  = t1;
