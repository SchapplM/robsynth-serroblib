% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRPR6
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
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRPR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:36:46
% EndTime: 2019-03-08 23:36:46
% DurationCPUTime: 0.15s
% Computational Cost: add. (235->45), mult. (515->95), div. (0->0), fcn. (357->10), ass. (0->42)
t210 = pkin(4) + pkin(5);
t194 = qJD(2) ^ 2;
t209 = t194 / 0.2e1;
t189 = sin(qJ(2));
t205 = qJD(1) * sin(pkin(6));
t175 = qJD(2) * pkin(8) + t189 * t205;
t188 = sin(qJ(3));
t192 = cos(qJ(3));
t204 = qJD(1) * cos(pkin(6));
t206 = t192 * t175 + t188 * t204;
t167 = qJD(3) * pkin(9) + t206;
t193 = cos(qJ(2));
t200 = t193 * t205;
t169 = -t200 + (-pkin(3) * t192 - pkin(9) * t188 - pkin(2)) * qJD(2);
t187 = sin(qJ(4));
t191 = cos(qJ(4));
t208 = t191 * t167 + t187 * t169;
t207 = -t188 * t175 + t192 * t204;
t203 = qJD(2) * t188;
t202 = t192 * qJD(2);
t201 = qJD(2) * qJD(3);
t181 = -qJD(4) + t202;
t160 = -t181 * qJ(5) + t208;
t199 = qJD(2) * t205;
t198 = -t187 * t167 + t191 * t169;
t197 = qJD(3) * pkin(3) + t207;
t196 = qJD(5) - t198;
t174 = t187 * qJD(3) + t191 * t203;
t195 = t174 * qJ(5) + t197;
t190 = cos(qJ(6));
t186 = sin(qJ(6));
t178 = qJD(6) + t181;
t176 = -qJD(2) * pkin(2) - t200;
t173 = -t191 * qJD(3) + t187 * t203;
t163 = t186 * t173 + t190 * t174;
t162 = -t190 * t173 + t186 * t174;
t161 = t173 * pkin(4) - t195;
t159 = t181 * pkin(4) + t196;
t158 = -t210 * t173 + t195;
t157 = t173 * pkin(10) + t160;
t156 = -t174 * pkin(10) + t210 * t181 + t196;
t1 = [qJD(1) ^ 2 / 0.2e1, t209, t193 * t199, -t189 * t199, t188 ^ 2 * t209, t188 * t194 * t192, t188 * t201, t192 * t201, qJD(3) ^ 2 / 0.2e1, t207 * qJD(3) - t176 * t202, -t206 * qJD(3) + t176 * t203, t174 ^ 2 / 0.2e1, -t174 * t173, -t174 * t181, t173 * t181, t181 ^ 2 / 0.2e1, -t173 * t197 - t198 * t181, -t174 * t197 + t208 * t181, t159 * t181 + t161 * t173, t159 * t174 - t160 * t173, -t160 * t181 - t161 * t174, t160 ^ 2 / 0.2e1 + t161 ^ 2 / 0.2e1 + t159 ^ 2 / 0.2e1, t163 ^ 2 / 0.2e1, -t163 * t162, t163 * t178, -t162 * t178, t178 ^ 2 / 0.2e1 (t190 * t156 - t186 * t157) * t178 + t158 * t162 -(t186 * t156 + t190 * t157) * t178 + t158 * t163;];
T_reg  = t1;
