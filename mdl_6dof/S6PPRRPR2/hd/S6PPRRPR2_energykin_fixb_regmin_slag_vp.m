% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PPRRPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:50:54
% EndTime: 2019-03-08 18:50:54
% DurationCPUTime: 0.14s
% Computational Cost: add. (178->40), mult. (454->90), div. (0->0), fcn. (356->12), ass. (0->42)
t180 = cos(pkin(6)) * qJD(1) + qJD(2);
t184 = sin(pkin(7));
t187 = cos(pkin(7));
t186 = cos(pkin(12));
t185 = sin(pkin(6));
t208 = qJD(1) * t185;
t202 = t186 * t208;
t214 = t180 * t184 + t187 * t202;
t190 = sin(qJ(3));
t193 = cos(qJ(3));
t183 = sin(pkin(12));
t203 = t183 * t208;
t213 = -t190 * t203 + t214 * t193;
t212 = -pkin(4) - pkin(10);
t194 = qJD(3) ^ 2;
t211 = t194 / 0.2e1;
t204 = t214 * t190 + t193 * t203;
t171 = qJD(3) * pkin(9) + t204;
t173 = t187 * t180 - t184 * t202;
t189 = sin(qJ(4));
t192 = cos(qJ(4));
t209 = t192 * t171 + t189 * t173;
t207 = qJD(3) * t192;
t206 = t189 * qJD(3);
t205 = qJD(3) * qJD(4);
t201 = -qJ(5) * t189 - pkin(3);
t200 = -t189 * t171 + t192 * t173;
t165 = -qJD(4) * qJ(5) - t209;
t198 = qJD(5) - t200;
t195 = qJD(1) ^ 2;
t191 = cos(qJ(6));
t188 = sin(qJ(6));
t181 = qJD(6) + t206;
t177 = t191 * qJD(4) - t188 * t207;
t176 = t188 * qJD(4) + t191 * t207;
t170 = -qJD(3) * pkin(3) - t213;
t167 = (-pkin(4) * t192 + t201) * qJD(3) - t213;
t166 = (t212 * t192 + t201) * qJD(3) - t213;
t164 = -qJD(4) * pkin(4) + t198;
t163 = pkin(5) * t207 - t165;
t162 = pkin(5) * t206 + t212 * qJD(4) + t198;
t1 = [t195 / 0.2e1, t180 ^ 2 / 0.2e1 + (t183 ^ 2 / 0.2e1 + t186 ^ 2 / 0.2e1) * t195 * t185 ^ 2, t211, t213 * qJD(3), -t204 * qJD(3), t189 ^ 2 * t211, t189 * t194 * t192, t189 * t205, t192 * t205, qJD(4) ^ 2 / 0.2e1, t200 * qJD(4) - t170 * t207, -t209 * qJD(4) + t170 * t206 (t164 * t189 - t165 * t192) * qJD(3), t164 * qJD(4) + t167 * t207, -t165 * qJD(4) - t167 * t206, t167 ^ 2 / 0.2e1 + t165 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1, t177 ^ 2 / 0.2e1, -t177 * t176, t177 * t181, -t176 * t181, t181 ^ 2 / 0.2e1 (t191 * t162 - t188 * t166) * t181 + t163 * t176 -(t188 * t162 + t191 * t166) * t181 + t163 * t177;];
T_reg  = t1;
