% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:39:50
% EndTime: 2019-03-09 00:39:50
% DurationCPUTime: 0.10s
% Computational Cost: add. (324->48), mult. (758->104), div. (0->0), fcn. (591->12), ass. (0->46)
t196 = qJD(2) ^ 2;
t210 = t196 / 0.2e1;
t209 = cos(qJ(4));
t189 = sin(qJ(4));
t190 = sin(qJ(3));
t194 = cos(qJ(3));
t175 = (t189 * t194 + t209 * t190) * qJD(2);
t184 = qJD(3) + qJD(4);
t191 = sin(qJ(2));
t205 = qJD(1) * sin(pkin(6));
t177 = qJD(2) * pkin(8) + t191 * t205;
t204 = qJD(1) * cos(pkin(6));
t181 = t194 * t204;
t170 = qJD(3) * pkin(3) + t181 + (-pkin(9) * qJD(2) - t177) * t190;
t202 = qJD(2) * t194;
t206 = t194 * t177 + t190 * t204;
t171 = pkin(9) * t202 + t206;
t198 = t209 * t170 - t189 * t171;
t158 = t184 * pkin(4) - t175 * pkin(10) + t198;
t203 = qJD(2) * t190;
t174 = t189 * t203 - t209 * t202;
t207 = t189 * t170 + t209 * t171;
t160 = -t174 * pkin(10) + t207;
t188 = sin(qJ(5));
t193 = cos(qJ(5));
t208 = t188 * t158 + t193 * t160;
t201 = qJD(2) * qJD(3);
t195 = cos(qJ(2));
t200 = t195 * t205;
t199 = qJD(2) * t205;
t164 = t193 * t174 + t188 * t175;
t197 = t193 * t158 - t188 * t160;
t173 = -t200 + (-pkin(3) * t194 - pkin(2)) * qJD(2);
t168 = t174 * pkin(4) + t173;
t192 = cos(qJ(6));
t187 = sin(qJ(6));
t183 = qJD(5) + t184;
t178 = -qJD(2) * pkin(2) - t200;
t165 = -t188 * t174 + t193 * t175;
t163 = qJD(6) + t164;
t162 = t192 * t165 + t187 * t183;
t161 = t187 * t165 - t192 * t183;
t156 = t164 * pkin(5) - t165 * pkin(11) + t168;
t155 = t183 * pkin(11) + t208;
t154 = -t183 * pkin(5) - t197;
t1 = [qJD(1) ^ 2 / 0.2e1, t210, t195 * t199, -t191 * t199, t190 ^ 2 * t210, t190 * t196 * t194, t190 * t201, t194 * t201, qJD(3) ^ 2 / 0.2e1 (-t190 * t177 + t181) * qJD(3) - t178 * t202, -t206 * qJD(3) + t178 * t203, t175 ^ 2 / 0.2e1, -t175 * t174, t175 * t184, -t174 * t184, t184 ^ 2 / 0.2e1, t173 * t174 + t198 * t184, t173 * t175 - t207 * t184, t165 ^ 2 / 0.2e1, -t165 * t164, t165 * t183, -t164 * t183, t183 ^ 2 / 0.2e1, t168 * t164 + t197 * t183, t168 * t165 - t208 * t183, t162 ^ 2 / 0.2e1, -t162 * t161, t162 * t163, -t161 * t163, t163 ^ 2 / 0.2e1 (-t187 * t155 + t192 * t156) * t163 + t154 * t161 -(t192 * t155 + t187 * t156) * t163 + t154 * t162;];
T_reg  = t1;
