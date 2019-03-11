% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:38:55
% EndTime: 2019-03-08 20:38:55
% DurationCPUTime: 0.13s
% Computational Cost: add. (280->49), mult. (692->101), div. (0->0), fcn. (542->12), ass. (0->43)
t207 = cos(qJ(5));
t193 = sin(qJ(2));
t204 = qJD(1) * sin(pkin(6));
t181 = qJD(2) * qJ(3) + t193 * t204;
t188 = cos(pkin(12));
t203 = qJD(1) * cos(pkin(6));
t183 = t188 * t203;
t186 = sin(pkin(12));
t167 = t183 + (-pkin(8) * qJD(2) - t181) * t186;
t173 = t188 * t181 + t186 * t203;
t201 = qJD(2) * t188;
t168 = pkin(8) * t201 + t173;
t192 = sin(qJ(4));
t195 = cos(qJ(4));
t205 = t192 * t167 + t195 * t168;
t158 = qJD(4) * pkin(9) + t205;
t196 = cos(qJ(2));
t197 = -t196 * t204 + qJD(3);
t175 = (-pkin(3) * t188 - pkin(2)) * qJD(2) + t197;
t202 = qJD(2) * t186;
t177 = t192 * t202 - t195 * t201;
t178 = (t186 * t195 + t188 * t192) * qJD(2);
t163 = t177 * pkin(4) - t178 * pkin(9) + t175;
t191 = sin(qJ(5));
t206 = t207 * t158 + t191 * t163;
t200 = qJD(2) * t204;
t199 = -t191 * t158 + t207 * t163;
t198 = t195 * t167 - t192 * t168;
t176 = qJD(5) + t177;
t157 = -qJD(4) * pkin(4) - t198;
t194 = cos(qJ(6));
t190 = sin(qJ(6));
t180 = -qJD(2) * pkin(2) + t197;
t174 = qJD(6) + t176;
t172 = -t186 * t181 + t183;
t171 = t191 * qJD(4) + t207 * t178;
t170 = -t207 * qJD(4) + t191 * t178;
t162 = -t190 * t170 + t194 * t171;
t161 = t194 * t170 + t190 * t171;
t155 = t170 * pkin(5) + t157;
t154 = -t170 * pkin(10) + t206;
t153 = t176 * pkin(5) - t171 * pkin(10) + t199;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t196 * t200, -t193 * t200, -t180 * t201, t180 * t202 (-t172 * t186 + t173 * t188) * qJD(2), t173 ^ 2 / 0.2e1 + t172 ^ 2 / 0.2e1 + t180 ^ 2 / 0.2e1, t178 ^ 2 / 0.2e1, -t178 * t177, t178 * qJD(4), -t177 * qJD(4), qJD(4) ^ 2 / 0.2e1, t198 * qJD(4) + t175 * t177, -t205 * qJD(4) + t175 * t178, t171 ^ 2 / 0.2e1, -t171 * t170, t171 * t176, -t170 * t176, t176 ^ 2 / 0.2e1, t157 * t170 + t199 * t176, t157 * t171 - t206 * t176, t162 ^ 2 / 0.2e1, -t162 * t161, t162 * t174, -t161 * t174, t174 ^ 2 / 0.2e1 (t194 * t153 - t190 * t154) * t174 + t155 * t161 -(t190 * t153 + t194 * t154) * t174 + t155 * t162;];
T_reg  = t1;
