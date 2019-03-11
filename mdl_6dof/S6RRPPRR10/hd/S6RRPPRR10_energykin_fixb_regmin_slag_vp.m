% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:37:00
% EndTime: 2019-03-09 09:37:00
% DurationCPUTime: 0.13s
% Computational Cost: add. (359->53), mult. (771->109), div. (0->0), fcn. (491->8), ass. (0->44)
t194 = qJD(1) ^ 2;
t207 = t194 / 0.2e1;
t206 = cos(qJ(5));
t205 = -pkin(2) - qJ(4);
t193 = cos(qJ(2));
t204 = t193 * t194;
t191 = sin(qJ(2));
t196 = -qJ(3) * t191 - pkin(1);
t172 = (t205 * t193 + t196) * qJD(1);
t201 = t191 * qJD(1);
t200 = pkin(7) * t201 + qJD(3);
t173 = pkin(3) * t201 + t205 * qJD(2) + t200;
t187 = sin(pkin(10));
t188 = cos(pkin(10));
t163 = -t187 * t172 + t188 * t173;
t202 = qJD(1) * t193;
t178 = t188 * qJD(2) - t187 * t202;
t159 = pkin(4) * t201 - t178 * pkin(8) + t163;
t164 = t188 * t172 + t187 * t173;
t177 = t187 * qJD(2) + t188 * t202;
t162 = -t177 * pkin(8) + t164;
t190 = sin(qJ(5));
t203 = t190 * t159 + t206 * t162;
t180 = -pkin(7) * t202 - qJD(2) * qJ(3);
t199 = qJD(1) * qJD(2);
t198 = t191 * t199;
t197 = t193 * t199;
t195 = t206 * t159 - t190 * t162;
t182 = qJD(5) + t201;
t175 = pkin(3) * t202 + qJD(4) - t180;
t168 = t177 * pkin(4) + t175;
t192 = cos(qJ(6));
t189 = sin(qJ(6));
t181 = qJD(6) + t182;
t179 = -qJD(2) * pkin(2) + t200;
t176 = (-pkin(2) * t193 + t196) * qJD(1);
t167 = -t190 * t177 + t206 * t178;
t166 = t206 * t177 + t190 * t178;
t160 = t166 * pkin(5) + t168;
t156 = -t189 * t166 + t192 * t167;
t155 = t192 * t166 + t189 * t167;
t154 = -t166 * pkin(9) + t203;
t153 = t182 * pkin(5) - t167 * pkin(9) + t195;
t1 = [t207, 0, 0, t191 ^ 2 * t207, t191 * t204, t198, t197, qJD(2) ^ 2 / 0.2e1, pkin(1) * t204 - pkin(7) * t198, -t194 * pkin(1) * t191 - pkin(7) * t197 (t179 * t191 - t180 * t193) * qJD(1), t179 * qJD(2) + t176 * t202, -t180 * qJD(2) - t176 * t201, t176 ^ 2 / 0.2e1 + t180 ^ 2 / 0.2e1 + t179 ^ 2 / 0.2e1, t163 * t201 + t175 * t177, -t164 * t201 + t175 * t178, -t163 * t178 - t164 * t177, t164 ^ 2 / 0.2e1 + t163 ^ 2 / 0.2e1 + t175 ^ 2 / 0.2e1, t167 ^ 2 / 0.2e1, -t167 * t166, t167 * t182, -t166 * t182, t182 ^ 2 / 0.2e1, t168 * t166 + t195 * t182, t168 * t167 - t203 * t182, t156 ^ 2 / 0.2e1, -t156 * t155, t156 * t181, -t155 * t181, t181 ^ 2 / 0.2e1 (t192 * t153 - t189 * t154) * t181 + t160 * t155 -(t189 * t153 + t192 * t154) * t181 + t160 * t156;];
T_reg  = t1;
