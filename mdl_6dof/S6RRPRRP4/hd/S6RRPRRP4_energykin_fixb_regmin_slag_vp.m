% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:55:38
% EndTime: 2019-03-09 11:55:39
% DurationCPUTime: 0.12s
% Computational Cost: add. (543->49), mult. (1261->102), div. (0->0), fcn. (916->8), ass. (0->44)
t194 = qJD(1) ^ 2;
t207 = t194 / 0.2e1;
t206 = cos(qJ(4));
t205 = pkin(7) + qJ(3);
t193 = cos(qJ(2));
t204 = t193 * t194;
t187 = sin(pkin(10));
t188 = cos(pkin(10));
t191 = sin(qJ(2));
t179 = (t187 * t193 + t188 * t191) * qJD(1);
t190 = sin(qJ(4));
t174 = t190 * qJD(2) + t206 * t179;
t200 = qJD(1) * t193;
t201 = qJD(1) * t191;
t178 = -t187 * t201 + t188 * t200;
t177 = qJD(4) - t178;
t184 = qJD(3) + (-pkin(2) * t193 - pkin(1)) * qJD(1);
t166 = -t178 * pkin(3) - t179 * pkin(8) + t184;
t182 = qJD(2) * pkin(2) - t205 * t201;
t183 = t205 * t200;
t171 = t187 * t182 + t188 * t183;
t169 = qJD(2) * pkin(8) + t171;
t196 = t206 * t166 - t190 * t169;
t158 = t177 * pkin(4) - t174 * pkin(9) + t196;
t173 = -t206 * qJD(2) + t190 * t179;
t202 = t190 * t166 + t206 * t169;
t160 = -t173 * pkin(9) + t202;
t189 = sin(qJ(5));
t192 = cos(qJ(5));
t203 = t189 * t158 + t192 * t160;
t199 = qJD(1) * qJD(2);
t198 = t191 * t199;
t197 = t193 * t199;
t170 = t188 * t182 - t187 * t183;
t195 = t192 * t158 - t189 * t160;
t168 = -qJD(2) * pkin(3) - t170;
t161 = t173 * pkin(4) + t168;
t175 = qJD(5) + t177;
t163 = -t189 * t173 + t192 * t174;
t162 = t192 * t173 + t189 * t174;
t156 = t162 * pkin(5) - t163 * qJ(6) + t161;
t155 = t175 * qJ(6) + t203;
t154 = -t175 * pkin(5) + qJD(6) - t195;
t1 = [t207, 0, 0, t191 ^ 2 * t207, t191 * t204, t198, t197, qJD(2) ^ 2 / 0.2e1, pkin(1) * t204 - pkin(7) * t198, -t194 * pkin(1) * t191 - pkin(7) * t197, -t170 * t179 + t171 * t178, t171 ^ 2 / 0.2e1 + t170 ^ 2 / 0.2e1 + t184 ^ 2 / 0.2e1, t174 ^ 2 / 0.2e1, -t174 * t173, t174 * t177, -t173 * t177, t177 ^ 2 / 0.2e1, t168 * t173 + t196 * t177, t168 * t174 - t202 * t177, t163 ^ 2 / 0.2e1, -t163 * t162, t163 * t175, -t162 * t175, t175 ^ 2 / 0.2e1, t161 * t162 + t195 * t175, t161 * t163 - t203 * t175, -t154 * t175 + t156 * t162, t154 * t163 - t155 * t162, t155 * t175 - t156 * t163, t155 ^ 2 / 0.2e1 + t156 ^ 2 / 0.2e1 + t154 ^ 2 / 0.2e1;];
T_reg  = t1;
