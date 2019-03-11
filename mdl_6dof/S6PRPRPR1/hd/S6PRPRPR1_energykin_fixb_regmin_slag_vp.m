% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% T_reg [1x21]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRPRPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:27:56
% EndTime: 2019-03-08 19:27:56
% DurationCPUTime: 0.08s
% Computational Cost: add. (198->40), mult. (467->89), div. (0->0), fcn. (341->12), ass. (0->40)
t202 = qJD(2) ^ 2;
t210 = t202 / 0.2e1;
t201 = cos(qJ(2));
t208 = qJD(1) * sin(pkin(6));
t183 = qJD(2) * pkin(2) + t201 * t208;
t192 = sin(pkin(11));
t195 = cos(pkin(11));
t198 = sin(qJ(2));
t204 = t198 * t208;
t176 = t192 * t183 + t195 * t204;
t174 = qJD(2) * pkin(8) + t176;
t189 = cos(pkin(6)) * qJD(1) + qJD(3);
t200 = cos(qJ(4));
t187 = t200 * t189;
t197 = sin(qJ(4));
t169 = qJD(4) * pkin(4) + t187 + (-qJ(5) * qJD(2) - t174) * t197;
t206 = qJD(2) * t200;
t209 = t200 * t174 + t197 * t189;
t170 = qJ(5) * t206 + t209;
t191 = sin(pkin(12));
t194 = cos(pkin(12));
t165 = t191 * t169 + t194 * t170;
t207 = qJD(2) * t197;
t205 = qJD(2) * qJD(4);
t203 = qJD(2) * t208;
t175 = t195 * t183 - t192 * t204;
t180 = -t191 * t207 + t194 * t206;
t164 = t194 * t169 - t191 * t170;
t171 = qJD(5) + (-pkin(4) * t200 - pkin(3)) * qJD(2) - t175;
t199 = cos(qJ(6));
t196 = sin(qJ(6));
t181 = (t191 * t200 + t194 * t197) * qJD(2);
t179 = qJD(6) - t180;
t178 = t196 * qJD(4) + t199 * t181;
t177 = -t199 * qJD(4) + t196 * t181;
t173 = -qJD(2) * pkin(3) - t175;
t166 = -t180 * pkin(5) - t181 * pkin(9) + t171;
t163 = qJD(4) * pkin(9) + t165;
t162 = -qJD(4) * pkin(5) - t164;
t1 = [qJD(1) ^ 2 / 0.2e1, t210, t201 * t203, -t198 * t203, t176 ^ 2 / 0.2e1 + t175 ^ 2 / 0.2e1 + t189 ^ 2 / 0.2e1, t197 ^ 2 * t210, t197 * t202 * t200, t197 * t205, t200 * t205, qJD(4) ^ 2 / 0.2e1 (-t197 * t174 + t187) * qJD(4) - t173 * t206, -t209 * qJD(4) + t173 * t207, -t164 * t181 + t165 * t180, t165 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1 + t171 ^ 2 / 0.2e1, t178 ^ 2 / 0.2e1, -t178 * t177, t178 * t179, -t177 * t179, t179 ^ 2 / 0.2e1 (-t196 * t163 + t199 * t166) * t179 + t162 * t177 -(t199 * t163 + t196 * t166) * t179 + t162 * t178;];
T_reg  = t1;
