% Calculate minimal parameter regressor of potential energy for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:35:32
% EndTime: 2019-03-09 08:35:32
% DurationCPUTime: 0.11s
% Computational Cost: add. (80->51), mult. (145->59), div. (0->0), fcn. (140->6), ass. (0->31)
t182 = sin(qJ(2));
t205 = qJ(3) * t182 + pkin(1);
t186 = cos(qJ(1));
t183 = sin(qJ(1));
t185 = cos(qJ(2));
t197 = t183 * t185;
t204 = pkin(3) * t197 + t186 * qJ(4);
t181 = sin(qJ(5));
t203 = pkin(5) * t181;
t202 = g(3) * t185;
t201 = t182 * pkin(2) + pkin(6);
t199 = t183 * t181;
t184 = cos(qJ(5));
t198 = t183 * t184;
t196 = t185 * t186;
t195 = t186 * t181;
t194 = t186 * t184;
t193 = pkin(2) * t197 + t205 * t183;
t192 = pkin(2) * t196 + t183 * pkin(7) + t205 * t186;
t191 = -t185 * qJ(3) + t201;
t190 = g(1) * t186 + g(2) * t183;
t189 = pkin(3) * t196 + t192;
t172 = t184 * pkin(5) + pkin(4);
t180 = -qJ(6) - pkin(8);
t188 = t172 * t182 - t180 * t185;
t187 = -t186 * pkin(7) + t193;
t174 = t182 * pkin(3);
t165 = g(1) * t183 - g(2) * t186;
t164 = g(3) * t182 + t190 * t185;
t163 = t190 * t182 - t202;
t1 = [0, -t190, t165, 0, 0, 0, 0, 0, -t164, t163, -t164, -t165, -t163, -g(1) * t192 - g(2) * t187 - g(3) * t191, -t163, t164, t165, -g(1) * (-t183 * qJ(4) + t189) - g(2) * (t187 + t204) - g(3) * (t174 + t191) 0, 0, 0, 0, 0, -g(1) * (t182 * t194 - t199) - g(2) * (t182 * t198 + t195) + t184 * t202, -g(1) * (-t182 * t195 - t198) - g(2) * (-t182 * t199 + t194) - t181 * t202, -t164, -g(1) * t189 - g(2) * (t193 + t204) - g(3) * (-t182 * t180 + t174 + (-qJ(3) - t172) * t185 + t201) + (-g(1) * t188 - g(2) * (-pkin(7) + t203)) * t186 + (-g(1) * (-qJ(4) - t203) - g(2) * t188) * t183;];
U_reg  = t1;
