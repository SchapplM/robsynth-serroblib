% Calculate minimal parameter regressor of potential energy for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRP8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:26:03
% EndTime: 2019-03-09 03:26:04
% DurationCPUTime: 0.10s
% Computational Cost: add. (98->48), mult. (130->60), div. (0->0), fcn. (133->8), ass. (0->33)
t184 = sin(qJ(3));
t202 = pkin(3) * t184;
t181 = qJ(3) + pkin(9);
t175 = sin(t181);
t201 = pkin(4) * t175;
t185 = sin(qJ(1));
t200 = g(1) * t185;
t176 = cos(t181);
t199 = g(3) * t176;
t183 = sin(qJ(5));
t198 = t185 * t183;
t186 = cos(qJ(5));
t197 = t185 * t186;
t188 = cos(qJ(1));
t196 = t188 * t183;
t195 = t188 * t186;
t194 = pkin(1) * t188 + qJ(2) * t185;
t187 = cos(qJ(3));
t193 = pkin(3) * t187 + pkin(2) + pkin(6);
t192 = t185 * t202 + t194;
t191 = -qJ(2) - t202;
t178 = t185 * pkin(1);
t182 = -qJ(4) - pkin(7);
t190 = -t185 * t182 + t178;
t171 = -g(2) * t188 + t200;
t167 = t175 * t198 - t195;
t169 = t175 * t196 + t197;
t189 = g(1) * t167 - g(2) * t169 + t183 * t199;
t172 = g(1) * t188 + g(2) * t185;
t170 = -t175 * t195 + t198;
t168 = t175 * t197 + t196;
t166 = -g(1) * t168 - g(2) * t170 - t186 * t199;
t1 = [0, -t172, t171, t172, -t171, -g(1) * t194 - g(2) * (-qJ(2) * t188 + t178) - g(3) * pkin(6), 0, 0, 0, 0, 0, -g(3) * t187 - t171 * t184, g(3) * t184 - t171 * t187, -t172, -g(1) * (-t188 * t182 + t192) - g(2) * (t188 * t191 + t190) - g(3) * t193, 0, 0, 0, 0, 0, t166, t189, t166, -g(3) * t175 + t171 * t176, -t189, -g(1) * (t168 * pkin(5) + t167 * qJ(6) + t185 * t201 + t192) - g(2) * (t170 * pkin(5) - t169 * qJ(6) + t190) - g(3) * (t175 * pkin(8) + t193) + (pkin(8) * t200 - g(3) * (pkin(5) * t186 + qJ(6) * t183 + pkin(4))) * t176 + (g(1) * t182 - g(2) * (pkin(8) * t176 + t191 - t201)) * t188;];
U_reg  = t1;
