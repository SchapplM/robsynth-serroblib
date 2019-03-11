% Calculate minimal parameter regressor of potential energy for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% U_reg [1x31]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:32:32
% EndTime: 2019-03-09 06:32:32
% DurationCPUTime: 0.08s
% Computational Cost: add. (99->49), mult. (136->64), div. (0->0), fcn. (146->8), ass. (0->34)
t185 = cos(qJ(3));
t200 = g(3) * t185;
t180 = qJ(4) + qJ(5);
t175 = sin(t180);
t183 = sin(qJ(1));
t199 = t183 * t175;
t176 = cos(t180);
t198 = t183 * t176;
t181 = sin(qJ(4));
t197 = t183 * t181;
t184 = cos(qJ(4));
t196 = t183 * t184;
t186 = cos(qJ(1));
t195 = t186 * t175;
t194 = t186 * t176;
t193 = t186 * t181;
t192 = t186 * t184;
t191 = t186 * pkin(1) + t183 * qJ(2);
t190 = pkin(4) * t181 + pkin(7);
t171 = g(1) * t183 - g(2) * t186;
t174 = t184 * pkin(4) + pkin(3);
t182 = sin(qJ(3));
t187 = -pkin(9) - pkin(8);
t189 = t174 * t182 + t185 * t187;
t166 = t182 * t199 - t194;
t168 = t182 * t195 + t198;
t188 = g(1) * t166 - g(2) * t168 + t175 * t200;
t178 = t183 * pkin(1);
t172 = g(1) * t186 + g(2) * t183;
t170 = -g(3) * t182 + t171 * t185;
t169 = -t182 * t194 + t199;
t167 = t182 * t198 + t195;
t165 = -g(1) * t167 - g(2) * t169 - t176 * t200;
t1 = [0, -t172, t171, t172, -t171, -g(1) * t191 - g(2) * (-t186 * qJ(2) + t178) - g(3) * pkin(6), 0, 0, 0, 0, 0, -t171 * t182 - t200, -t170, 0, 0, 0, 0, 0, -g(1) * (t182 * t196 + t193) - g(2) * (-t182 * t192 + t197) - t184 * t200, -g(1) * (-t182 * t197 + t192) - g(2) * (t182 * t193 + t196) + t181 * t200, 0, 0, 0, 0, 0, t165, t188, t165, t170, -t188, -g(1) * (t167 * pkin(5) + t166 * qJ(6) + t191) - g(2) * (t169 * pkin(5) - t168 * qJ(6) + t178) - g(3) * (-t182 * t187 + pkin(2) + pkin(6)) - (pkin(5) * t176 + qJ(6) * t175 + t174) * t200 + (-g(1) * t189 - g(2) * t190) * t183 + (-g(1) * t190 - g(2) * (-qJ(2) - t189)) * t186;];
U_reg  = t1;
