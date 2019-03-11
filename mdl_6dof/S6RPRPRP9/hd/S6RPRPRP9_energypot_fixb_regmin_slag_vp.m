% Calculate minimal parameter regressor of potential energy for
% S6RPRPRP9
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
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRP9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:29:27
% EndTime: 2019-03-09 03:29:27
% DurationCPUTime: 0.12s
% Computational Cost: add. (114->58), mult. (160->73), div. (0->0), fcn. (167->8), ass. (0->39)
t207 = pkin(2) + pkin(6);
t190 = cos(qJ(1));
t206 = g(2) * t190;
t189 = cos(qJ(3));
t205 = g(3) * t189;
t185 = cos(pkin(9));
t175 = t185 * pkin(4) + pkin(3);
t187 = sin(qJ(3));
t204 = t175 * t187;
t183 = pkin(9) + qJ(5);
t176 = sin(t183);
t188 = sin(qJ(1));
t203 = t188 * t176;
t177 = cos(t183);
t202 = t188 * t177;
t184 = sin(pkin(9));
t201 = t188 * t184;
t200 = t188 * t185;
t199 = t190 * t176;
t198 = t190 * t177;
t197 = t190 * t184;
t196 = t190 * t185;
t195 = t190 * pkin(1) + t188 * qJ(2);
t194 = t190 * pkin(7) + t195;
t180 = t188 * pkin(1);
t193 = -t190 * qJ(2) + t180;
t172 = g(1) * t188 - t206;
t192 = pkin(3) * t187 - qJ(4) * t189;
t167 = t187 * t203 - t198;
t169 = t187 * t199 + t202;
t191 = g(1) * t167 - g(2) * t169 + t176 * t205;
t186 = -pkin(8) - qJ(4);
t179 = t188 * pkin(7);
t173 = g(1) * t190 + g(2) * t188;
t171 = -g(3) * t187 + t172 * t189;
t170 = -t187 * t198 + t203;
t168 = t187 * t202 + t199;
t166 = -g(1) * t168 - g(2) * t170 - t177 * t205;
t1 = [0, -t173, t172, t173, -t172, -g(3) * pkin(6) - g(1) * t195 - g(2) * t193, 0, 0, 0, 0, 0, -t172 * t187 - t205, -t171, -g(1) * (t187 * t200 + t197) - g(2) * (-t187 * t196 + t201) - t185 * t205, -g(1) * (-t187 * t201 + t196) - g(2) * (t187 * t197 + t200) + t184 * t205, t171, -g(1) * (t192 * t188 + t194) - g(2) * (t179 + t180) - g(3) * (t189 * pkin(3) + t187 * qJ(4) + t207) - (-qJ(2) - t192) * t206, 0, 0, 0, 0, 0, t166, t191, t166, t171, -t191, -g(1) * (pkin(4) * t197 + t168 * pkin(5) + t167 * qJ(6) + t188 * t204 + t194) - g(2) * (pkin(4) * t201 + t170 * pkin(5) - t169 * qJ(6) - t190 * t204 + t179 + t193) - g(3) * (-t187 * t186 + t207) + (-g(3) * (pkin(5) * t177 + qJ(6) * t176 + t175) - t172 * t186) * t189;];
U_reg  = t1;
