% Calculate minimal parameter regressor of potential energy for
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% U_reg [1x31]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:09:19
% EndTime: 2019-03-09 04:09:19
% DurationCPUTime: 0.10s
% Computational Cost: add. (81->48), mult. (109->69), div. (0->0), fcn. (115->10), ass. (0->32)
t191 = sin(qJ(3));
t193 = cos(qJ(3));
t211 = pkin(3) * t191 - qJ(4) * t193;
t209 = g(3) * t193;
t188 = pkin(10) + qJ(5);
t184 = qJ(6) + t188;
t180 = sin(t184);
t192 = sin(qJ(1));
t207 = t192 * t180;
t181 = cos(t184);
t206 = t192 * t181;
t182 = sin(t188);
t205 = t192 * t182;
t183 = cos(t188);
t204 = t192 * t183;
t189 = sin(pkin(10));
t203 = t192 * t189;
t190 = cos(pkin(10));
t202 = t192 * t190;
t194 = cos(qJ(1));
t201 = t194 * t180;
t200 = t194 * t181;
t199 = t194 * t182;
t198 = t194 * t183;
t197 = t194 * t189;
t196 = t194 * t190;
t195 = t194 * pkin(1) + t192 * qJ(2);
t178 = g(1) * t192 - g(2) * t194;
t186 = t192 * pkin(1);
t179 = g(1) * t194 + g(2) * t192;
t177 = -g(3) * t191 + t178 * t193;
t1 = [0, -t179, t178, t179, -t178, -g(1) * t195 - g(2) * (-t194 * qJ(2) + t186) - g(3) * pkin(6), 0, 0, 0, 0, 0, -t178 * t191 - t209, -t177, -g(1) * (t191 * t202 + t197) - g(2) * (-t191 * t196 + t203) - t190 * t209, -g(1) * (-t191 * t203 + t196) - g(2) * (t191 * t197 + t202) + t189 * t209, t177, -g(1) * (t211 * t192 + t195) - g(2) * (t192 * pkin(7) + t186) - g(3) * (t193 * pkin(3) + t191 * qJ(4) + pkin(2) + pkin(6)) + (-g(1) * pkin(7) - g(2) * (-qJ(2) - t211)) * t194, 0, 0, 0, 0, 0, -g(1) * (t191 * t204 + t199) - g(2) * (-t191 * t198 + t205) - t183 * t209, -g(1) * (-t191 * t205 + t198) - g(2) * (t191 * t199 + t204) + t182 * t209, 0, 0, 0, 0, 0, -g(1) * (t191 * t206 + t201) - g(2) * (-t191 * t200 + t207) - t181 * t209, -g(1) * (-t191 * t207 + t200) - g(2) * (t191 * t201 + t206) + t180 * t209;];
U_reg  = t1;
