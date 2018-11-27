% Calculate minimal parameter regressor of potential energy for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% 
% Output:
% U_reg [1x45]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S7RRRRRRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_energypot_fixb_regmin_slag_vp: qJ has to be [7x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S7RRRRRRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_energypot_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 20:34:13
% EndTime: 2018-11-26 20:34:13
% DurationCPUTime: 0.13s
% Computational Cost: add. (158->53), mult. (408->99), div. (0->0), fcn. (540->14), ass. (0->42)
t217 = sin(qJ(3));
t218 = sin(qJ(2));
t234 = t217 * t218;
t219 = sin(qJ(1));
t233 = t218 * t219;
t224 = cos(qJ(3));
t232 = t218 * t224;
t226 = cos(qJ(1));
t231 = t218 * t226;
t225 = cos(qJ(2));
t230 = t219 * t225;
t229 = t226 * t217;
t228 = t226 * t224;
t227 = g(1) * t226 + g(2) * t219;
t223 = cos(qJ(4));
t222 = cos(qJ(5));
t221 = cos(qJ(6));
t220 = cos(qJ(7));
t216 = sin(qJ(4));
t215 = sin(qJ(5));
t214 = sin(qJ(6));
t213 = sin(qJ(7));
t212 = -t219 * t217 + t225 * t228;
t211 = -t219 * t224 - t225 * t229;
t210 = t224 * t230 + t229;
t209 = -t217 * t230 + t228;
t208 = -t225 * t216 + t223 * t232;
t207 = t216 * t232 + t225 * t223;
t206 = t212 * t223 + t216 * t231;
t205 = t212 * t216 - t223 * t231;
t204 = t210 * t223 + t216 * t233;
t203 = t210 * t216 - t223 * t233;
t202 = t208 * t222 - t215 * t234;
t201 = t208 * t215 + t222 * t234;
t200 = t206 * t222 + t211 * t215;
t199 = t206 * t215 - t211 * t222;
t198 = t204 * t222 + t209 * t215;
t197 = t204 * t215 - t209 * t222;
t196 = t202 * t221 + t207 * t214;
t195 = t200 * t221 + t205 * t214;
t194 = t198 * t221 + t203 * t214;
t1 = [0, -t227, g(1) * t219 - g(2) * t226, 0, 0, 0, 0, 0, -g(3) * t218 - t227 * t225, -g(3) * t225 + t227 * t218, 0, 0, 0, 0, 0, -g(1) * t212 - g(2) * t210 - g(3) * t232, -g(1) * t211 - g(2) * t209 + g(3) * t234, 0, 0, 0, 0, 0, -g(1) * t206 - g(2) * t204 - g(3) * t208, g(1) * t205 + g(2) * t203 + g(3) * t207, 0, 0, 0, 0, 0, -g(1) * t200 - g(2) * t198 - g(3) * t202, g(1) * t199 + g(2) * t197 + g(3) * t201, 0, 0, 0, 0, 0, -g(1) * t195 - g(2) * t194 - g(3) * t196, -g(1) * (-t200 * t214 + t205 * t221) - g(2) * (-t198 * t214 + t203 * t221) - g(3) * (-t202 * t214 + t207 * t221) 0, 0, 0, 0, 0, -g(1) * (t195 * t220 - t199 * t213) - g(2) * (t194 * t220 - t197 * t213) - g(3) * (t196 * t220 - t201 * t213) -g(1) * (-t195 * t213 - t199 * t220) - g(2) * (-t194 * t213 - t197 * t220) - g(3) * (-t196 * t213 - t201 * t220);];
U_reg  = t1;
