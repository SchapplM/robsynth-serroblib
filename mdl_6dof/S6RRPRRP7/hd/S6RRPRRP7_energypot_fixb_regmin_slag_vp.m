% Calculate minimal parameter regressor of potential energy for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:20:05
% EndTime: 2019-03-09 12:20:05
% DurationCPUTime: 0.10s
% Computational Cost: add. (107->54), mult. (240->73), div. (0->0), fcn. (279->8), ass. (0->33)
t221 = sin(qJ(4));
t222 = sin(qJ(2));
t225 = cos(qJ(4));
t226 = cos(qJ(2));
t240 = t226 * t221 - t222 * t225;
t239 = g(3) * t240;
t227 = cos(qJ(1));
t237 = t222 * t227;
t223 = sin(qJ(1));
t236 = t223 * t226;
t234 = t226 * t227;
t233 = pkin(2) * t236 + (qJ(3) * t222 + pkin(1)) * t223;
t232 = t227 * pkin(1) + pkin(2) * t234 + t223 * pkin(7) + qJ(3) * t237;
t231 = t222 * pkin(2) - t226 * qJ(3) + pkin(6);
t230 = g(1) * t227 + g(2) * t223;
t206 = t222 * t221 + t226 * t225;
t202 = t240 * t223;
t204 = t221 * t234 - t225 * t237;
t229 = g(1) * t204 + g(2) * t202 + g(3) * t206;
t203 = t206 * t223;
t220 = sin(qJ(5));
t224 = cos(qJ(5));
t196 = t203 * t220 - t227 * t224;
t205 = t206 * t227;
t198 = t205 * t220 + t223 * t224;
t228 = g(1) * t198 + g(2) * t196 - t220 * t239;
t208 = g(1) * t223 - g(2) * t227;
t201 = -g(3) * t222 - t230 * t226;
t200 = -g(3) * t226 + t230 * t222;
t199 = t205 * t224 - t223 * t220;
t197 = t203 * t224 + t227 * t220;
t195 = -g(1) * t199 - g(2) * t197 + t224 * t239;
t1 = [0, -t230, t208, 0, 0, 0, 0, 0, t201, t200, t201, -t208, -t200, -g(1) * t232 - g(2) * (-t227 * pkin(7) + t233) - g(3) * t231, 0, 0, 0, 0, 0, -g(1) * t205 - g(2) * t203 + t239, t229, 0, 0, 0, 0, 0, t195, t228, t195, -t229, -t228, -g(1) * (pkin(3) * t234 + t205 * pkin(4) + t199 * pkin(5) - t223 * pkin(8) + t204 * pkin(9) + t198 * qJ(6) + t232) - g(2) * (pkin(3) * t236 + t203 * pkin(4) + t197 * pkin(5) + t202 * pkin(9) + t196 * qJ(6) + (-pkin(7) + pkin(8)) * t227 + t233) - g(3) * (t222 * pkin(3) + t206 * pkin(9) + t231) + (pkin(5) * t224 + qJ(6) * t220 + pkin(4)) * t239;];
U_reg  = t1;
