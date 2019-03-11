% Calculate minimal parameter regressor of potential energy for
% S6RRPRRP12
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
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP12_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP12_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:53:31
% EndTime: 2019-03-09 12:53:31
% DurationCPUTime: 0.12s
% Computational Cost: add. (110->53), mult. (162->66), div. (0->0), fcn. (172->8), ass. (0->37)
t214 = sin(qJ(2));
t217 = cos(qJ(2));
t238 = pkin(2) * t217 + qJ(3) * t214 + pkin(1);
t236 = g(3) * t217;
t235 = t214 * pkin(2) + pkin(6);
t212 = qJ(4) + qJ(5);
t206 = sin(t212);
t215 = sin(qJ(1));
t233 = t215 * t206;
t207 = cos(t212);
t232 = t215 * t207;
t213 = sin(qJ(4));
t231 = t215 * t213;
t216 = cos(qJ(4));
t230 = t215 * t216;
t218 = cos(qJ(1));
t229 = t218 * t206;
t228 = t218 * t207;
t227 = t218 * t213;
t226 = t218 * t216;
t225 = t214 * t231;
t224 = t214 * t227;
t223 = t215 * pkin(7) + t238 * t218;
t222 = g(1) * t218 + g(2) * t215;
t221 = -t218 * pkin(7) + t238 * t215;
t193 = -t214 * t228 + t233;
t195 = t214 * t232 + t229;
t220 = g(1) * t193 - g(2) * t195 + t207 * t236;
t219 = -pkin(9) - pkin(8);
t205 = t216 * pkin(4) + pkin(3);
t200 = g(1) * t215 - g(2) * t218;
t198 = g(3) * t214 + t222 * t217;
t197 = t222 * t214 - t236;
t196 = t214 * t233 - t228;
t194 = t214 * t229 + t232;
t192 = -g(1) * t194 - g(2) * t196 + t206 * t236;
t1 = [0, -t222, t200, 0, 0, 0, 0, 0, -t198, t197, -t200, t198, -t197, -g(1) * t223 - g(2) * t221 - g(3) * (-t217 * qJ(3) + t235) 0, 0, 0, 0, 0, -g(1) * (t224 + t230) - g(2) * (t225 - t226) + t213 * t236, -g(1) * (t214 * t226 - t231) - g(2) * (t214 * t230 + t227) + t216 * t236, 0, 0, 0, 0, 0, t192, t220, t192, -t198, -t220, -g(1) * (pkin(4) * t224 + t194 * pkin(5) + t193 * qJ(6) + t215 * t205 + t223) - g(2) * (pkin(4) * t225 + t196 * pkin(5) - t195 * qJ(6) - t218 * t205 + t221) - g(3) * (-t214 * t219 + t235) + (-g(3) * (-pkin(4) * t213 - pkin(5) * t206 + qJ(6) * t207 - qJ(3)) + t222 * t219) * t217;];
U_reg  = t1;
