% Calculate minimal parameter regressor of potential energy for
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:49:16
% EndTime: 2019-03-08 21:49:16
% DurationCPUTime: 0.13s
% Computational Cost: add. (197->66), mult. (467->93), div. (0->0), fcn. (586->10), ass. (0->43)
t244 = pkin(4) + pkin(8);
t221 = sin(pkin(6));
t243 = pkin(7) * t221;
t225 = sin(qJ(3));
t242 = t221 * t225;
t226 = sin(qJ(2));
t241 = t221 * t226;
t228 = cos(qJ(3));
t240 = t221 * t228;
t229 = cos(qJ(2));
t239 = t221 * t229;
t223 = cos(pkin(6));
t238 = t223 * t226;
t237 = t223 * t229;
t209 = -t223 * t228 + t225 * t241;
t210 = t223 * t225 + t226 * t240;
t236 = pkin(2) * t241 + t210 * pkin(3) + t223 * pkin(7) + t209 * qJ(4) + qJ(1);
t220 = sin(pkin(10));
t222 = cos(pkin(10));
t208 = -t220 * t238 + t222 * t229;
t198 = t208 * t225 - t220 * t240;
t199 = t208 * t228 + t220 * t242;
t235 = t222 * pkin(1) + t208 * pkin(2) + t199 * pkin(3) + t198 * qJ(4) + t220 * t243;
t206 = t220 * t229 + t222 * t238;
t196 = t206 * t225 + t222 * t240;
t205 = t220 * t226 - t222 * t237;
t224 = sin(qJ(5));
t227 = cos(qJ(5));
t188 = -t196 * t227 + t205 * t224;
t207 = t220 * t237 + t222 * t226;
t190 = -t198 * t227 + t207 * t224;
t200 = t209 * t227 + t224 * t239;
t234 = g(1) * t190 + g(2) * t188 - g(3) * t200;
t233 = g(1) * t198 + g(2) * t196 + g(3) * t209;
t197 = t206 * t228 - t222 * t242;
t232 = g(1) * t199 + g(2) * t197 + g(3) * t210;
t231 = -g(1) * t207 - g(2) * t205 + g(3) * t239;
t230 = t220 * pkin(1) + t206 * pkin(2) + t197 * pkin(3) + t196 * qJ(4) - t222 * t243;
t201 = t209 * t224 - t227 * t239;
t191 = t198 * t224 + t207 * t227;
t189 = t196 * t224 + t205 * t227;
t186 = -g(1) * t191 - g(2) * t189 - g(3) * t201;
t1 = [-g(3) * qJ(1), 0, -g(1) * t208 - g(2) * t206 - g(3) * t241, -t231, 0, 0, 0, 0, 0, -t232, t233, t231, t232, -t233, -g(1) * (t207 * pkin(8) + t235) - g(2) * (t205 * pkin(8) + t230) - g(3) * (-pkin(8) * t239 + t236) 0, 0, 0, 0, 0, t186, t234, t186, -t232, -t234, -g(1) * (t191 * pkin(5) + t199 * pkin(9) + t190 * qJ(6) + t244 * t207 + t235) - g(2) * (t189 * pkin(5) + t197 * pkin(9) + t188 * qJ(6) + t244 * t205 + t230) - g(3) * (t201 * pkin(5) + t210 * pkin(9) - t200 * qJ(6) - t244 * t239 + t236);];
U_reg  = t1;
