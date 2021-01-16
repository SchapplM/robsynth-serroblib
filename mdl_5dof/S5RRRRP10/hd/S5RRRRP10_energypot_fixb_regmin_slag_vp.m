% Calculate minimal parameter regressor of potential energy for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:35:13
% EndTime: 2021-01-16 00:35:14
% DurationCPUTime: 0.20s
% Computational Cost: add. (131->54), mult. (259->93), div. (0->0), fcn. (315->10), ass. (0->41)
t215 = cos(pkin(5));
t223 = cos(qJ(2));
t224 = cos(qJ(1));
t229 = t224 * t223;
t219 = sin(qJ(2));
t220 = sin(qJ(1));
t234 = t220 * t219;
t204 = t215 * t234 - t229;
t230 = t224 * t219;
t233 = t220 * t223;
t205 = -t215 * t230 - t233;
t218 = sin(qJ(3));
t222 = cos(qJ(3));
t214 = sin(pkin(5));
t236 = t218 * t214;
t239 = t214 * t222;
t242 = g(1) * (t204 * t218 + t220 * t239) - g(2) * (-t205 * t218 + t224 * t239) - g(3) * (-t215 * t222 + t219 * t236);
t240 = t214 * t219;
t238 = t214 * t223;
t217 = sin(qJ(4));
t237 = t217 * t223;
t235 = t219 * t222;
t221 = cos(qJ(4));
t232 = t221 * t223;
t231 = t222 * t223;
t213 = t221 * pkin(4) + pkin(3);
t216 = qJ(5) + pkin(9);
t228 = t213 * t218 - t216 * t222 + pkin(7);
t203 = t215 * t235 - t236;
t227 = t203 * t217 + t215 * t232;
t212 = t217 * pkin(4) + pkin(8);
t207 = t217 * t231 - t221 * t219;
t206 = t215 * t233 + t230;
t202 = t214 * t235 + t215 * t218;
t200 = t213 * t222 + t216 * t218 + pkin(2);
t199 = -t203 * t220 + t222 * t229;
t198 = t200 * t223 + t212 * t219 + pkin(1);
t197 = t214 * t228 + (-t200 * t219 + t212 * t223) * t215;
t196 = -g(1) * (t199 * t221 + t206 * t217) - g(2) * ((t203 * t221 - t215 * t237) * t224 + t220 * (t217 * t219 + t221 * t231)) - g(3) * (t202 * t221 - t214 * t237);
t195 = -g(1) * (-t224 * t207 + t227 * t220) - g(2) * (-t220 * t207 - t227 * t224) - g(3) * (-t202 * t217 - t214 * t232);
t1 = [0, -g(1) * t224 - g(2) * t220, g(1) * t220 - g(2) * t224, 0, 0, 0, 0, 0, g(1) * t204 + g(2) * t205 - g(3) * t240, g(1) * t206 - g(2) * (t215 * t229 - t234) - g(3) * t238, 0, 0, 0, 0, 0, -g(1) * t199 - g(2) * (-t205 * t222 - t224 * t236) - g(3) * t202, -t242, 0, 0, 0, 0, 0, t196, t195, t196, t195, t242, -g(1) * (t197 * t220 + t198 * t224) - g(2) * (-t197 * t224 + t198 * t220) + (-t200 * t240 + t212 * t238 - t228 * t215 - pkin(6)) * g(3);];
U_reg = t1;
