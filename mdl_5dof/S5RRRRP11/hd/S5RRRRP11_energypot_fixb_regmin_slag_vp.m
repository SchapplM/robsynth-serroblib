% Calculate minimal parameter regressor of potential energy for
% S5RRRRP11
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
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP11_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP11_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:18:31
% EndTime: 2019-12-31 22:18:31
% DurationCPUTime: 0.13s
% Computational Cost: add. (136->55), mult. (330->86), div. (0->0), fcn. (420->10), ass. (0->39)
t215 = sin(pkin(5));
t219 = sin(qJ(2));
t235 = t215 * t219;
t220 = sin(qJ(1));
t234 = t215 * t220;
t222 = cos(qJ(3));
t233 = t215 * t222;
t223 = cos(qJ(2));
t232 = t215 * t223;
t224 = cos(qJ(1));
t231 = t215 * t224;
t230 = t220 * t219;
t229 = t220 * t223;
t228 = t224 * t219;
t227 = t224 * t223;
t216 = cos(pkin(5));
t209 = t216 * t228 + t229;
t218 = sin(qJ(3));
t201 = t209 * t222 - t218 * t231;
t208 = -t216 * t227 + t230;
t217 = sin(qJ(4));
t221 = cos(qJ(4));
t194 = t201 * t217 - t208 * t221;
t211 = -t216 * t230 + t227;
t203 = t211 * t222 + t218 * t234;
t210 = t216 * t229 + t228;
t196 = t203 * t217 - t210 * t221;
t207 = t216 * t218 + t219 * t233;
t198 = t207 * t217 + t221 * t232;
t226 = g(1) * t196 + g(2) * t194 + g(3) * t198;
t200 = t209 * t218 + t222 * t231;
t202 = t211 * t218 - t220 * t233;
t206 = -t216 * t222 + t218 * t235;
t225 = g(1) * t202 + g(2) * t200 + g(3) * t206;
t199 = t207 * t221 - t217 * t232;
t197 = t203 * t221 + t210 * t217;
t195 = t201 * t221 + t208 * t217;
t193 = -g(1) * t197 - g(2) * t195 - g(3) * t199;
t1 = [0, -g(1) * t224 - g(2) * t220, g(1) * t220 - g(2) * t224, 0, 0, 0, 0, 0, -g(1) * t211 - g(2) * t209 - g(3) * t235, g(1) * t210 + g(2) * t208 - g(3) * t232, 0, 0, 0, 0, 0, -g(1) * t203 - g(2) * t201 - g(3) * t207, t225, 0, 0, 0, 0, 0, t193, t226, t193, -t225, -t226, -g(1) * (t224 * pkin(1) + t211 * pkin(2) + t203 * pkin(3) + t197 * pkin(4) + pkin(7) * t234 + t210 * pkin(8) + t202 * pkin(9) + t196 * qJ(5)) - g(2) * (t220 * pkin(1) + t209 * pkin(2) + t201 * pkin(3) + t195 * pkin(4) - pkin(7) * t231 + t208 * pkin(8) + t200 * pkin(9) + t194 * qJ(5)) - g(3) * (t207 * pkin(3) + t199 * pkin(4) + t216 * pkin(7) + t206 * pkin(9) + t198 * qJ(5) + pkin(6) + (pkin(2) * t219 - pkin(8) * t223) * t215);];
U_reg = t1;
