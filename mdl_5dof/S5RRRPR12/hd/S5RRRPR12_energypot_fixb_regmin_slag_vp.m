% Calculate minimal parameter regressor of potential energy for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR12_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR12_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:40:34
% EndTime: 2019-12-31 21:40:34
% DurationCPUTime: 0.11s
% Computational Cost: add. (122->59), mult. (268->98), div. (0->0), fcn. (338->12), ass. (0->34)
t220 = sin(pkin(5));
t224 = sin(qJ(2));
t238 = t220 * t224;
t225 = sin(qJ(1));
t237 = t220 * t225;
t226 = cos(qJ(3));
t236 = t220 * t226;
t227 = cos(qJ(2));
t235 = t220 * t227;
t228 = cos(qJ(1));
t234 = t220 * t228;
t233 = t225 * t224;
t232 = t225 * t227;
t231 = t228 * t224;
t230 = t228 * t227;
t222 = cos(pkin(5));
t210 = t222 * t231 + t232;
t223 = sin(qJ(3));
t203 = t210 * t223 + t226 * t234;
t212 = -t222 * t233 + t230;
t205 = t212 * t223 - t225 * t236;
t207 = -t222 * t226 + t223 * t238;
t229 = g(1) * t205 + g(2) * t203 + g(3) * t207;
t221 = cos(pkin(10));
t219 = sin(pkin(10));
t218 = pkin(10) + qJ(5);
t217 = cos(t218);
t216 = sin(t218);
t211 = t222 * t232 + t231;
t209 = -t222 * t230 + t233;
t208 = t222 * t223 + t224 * t236;
t206 = t212 * t226 + t223 * t237;
t204 = t210 * t226 - t223 * t234;
t1 = [0, -g(1) * t228 - g(2) * t225, g(1) * t225 - g(2) * t228, 0, 0, 0, 0, 0, -g(1) * t212 - g(2) * t210 - g(3) * t238, g(1) * t211 + g(2) * t209 - g(3) * t235, 0, 0, 0, 0, 0, -g(1) * t206 - g(2) * t204 - g(3) * t208, t229, -g(1) * (t206 * t221 + t211 * t219) - g(2) * (t204 * t221 + t209 * t219) - g(3) * (t208 * t221 - t219 * t235), -g(1) * (-t206 * t219 + t211 * t221) - g(2) * (-t204 * t219 + t209 * t221) - g(3) * (-t208 * t219 - t221 * t235), -t229, -g(1) * (t228 * pkin(1) + t212 * pkin(2) + t206 * pkin(3) + pkin(7) * t237 + t211 * pkin(8) + t205 * qJ(4)) - g(2) * (t225 * pkin(1) + t210 * pkin(2) + t204 * pkin(3) - pkin(7) * t234 + t209 * pkin(8) + t203 * qJ(4)) - g(3) * (t208 * pkin(3) + t222 * pkin(7) + t207 * qJ(4) + pkin(6) + (pkin(2) * t224 - pkin(8) * t227) * t220), 0, 0, 0, 0, 0, -g(1) * (t206 * t217 + t211 * t216) - g(2) * (t204 * t217 + t209 * t216) - g(3) * (t208 * t217 - t216 * t235), -g(1) * (-t206 * t216 + t211 * t217) - g(2) * (-t204 * t216 + t209 * t217) - g(3) * (-t208 * t216 - t217 * t235);];
U_reg = t1;
