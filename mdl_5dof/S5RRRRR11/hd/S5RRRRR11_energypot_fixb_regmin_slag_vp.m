% Calculate minimal parameter regressor of potential energy for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x31]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR11_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR11_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:43:41
% EndTime: 2019-12-31 22:43:41
% DurationCPUTime: 0.09s
% Computational Cost: add. (84->41), mult. (184->77), div. (0->0), fcn. (238->12), ass. (0->29)
t215 = sin(pkin(5));
t219 = sin(qJ(2));
t232 = t215 * t219;
t222 = cos(qJ(3));
t231 = t215 * t222;
t223 = cos(qJ(2));
t230 = t215 * t223;
t224 = cos(qJ(1));
t229 = t215 * t224;
t220 = sin(qJ(1));
t228 = t220 * t219;
t227 = t220 * t223;
t226 = t224 * t219;
t225 = t224 * t223;
t221 = cos(qJ(4));
t218 = sin(qJ(3));
t217 = sin(qJ(4));
t216 = cos(pkin(5));
t214 = qJ(4) + qJ(5);
t213 = cos(t214);
t212 = sin(t214);
t211 = -t216 * t228 + t225;
t210 = t216 * t227 + t226;
t209 = t216 * t226 + t227;
t208 = -t216 * t225 + t228;
t207 = t216 * t218 + t219 * t231;
t206 = t220 * t215 * t218 + t211 * t222;
t205 = t209 * t222 - t218 * t229;
t1 = [0, -g(1) * t224 - g(2) * t220, g(1) * t220 - g(2) * t224, 0, 0, 0, 0, 0, -g(1) * t211 - g(2) * t209 - g(3) * t232, g(1) * t210 + g(2) * t208 - g(3) * t230, 0, 0, 0, 0, 0, -g(1) * t206 - g(2) * t205 - g(3) * t207, -g(1) * (-t211 * t218 + t220 * t231) - g(2) * (-t209 * t218 - t222 * t229) - g(3) * (t216 * t222 - t218 * t232), 0, 0, 0, 0, 0, -g(1) * (t206 * t221 + t210 * t217) - g(2) * (t205 * t221 + t208 * t217) - g(3) * (t207 * t221 - t217 * t230), -g(1) * (-t206 * t217 + t210 * t221) - g(2) * (-t205 * t217 + t208 * t221) - g(3) * (-t207 * t217 - t221 * t230), 0, 0, 0, 0, 0, -g(1) * (t206 * t213 + t210 * t212) - g(2) * (t205 * t213 + t208 * t212) - g(3) * (t207 * t213 - t212 * t230), -g(1) * (-t206 * t212 + t210 * t213) - g(2) * (-t205 * t212 + t208 * t213) - g(3) * (-t207 * t212 - t213 * t230);];
U_reg = t1;
