% Calculate minimal parameter regressor of potential energy for
% S5RRRPR10
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
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:29:50
% EndTime: 2019-12-31 21:29:50
% DurationCPUTime: 0.12s
% Computational Cost: add. (88->50), mult. (173->83), div. (0->0), fcn. (210->12), ass. (0->35)
t226 = sin(qJ(1));
t230 = cos(qJ(1));
t245 = -g(1) * t226 + g(2) * t230;
t220 = sin(pkin(5));
t225 = sin(qJ(2));
t242 = t220 * t225;
t241 = t220 * t226;
t228 = cos(qJ(3));
t240 = t220 * t228;
t229 = cos(qJ(2));
t239 = t220 * t229;
t238 = t220 * t230;
t221 = cos(pkin(5));
t224 = sin(qJ(3));
t237 = t221 * t224;
t236 = t226 * t225;
t235 = t226 * t229;
t234 = t230 * t225;
t233 = t230 * t229;
t211 = -t221 * t233 + t236;
t213 = t221 * t235 + t234;
t231 = -g(1) * t213 - g(2) * t211 + g(3) * t239;
t227 = cos(qJ(5));
t223 = sin(qJ(5));
t222 = -qJ(4) - pkin(8);
t219 = qJ(3) + pkin(10);
t218 = cos(t219);
t217 = sin(t219);
t216 = t228 * pkin(3) + pkin(2);
t214 = -t221 * t236 + t233;
t212 = t221 * t234 + t235;
t210 = t221 * t217 + t218 * t242;
t209 = t214 * t218 + t217 * t241;
t208 = t212 * t218 - t217 * t238;
t1 = [0, -g(1) * t230 - g(2) * t226, -t245, 0, 0, 0, 0, 0, -g(1) * t214 - g(2) * t212 - g(3) * t242, -t231, 0, 0, 0, 0, 0, -g(1) * (t214 * t228 + t224 * t241) - g(2) * (t212 * t228 - t224 * t238) - g(3) * (t225 * t240 + t237), -g(1) * (-t214 * t224 + t226 * t240) - g(2) * (-t212 * t224 - t228 * t238) - g(3) * (t221 * t228 - t224 * t242), t231, -g(1) * (t230 * pkin(1) - t213 * t222 + t214 * t216) - g(2) * (t226 * pkin(1) - t211 * t222 + t212 * t216) - g(3) * (pkin(3) * t237 + t221 * pkin(7) + pkin(6)) + (-g(3) * (t216 * t225 + t222 * t229) + t245 * (pkin(3) * t224 + pkin(7))) * t220, 0, 0, 0, 0, 0, -g(1) * (t209 * t227 + t213 * t223) - g(2) * (t208 * t227 + t211 * t223) - g(3) * (t210 * t227 - t223 * t239), -g(1) * (-t209 * t223 + t213 * t227) - g(2) * (-t208 * t223 + t211 * t227) - g(3) * (-t210 * t223 - t227 * t239);];
U_reg = t1;
