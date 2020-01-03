% Calculate minimal parameter regressor of potential energy for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:26:38
% EndTime: 2019-12-31 20:26:39
% DurationCPUTime: 0.18s
% Computational Cost: add. (93->46), mult. (226->89), div. (0->0), fcn. (287->12), ass. (0->37)
t245 = pkin(7) + qJ(3);
t224 = sin(pkin(5));
t229 = sin(qJ(2));
t244 = t224 * t229;
t230 = sin(qJ(1));
t243 = t224 * t230;
t234 = cos(qJ(1));
t242 = t224 * t234;
t241 = t230 * t229;
t233 = cos(qJ(2));
t240 = t230 * t233;
t239 = t234 * t229;
t238 = t234 * t233;
t237 = g(1) * t230 - g(2) * t234;
t223 = sin(pkin(10));
t225 = cos(pkin(10));
t236 = t233 * t223 + t229 * t225;
t235 = t229 * t223 - t233 * t225;
t232 = cos(qJ(4));
t231 = cos(qJ(5));
t228 = sin(qJ(4));
t227 = sin(qJ(5));
t226 = cos(pkin(5));
t222 = t233 * pkin(2) + pkin(1);
t219 = t226 * t229 * pkin(2) - t245 * t224;
t218 = t236 * t226;
t217 = t235 * t226;
t216 = t236 * t224;
t215 = t235 * t224;
t214 = t216 * t232 + t226 * t228;
t213 = -t230 * t218 - t234 * t235;
t212 = -t230 * t217 + t234 * t236;
t211 = t234 * t218 - t230 * t235;
t210 = t234 * t217 + t230 * t236;
t209 = t213 * t232 + t228 * t243;
t208 = t211 * t232 - t228 * t242;
t1 = [0, -g(1) * t234 - g(2) * t230, t237, 0, 0, 0, 0, 0, -g(1) * (-t226 * t241 + t238) - g(2) * (t226 * t239 + t240) - g(3) * t244, -g(1) * (-t226 * t240 - t239) - g(2) * (t226 * t238 - t241) - g(3) * t224 * t233, -g(3) * t226 - t237 * t224, -g(1) * (-t230 * t219 + t234 * t222) - g(2) * (t234 * t219 + t230 * t222) - g(3) * (pkin(2) * t244 + t245 * t226 + pkin(6)), 0, 0, 0, 0, 0, -g(1) * t209 - g(2) * t208 - g(3) * t214, -g(1) * (-t213 * t228 + t232 * t243) - g(2) * (-t211 * t228 - t232 * t242) - g(3) * (-t216 * t228 + t226 * t232), 0, 0, 0, 0, 0, -g(1) * (t209 * t231 + t212 * t227) - g(2) * (t208 * t231 + t210 * t227) - g(3) * (t214 * t231 + t215 * t227), -g(1) * (-t209 * t227 + t212 * t231) - g(2) * (-t208 * t227 + t210 * t231) - g(3) * (-t214 * t227 + t215 * t231);];
U_reg = t1;
