% Calculate minimal parameter regressor of potential energy for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP13_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP13_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:01:32
% EndTime: 2019-03-09 13:01:32
% DurationCPUTime: 0.17s
% Computational Cost: add. (130->66), mult. (292->97), div. (0->0), fcn. (351->10), ass. (0->43)
t242 = sin(qJ(1));
t266 = g(1) * t242;
t246 = cos(qJ(1));
t265 = g(2) * t246;
t237 = cos(pkin(6));
t241 = sin(qJ(2));
t255 = t246 * t241;
t245 = cos(qJ(2));
t256 = t242 * t245;
t225 = t237 * t255 + t256;
t239 = sin(qJ(5));
t264 = t225 * t239;
t254 = t246 * t245;
t257 = t242 * t241;
t227 = -t237 * t257 + t254;
t263 = t227 * t239;
t236 = sin(pkin(6));
t262 = t236 * t241;
t261 = t236 * t242;
t260 = t236 * t245;
t259 = t236 * t246;
t258 = t239 * t241;
t253 = pkin(2) * t262 + t237 * pkin(8) + pkin(7);
t252 = -t265 + t266;
t224 = -t237 * t254 + t257;
t251 = t242 * pkin(1) + t225 * pkin(2) + t224 * qJ(3);
t226 = t237 * t256 + t255;
t250 = t246 * pkin(1) + t227 * pkin(2) + pkin(8) * t261 + t226 * qJ(3);
t240 = sin(qJ(4));
t244 = cos(qJ(4));
t215 = -t226 * t244 + t240 * t261;
t217 = t224 * t244 + t240 * t259;
t222 = t237 * t240 + t244 * t260;
t249 = g(1) * t215 - g(2) * t217 + g(3) * t222;
t248 = -g(1) * t226 - g(2) * t224 + g(3) * t260;
t247 = g(1) * t227 + g(2) * t225 + g(3) * t262;
t243 = cos(qJ(5));
t238 = -qJ(6) - pkin(10);
t232 = t243 * pkin(5) + pkin(4);
t223 = t237 * t244 - t240 * t260;
t218 = t224 * t240 - t244 * t259;
t216 = t226 * t240 + t244 * t261;
t1 = [0, -g(1) * t246 - g(2) * t242, t252, 0, 0, 0, 0, 0, -t247, -t248, -g(3) * t237 - t252 * t236, t247, t248, -g(1) * t250 - g(2) * (-pkin(8) * t259 + t251) - g(3) * (-qJ(3) * t260 + t253) 0, 0, 0, 0, 0, -g(1) * t216 - g(2) * t218 - g(3) * t223, t249, 0, 0, 0, 0, 0, -g(1) * (t216 * t243 + t263) - g(2) * (t218 * t243 + t264) - g(3) * (t223 * t243 + t236 * t258) -g(1) * (-t216 * t239 + t227 * t243) - g(2) * (-t218 * t239 + t225 * t243) - g(3) * (-t223 * t239 + t243 * t262) -t249, -g(1) * (pkin(5) * t263 + pkin(9) * t227 - t215 * t238 + t216 * t232 + t250) - g(2) * (pkin(5) * t264 + t225 * pkin(9) + t217 * t238 + t218 * t232 + t251) - g(3) * (t237 * pkin(3) - t222 * t238 + t223 * t232 + t253) + (-pkin(3) * t266 - g(3) * (pkin(5) * t258 + pkin(9) * t241 - qJ(3) * t245) - (-pkin(3) - pkin(8)) * t265) * t236;];
U_reg  = t1;
