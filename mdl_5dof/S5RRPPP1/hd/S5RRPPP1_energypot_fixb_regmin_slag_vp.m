% Calculate minimal parameter regressor of potential energy for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:24:30
% EndTime: 2019-12-31 19:24:30
% DurationCPUTime: 0.13s
% Computational Cost: add. (141->55), mult. (346->75), div. (0->0), fcn. (400->8), ass. (0->40)
t233 = sin(qJ(2));
t232 = cos(pkin(5));
t234 = sin(qJ(1));
t252 = t234 * t232;
t230 = sin(pkin(5));
t236 = cos(qJ(1));
t256 = t230 * t236;
t214 = t233 * t256 + t252;
t260 = t233 * t252 + t256;
t259 = t233 * pkin(2) + pkin(6);
t258 = t230 * t234;
t235 = cos(qJ(2));
t257 = t230 * t235;
t255 = t232 * t235;
t229 = sin(pkin(8));
t254 = t233 * t229;
t231 = cos(pkin(8));
t253 = t233 * t231;
t251 = t234 * t235;
t250 = t235 * t236;
t249 = t236 * t232;
t247 = t233 * t258;
t246 = qJ(3) * t257;
t244 = g(1) * t236 + g(2) * t234;
t243 = t236 * pkin(1) + pkin(2) * t250 + t234 * pkin(7) + t214 * qJ(3);
t211 = -t231 * t255 + t254;
t212 = t229 * t255 + t253;
t242 = t212 * pkin(3) + t211 * qJ(4) + t259;
t206 = t229 * t251 + t260 * t231;
t208 = -t231 * t258 + (t229 * t235 + t232 * t253) * t236;
t241 = g(1) * t208 + g(2) * t206 + g(3) * t211;
t207 = -t260 * t229 + t231 * t251;
t209 = t229 * t258 + t231 * t250 - t249 * t254;
t240 = g(1) * t209 + g(2) * t207 + g(3) * t212;
t239 = qJ(3) * t247 + pkin(2) * t251 + t234 * pkin(1) + (-qJ(3) * t232 - pkin(7)) * t236;
t238 = t209 * pkin(3) + t208 * qJ(4) + t243;
t237 = t207 * pkin(3) + t206 * qJ(4) + t239;
t213 = t247 - t249;
t203 = -g(1) * t214 - g(2) * t213 + g(3) * t257;
t1 = [0, -t244, g(1) * t234 - g(2) * t236, 0, 0, 0, 0, 0, -g(3) * t233 - t244 * t235, -g(3) * t235 + t244 * t233, -t240, t241, t203, -g(1) * t243 - g(2) * t239 - g(3) * (-t246 + t259), t203, t240, -t241, -g(1) * t238 - g(2) * t237 - g(3) * (t242 - t246), t203, -t241, -t240, -g(1) * (t214 * pkin(4) + t209 * qJ(5) + t238) - g(2) * (t213 * pkin(4) + t207 * qJ(5) + t237) - g(3) * (t212 * qJ(5) + (-pkin(4) - qJ(3)) * t257 + t242);];
U_reg = t1;
