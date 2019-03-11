% Calculate minimal parameter regressor of potential energy for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:12:06
% EndTime: 2019-03-09 08:12:06
% DurationCPUTime: 0.13s
% Computational Cost: add. (120->52), mult. (130->69), div. (0->0), fcn. (129->10), ass. (0->37)
t229 = cos(qJ(1));
t222 = qJ(2) + pkin(9);
t218 = cos(t222);
t243 = t218 * t229;
t216 = sin(t222);
t245 = qJ(4) * t216;
t248 = pkin(3) * t243 + t229 * t245;
t247 = g(3) * t218;
t226 = sin(qJ(2));
t246 = t226 * pkin(2) + pkin(6);
t227 = sin(qJ(1));
t244 = t218 * t227;
t221 = pkin(10) + qJ(6);
t215 = sin(t221);
t242 = t227 * t215;
t217 = cos(t221);
t241 = t227 * t217;
t223 = sin(pkin(10));
t240 = t227 * t223;
t224 = cos(pkin(10));
t239 = t227 * t224;
t238 = t229 * t215;
t237 = t229 * t217;
t236 = t229 * t223;
t235 = t229 * t224;
t228 = cos(qJ(2));
t214 = t228 * pkin(2) + pkin(1);
t225 = -pkin(7) - qJ(3);
t234 = t227 * t214 + t229 * t225;
t211 = t229 * t214;
t233 = -t227 * t225 + t211;
t232 = pkin(3) * t244 + t227 * t245 + t234;
t231 = g(1) * t229 + g(2) * t227;
t230 = t216 * pkin(3) - t218 * qJ(4) + t246;
t207 = g(1) * t227 - g(2) * t229;
t204 = g(3) * t216 + t231 * t218;
t1 = [0, -t231, t207, 0, 0, 0, 0, 0, -g(3) * t226 - t231 * t228, -g(3) * t228 + t231 * t226, -t207, -g(1) * t233 - g(2) * t234 - g(3) * t246, -t207, t204, -t231 * t216 + t247, -g(1) * (t233 + t248) - g(2) * t232 - g(3) * t230, -g(1) * (t216 * t236 + t239) - g(2) * (t216 * t240 - t235) + t223 * t247, -g(1) * (t216 * t235 - t240) - g(2) * (t216 * t239 + t236) + t224 * t247, -t204, -g(1) * (qJ(5) * t243 + t211 + (pkin(4) - t225) * t227 + t248) - g(2) * (-t229 * pkin(4) + qJ(5) * t244 + t232) - g(3) * (t216 * qJ(5) + t230) 0, 0, 0, 0, 0, -g(1) * (t216 * t238 + t241) - g(2) * (t216 * t242 - t237) + t215 * t247, -g(1) * (t216 * t237 - t242) - g(2) * (t216 * t241 + t238) + t217 * t247;];
U_reg  = t1;
