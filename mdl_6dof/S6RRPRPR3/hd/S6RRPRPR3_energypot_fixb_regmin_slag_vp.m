% Calculate minimal parameter regressor of potential energy for
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:19:06
% EndTime: 2019-03-09 10:19:07
% DurationCPUTime: 0.10s
% Computational Cost: add. (98->44), mult. (98->62), div. (0->0), fcn. (100->10), ass. (0->32)
t237 = cos(qJ(4));
t225 = t237 * pkin(4) + pkin(3);
t231 = qJ(2) + pkin(10);
t227 = sin(t231);
t228 = cos(t231);
t232 = -qJ(5) - pkin(8);
t254 = t225 * t228 - t227 * t232;
t253 = g(3) * t227;
t235 = sin(qJ(2));
t252 = t235 * pkin(2) + pkin(6);
t229 = qJ(4) + pkin(11) + qJ(6);
t222 = sin(t229);
t236 = sin(qJ(1));
t249 = t236 * t222;
t223 = cos(t229);
t248 = t236 * t223;
t234 = sin(qJ(4));
t247 = t236 * t234;
t246 = t236 * t237;
t239 = cos(qJ(1));
t245 = t239 * t222;
t244 = t239 * t223;
t243 = t239 * t234;
t242 = t239 * t237;
t238 = cos(qJ(2));
t226 = t238 * pkin(2) + pkin(1);
t233 = -pkin(7) - qJ(3);
t241 = t236 * t226 + t239 * t233;
t240 = g(1) * t239 + g(2) * t236;
t221 = t239 * t226;
t219 = g(1) * t236 - g(2) * t239;
t1 = [0, -t240, t219, 0, 0, 0, 0, 0, -g(3) * t235 - t240 * t238, -g(3) * t238 + t240 * t235, -t219, -g(1) * (-t236 * t233 + t221) - g(2) * t241 - g(3) * t252, 0, 0, 0, 0, 0, -g(1) * (t228 * t242 + t247) - g(2) * (t228 * t246 - t243) - t237 * t253, -g(1) * (-t228 * t243 + t246) - g(2) * (-t228 * t247 - t242) + t234 * t253, g(3) * t228 - t240 * t227, -g(1) * (t254 * t239 + t221) - g(2) * (-pkin(4) * t243 + t241) - g(3) * (t227 * t225 + t228 * t232 + t252) + (-g(1) * (pkin(4) * t234 - t233) - g(2) * t254) * t236, 0, 0, 0, 0, 0, -g(1) * (t228 * t244 + t249) - g(2) * (t228 * t248 - t245) - t223 * t253, -g(1) * (-t228 * t245 + t248) - g(2) * (-t228 * t249 - t244) + t222 * t253;];
U_reg  = t1;
