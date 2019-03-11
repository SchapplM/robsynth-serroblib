% Calculate minimal parameter regressor of potential energy for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:29:55
% EndTime: 2019-03-08 20:29:55
% DurationCPUTime: 0.18s
% Computational Cost: add. (146->55), mult. (332->103), div. (0->0), fcn. (432->14), ass. (0->37)
t244 = pkin(7) + qJ(3);
t227 = sin(pkin(6));
t232 = sin(qJ(4));
t243 = t227 * t232;
t233 = sin(qJ(2));
t242 = t227 * t233;
t235 = cos(qJ(4));
t241 = t227 * t235;
t230 = cos(pkin(6));
t240 = t230 * t233;
t236 = cos(qJ(2));
t239 = t230 * t236;
t225 = sin(pkin(12));
t228 = cos(pkin(12));
t238 = t236 * t225 + t233 * t228;
t237 = t233 * t225 - t236 * t228;
t234 = cos(qJ(5));
t231 = sin(qJ(5));
t229 = cos(pkin(11));
t226 = sin(pkin(11));
t224 = qJ(5) + qJ(6);
t223 = cos(t224);
t222 = sin(t224);
t221 = t236 * pkin(2) + pkin(1);
t218 = pkin(2) * t240 - t227 * t244;
t217 = t238 * t230;
t216 = t237 * t230;
t215 = t238 * t227;
t214 = t237 * t227;
t213 = t215 * t235 + t230 * t232;
t212 = -t226 * t217 - t229 * t237;
t211 = -t226 * t216 + t229 * t238;
t210 = t229 * t217 - t226 * t237;
t209 = t229 * t216 + t226 * t238;
t208 = t212 * t235 + t226 * t243;
t207 = t210 * t235 - t229 * t243;
t1 = [-g(3) * qJ(1), 0, -g(1) * (-t226 * t240 + t229 * t236) - g(2) * (t226 * t236 + t229 * t240) - g(3) * t242, -g(1) * (-t226 * t239 - t229 * t233) - g(2) * (-t226 * t233 + t229 * t239) - g(3) * t227 * t236, -g(1) * (-t226 * t218 + t229 * t221) - g(2) * (t229 * t218 + t226 * t221) - g(3) * (pkin(2) * t242 + t230 * t244 + qJ(1)) 0, 0, 0, 0, 0, -g(1) * t208 - g(2) * t207 - g(3) * t213, -g(1) * (-t212 * t232 + t226 * t241) - g(2) * (-t210 * t232 - t229 * t241) - g(3) * (-t215 * t232 + t230 * t235) 0, 0, 0, 0, 0, -g(1) * (t208 * t234 + t211 * t231) - g(2) * (t207 * t234 + t209 * t231) - g(3) * (t213 * t234 + t214 * t231) -g(1) * (-t208 * t231 + t211 * t234) - g(2) * (-t207 * t231 + t209 * t234) - g(3) * (-t213 * t231 + t214 * t234) 0, 0, 0, 0, 0, -g(1) * (t208 * t223 + t211 * t222) - g(2) * (t207 * t223 + t209 * t222) - g(3) * (t213 * t223 + t214 * t222) -g(1) * (-t208 * t222 + t211 * t223) - g(2) * (-t207 * t222 + t209 * t223) - g(3) * (-t213 * t222 + t214 * t223);];
U_reg  = t1;
