% Calculate minimal parameter regressor of potential energy for
% S6RRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:50:33
% EndTime: 2019-03-09 11:50:33
% DurationCPUTime: 0.10s
% Computational Cost: add. (101->46), mult. (103->64), div. (0->0), fcn. (105->10), ass. (0->33)
t216 = qJ(4) + qJ(5);
t212 = cos(t216);
t221 = cos(qJ(4));
t202 = t221 * pkin(4) + pkin(5) * t212 + pkin(3);
t215 = qJ(2) + pkin(10);
t209 = sin(t215);
t210 = cos(t215);
t214 = -qJ(6) - pkin(9) - pkin(8);
t238 = t202 * t210 - t209 * t214;
t237 = g(3) * t209;
t219 = sin(qJ(2));
t236 = t219 * pkin(2) + pkin(6);
t211 = sin(t216);
t220 = sin(qJ(1));
t233 = t220 * t211;
t232 = t220 * t212;
t218 = sin(qJ(4));
t231 = t220 * t218;
t230 = t220 * t221;
t223 = cos(qJ(1));
t229 = t223 * t211;
t228 = t223 * t212;
t227 = t223 * t218;
t226 = t223 * t221;
t222 = cos(qJ(2));
t208 = t222 * pkin(2) + pkin(1);
t217 = -pkin(7) - qJ(3);
t225 = t220 * t208 + t223 * t217;
t224 = g(1) * t223 + g(2) * t220;
t206 = t223 * t208;
t204 = g(1) * t220 - g(2) * t223;
t203 = t218 * pkin(4) + pkin(5) * t211;
t1 = [0, -t224, t204, 0, 0, 0, 0, 0, -g(3) * t219 - t224 * t222, -g(3) * t222 + t224 * t219, -t204, -g(1) * (-t220 * t217 + t206) - g(2) * t225 - g(3) * t236, 0, 0, 0, 0, 0, -g(1) * (t210 * t226 + t231) - g(2) * (t210 * t230 - t227) - t221 * t237, -g(1) * (-t210 * t227 + t230) - g(2) * (-t210 * t231 - t226) + t218 * t237, 0, 0, 0, 0, 0, -g(1) * (t210 * t228 + t233) - g(2) * (t210 * t232 - t229) - t212 * t237, -g(1) * (-t210 * t229 + t232) - g(2) * (-t210 * t233 - t228) + t211 * t237, g(3) * t210 - t224 * t209, -g(1) * (t238 * t223 + t206) - g(2) * (-t223 * t203 + t225) - g(3) * (t209 * t202 + t210 * t214 + t236) + (-g(1) * (t203 - t217) - g(2) * t238) * t220;];
U_reg  = t1;
