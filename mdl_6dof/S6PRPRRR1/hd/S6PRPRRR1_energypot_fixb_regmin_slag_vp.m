% Calculate minimal parameter regressor of potential energy for
% S6PRPRRR1
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
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:25:05
% EndTime: 2019-03-08 20:25:05
% DurationCPUTime: 0.16s
% Computational Cost: add. (138->55), mult. (280->105), div. (0->0), fcn. (360->14), ass. (0->39)
t242 = pkin(7) + qJ(3);
t222 = sin(pkin(11));
t223 = sin(pkin(6));
t241 = t222 * t223;
t225 = cos(pkin(11));
t240 = t223 * t225;
t228 = sin(qJ(4));
t239 = t223 * t228;
t229 = sin(qJ(2));
t238 = t223 * t229;
t231 = cos(qJ(4));
t237 = t223 * t231;
t226 = cos(pkin(6));
t236 = t226 * t229;
t232 = cos(qJ(2));
t235 = t226 * t232;
t221 = sin(pkin(12));
t224 = cos(pkin(12));
t234 = t232 * t221 + t229 * t224;
t233 = t229 * t221 - t232 * t224;
t230 = cos(qJ(6));
t227 = sin(qJ(6));
t220 = qJ(4) + qJ(5);
t219 = cos(t220);
t218 = sin(t220);
t217 = t232 * pkin(2) + pkin(1);
t214 = pkin(2) * t236 - t242 * t223;
t213 = t234 * t226;
t212 = t233 * t226;
t211 = t234 * t223;
t210 = t233 * t223;
t209 = t211 * t219 + t226 * t218;
t208 = -t222 * t213 - t225 * t233;
t207 = -t222 * t212 + t225 * t234;
t206 = t225 * t213 - t222 * t233;
t205 = t225 * t212 + t222 * t234;
t204 = t208 * t219 + t218 * t241;
t203 = t206 * t219 - t218 * t240;
t1 = [-g(3) * qJ(1), 0, -g(1) * (-t222 * t236 + t225 * t232) - g(2) * (t222 * t232 + t225 * t236) - g(3) * t238, -g(1) * (-t222 * t235 - t225 * t229) - g(2) * (-t222 * t229 + t225 * t235) - g(3) * t223 * t232, -g(1) * (-t222 * t214 + t225 * t217) - g(2) * (t225 * t214 + t222 * t217) - g(3) * (pkin(2) * t238 + t242 * t226 + qJ(1)) 0, 0, 0, 0, 0, -g(1) * (t208 * t231 + t222 * t239) - g(2) * (t206 * t231 - t225 * t239) - g(3) * (t211 * t231 + t226 * t228) -g(1) * (-t208 * t228 + t222 * t237) - g(2) * (-t206 * t228 - t225 * t237) - g(3) * (-t211 * t228 + t226 * t231) 0, 0, 0, 0, 0, -g(1) * t204 - g(2) * t203 - g(3) * t209, -g(1) * (-t208 * t218 + t219 * t241) - g(2) * (-t206 * t218 - t219 * t240) - g(3) * (-t211 * t218 + t226 * t219) 0, 0, 0, 0, 0, -g(1) * (t204 * t230 + t207 * t227) - g(2) * (t203 * t230 + t205 * t227) - g(3) * (t209 * t230 + t210 * t227) -g(1) * (-t204 * t227 + t207 * t230) - g(2) * (-t203 * t227 + t205 * t230) - g(3) * (-t209 * t227 + t210 * t230);];
U_reg  = t1;
