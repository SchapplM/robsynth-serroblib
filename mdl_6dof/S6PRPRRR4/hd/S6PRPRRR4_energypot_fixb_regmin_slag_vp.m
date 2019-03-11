% Calculate minimal parameter regressor of potential energy for
% S6PRPRRR4
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
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:38:55
% EndTime: 2019-03-08 20:38:56
% DurationCPUTime: 0.14s
% Computational Cost: add. (152->63), mult. (257->108), div. (0->0), fcn. (323->14), ass. (0->31)
t223 = sin(pkin(11));
t224 = sin(pkin(6));
t238 = t223 * t224;
t226 = cos(pkin(11));
t237 = t224 * t226;
t229 = sin(qJ(2));
t236 = t224 * t229;
t231 = cos(qJ(2));
t235 = t224 * t231;
t227 = cos(pkin(6));
t234 = t227 * t229;
t233 = t227 * t231;
t211 = t223 * t229 - t226 * t233;
t213 = t223 * t233 + t226 * t229;
t232 = -g(1) * t213 - g(2) * t211 + g(3) * t235;
t230 = cos(qJ(5));
t228 = sin(qJ(5));
t225 = cos(pkin(12));
t222 = sin(pkin(12));
t221 = qJ(5) + qJ(6);
t220 = pkin(12) + qJ(4);
t219 = cos(t221);
t218 = sin(t221);
t217 = cos(t220);
t216 = sin(t220);
t214 = -t223 * t234 + t226 * t231;
t212 = t223 * t231 + t226 * t234;
t210 = t227 * t216 + t217 * t236;
t209 = t214 * t217 + t216 * t238;
t208 = t212 * t217 - t216 * t237;
t1 = [-g(3) * qJ(1), 0, -g(1) * t214 - g(2) * t212 - g(3) * t236, -t232, -g(1) * (t214 * t225 + t222 * t238) - g(2) * (t212 * t225 - t222 * t237) - g(3) * (t227 * t222 + t225 * t236) -g(1) * (-t214 * t222 + t225 * t238) - g(2) * (-t212 * t222 - t225 * t237) - g(3) * (-t222 * t236 + t227 * t225) t232, -g(1) * (t226 * pkin(1) + t214 * pkin(2) + pkin(7) * t238 + t213 * qJ(3)) - g(2) * (t223 * pkin(1) + t212 * pkin(2) - pkin(7) * t237 + t211 * qJ(3)) - g(3) * (t227 * pkin(7) + qJ(1) + (pkin(2) * t229 - qJ(3) * t231) * t224) 0, 0, 0, 0, 0, -g(1) * t209 - g(2) * t208 - g(3) * t210, -g(1) * (-t214 * t216 + t217 * t238) - g(2) * (-t212 * t216 - t217 * t237) - g(3) * (-t216 * t236 + t227 * t217) 0, 0, 0, 0, 0, -g(1) * (t209 * t230 + t213 * t228) - g(2) * (t208 * t230 + t211 * t228) - g(3) * (t210 * t230 - t228 * t235) -g(1) * (-t209 * t228 + t213 * t230) - g(2) * (-t208 * t228 + t211 * t230) - g(3) * (-t210 * t228 - t230 * t235) 0, 0, 0, 0, 0, -g(1) * (t209 * t219 + t213 * t218) - g(2) * (t208 * t219 + t211 * t218) - g(3) * (t210 * t219 - t218 * t235) -g(1) * (-t209 * t218 + t213 * t219) - g(2) * (-t208 * t218 + t211 * t219) - g(3) * (-t210 * t218 - t219 * t235);];
U_reg  = t1;
