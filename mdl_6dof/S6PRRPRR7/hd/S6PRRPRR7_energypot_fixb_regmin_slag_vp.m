% Calculate minimal parameter regressor of potential energy for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:34:10
% EndTime: 2019-03-08 22:34:11
% DurationCPUTime: 0.10s
% Computational Cost: add. (132->60), mult. (293->97), div. (0->0), fcn. (369->12), ass. (0->34)
t223 = sin(pkin(6));
t241 = pkin(7) * t223;
t227 = sin(qJ(3));
t240 = t223 * t227;
t228 = sin(qJ(2));
t239 = t223 * t228;
t230 = cos(qJ(3));
t238 = t223 * t230;
t231 = cos(qJ(2));
t237 = t223 * t231;
t225 = cos(pkin(6));
t236 = t225 * t228;
t235 = t225 * t231;
t222 = sin(pkin(11));
t224 = cos(pkin(11));
t211 = t222 * t231 + t224 * t236;
t206 = t211 * t227 + t224 * t238;
t213 = -t222 * t236 + t224 * t231;
t208 = t213 * t227 - t222 * t238;
t214 = -t225 * t230 + t227 * t239;
t234 = g(1) * t208 + g(2) * t206 + g(3) * t214;
t207 = t211 * t230 - t224 * t240;
t209 = t213 * t230 + t222 * t240;
t215 = t225 * t227 + t228 * t238;
t233 = g(1) * t209 + g(2) * t207 + g(3) * t215;
t210 = t222 * t228 - t224 * t235;
t212 = t222 * t235 + t224 * t228;
t232 = -g(1) * t212 - g(2) * t210 + g(3) * t237;
t229 = cos(qJ(5));
t226 = sin(qJ(5));
t221 = qJ(5) + qJ(6);
t220 = cos(t221);
t219 = sin(t221);
t1 = [-g(3) * qJ(1), 0, -g(1) * t213 - g(2) * t211 - g(3) * t239, -t232, 0, 0, 0, 0, 0, -t233, t234, t232, t233, -t234, -g(1) * (t224 * pkin(1) + t213 * pkin(2) + t209 * pkin(3) + t212 * pkin(8) + t208 * qJ(4) + t222 * t241) - g(2) * (t222 * pkin(1) + t211 * pkin(2) + t207 * pkin(3) + t210 * pkin(8) + t206 * qJ(4) - t224 * t241) - g(3) * (t215 * pkin(3) + t225 * pkin(7) + t214 * qJ(4) + qJ(1) + (pkin(2) * t228 - pkin(8) * t231) * t223) 0, 0, 0, 0, 0, -g(1) * (t208 * t226 + t212 * t229) - g(2) * (t206 * t226 + t210 * t229) - g(3) * (t214 * t226 - t229 * t237) -g(1) * (t208 * t229 - t212 * t226) - g(2) * (t206 * t229 - t210 * t226) - g(3) * (t214 * t229 + t226 * t237) 0, 0, 0, 0, 0, -g(1) * (t208 * t219 + t212 * t220) - g(2) * (t206 * t219 + t210 * t220) - g(3) * (t214 * t219 - t220 * t237) -g(1) * (t208 * t220 - t212 * t219) - g(2) * (t206 * t220 - t210 * t219) - g(3) * (t214 * t220 + t219 * t237);];
U_reg  = t1;
