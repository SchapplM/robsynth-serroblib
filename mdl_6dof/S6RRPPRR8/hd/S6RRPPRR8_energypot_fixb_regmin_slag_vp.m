% Calculate minimal parameter regressor of potential energy for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:26:07
% EndTime: 2019-03-09 09:26:07
% DurationCPUTime: 0.15s
% Computational Cost: add. (104->54), mult. (210->81), div. (0->0), fcn. (244->10), ass. (0->28)
t230 = sin(qJ(2));
t245 = qJ(3) * t230 + pkin(1);
t244 = g(3) * t230;
t231 = sin(qJ(1));
t233 = cos(qJ(2));
t242 = t231 * t233;
t227 = sin(pkin(10));
t234 = cos(qJ(1));
t241 = t234 * t227;
t228 = cos(pkin(10));
t240 = t234 * t228;
t239 = t231 * pkin(7) + (pkin(2) * t233 + t245) * t234;
t238 = t230 * pkin(2) - t233 * qJ(3) + pkin(6);
t237 = g(1) * t234 + g(2) * t231;
t236 = pkin(2) * t242 - t234 * pkin(7) + t245 * t231;
t210 = t227 * t242 + t240;
t212 = -t231 * t228 + t233 * t241;
t235 = g(1) * t212 + g(2) * t210 + t227 * t244;
t232 = cos(qJ(5));
t229 = sin(qJ(5));
t226 = qJ(5) + qJ(6);
t220 = cos(t226);
t219 = sin(t226);
t213 = t231 * t227 + t233 * t240;
t211 = t228 * t242 - t241;
t209 = -g(3) * t233 + t237 * t230;
t208 = -g(1) * t213 - g(2) * t211 - t228 * t244;
t1 = [0, -t237, g(1) * t231 - g(2) * t234, 0, 0, 0, 0, 0, -t237 * t233 - t244, t209, t208, t235, -t209, -g(1) * t239 - g(2) * t236 - g(3) * t238, t208, -t209, -t235, -g(1) * (t213 * pkin(3) + t212 * qJ(4) + t239) - g(2) * (t211 * pkin(3) + t210 * qJ(4) + t236) - g(3) * ((pkin(3) * t228 + qJ(4) * t227) * t230 + t238) 0, 0, 0, 0, 0, -g(1) * (t212 * t229 + t213 * t232) - g(2) * (t210 * t229 + t211 * t232) - (t227 * t229 + t228 * t232) * t244, -g(1) * (t212 * t232 - t213 * t229) - g(2) * (t210 * t232 - t211 * t229) - (t227 * t232 - t228 * t229) * t244, 0, 0, 0, 0, 0, -g(1) * (t212 * t219 + t213 * t220) - g(2) * (t210 * t219 + t211 * t220) - (t219 * t227 + t220 * t228) * t244, -g(1) * (t212 * t220 - t213 * t219) - g(2) * (t210 * t220 - t211 * t219) - (-t219 * t228 + t220 * t227) * t244;];
U_reg  = t1;
