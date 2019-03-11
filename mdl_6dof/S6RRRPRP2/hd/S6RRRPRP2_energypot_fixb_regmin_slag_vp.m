% Calculate minimal parameter regressor of potential energy for
% S6RRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:37:17
% EndTime: 2019-03-09 16:37:17
% DurationCPUTime: 0.07s
% Computational Cost: add. (139->42), mult. (127->56), div. (0->0), fcn. (133->10), ass. (0->32)
t222 = qJ(2) + qJ(3);
t217 = pkin(10) + t222;
t214 = cos(t217);
t239 = pkin(4) * t214;
t213 = sin(t217);
t238 = g(3) * t213;
t223 = sin(qJ(5));
t225 = sin(qJ(1));
t237 = t225 * t223;
t226 = cos(qJ(5));
t236 = t225 * t226;
t228 = cos(qJ(1));
t235 = t228 * t223;
t234 = t228 * t226;
t219 = cos(t222);
t227 = cos(qJ(2));
t210 = t227 * pkin(2) + pkin(3) * t219 + pkin(1);
t221 = -qJ(4) - pkin(8) - pkin(7);
t233 = t225 * t210 + t228 * t221;
t218 = sin(t222);
t224 = sin(qJ(2));
t232 = t224 * pkin(2) + pkin(3) * t218 + pkin(6);
t231 = t228 * t210 - t225 * t221;
t230 = g(1) * t228 + g(2) * t225;
t204 = t214 * t237 + t234;
t206 = t214 * t235 - t236;
t229 = g(1) * t206 + g(2) * t204 + t223 * t238;
t211 = g(1) * t225 - g(2) * t228;
t207 = t214 * t234 + t237;
t205 = t214 * t236 - t235;
t203 = -g(1) * t207 - g(2) * t205 - t226 * t238;
t1 = [0, -t230, t211, 0, 0, 0, 0, 0, -g(3) * t224 - t230 * t227, -g(3) * t227 + t230 * t224, 0, 0, 0, 0, 0, -g(3) * t218 - t230 * t219, -g(3) * t219 + t230 * t218, -t211, -g(1) * t231 - g(2) * t233 - g(3) * t232, 0, 0, 0, 0, 0, t203, t229, t203, g(3) * t214 - t230 * t213, -t229, -g(1) * (t207 * pkin(5) + t206 * qJ(6) + t228 * t239 + t231) - g(2) * (t205 * pkin(5) + t204 * qJ(6) + t225 * t239 + t233) - g(3) * (-t214 * pkin(9) + t232) + (-g(3) * (pkin(5) * t226 + qJ(6) * t223 + pkin(4)) - t230 * pkin(9)) * t213;];
U_reg  = t1;
