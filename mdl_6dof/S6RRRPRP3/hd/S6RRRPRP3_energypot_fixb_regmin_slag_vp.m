% Calculate minimal parameter regressor of potential energy for
% S6RRRPRP3
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
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:41:30
% EndTime: 2019-03-09 16:41:31
% DurationCPUTime: 0.12s
% Computational Cost: add. (164->57), mult. (161->74), div. (0->0), fcn. (171->10), ass. (0->39)
t217 = qJ(2) + qJ(3);
t213 = sin(t217);
t241 = g(3) * t213;
t221 = sin(qJ(2));
t240 = t221 * pkin(2) + pkin(6);
t219 = cos(pkin(10));
t207 = t219 * pkin(4) + pkin(3);
t214 = cos(t217);
t239 = t207 * t214;
t216 = pkin(10) + qJ(5);
t211 = sin(t216);
t222 = sin(qJ(1));
t238 = t222 * t211;
t212 = cos(t216);
t237 = t222 * t212;
t218 = sin(pkin(10));
t236 = t222 * t218;
t235 = t222 * t219;
t224 = cos(qJ(1));
t234 = t224 * t211;
t233 = t224 * t212;
t232 = t224 * t218;
t231 = t224 * t219;
t223 = cos(qJ(2));
t209 = t223 * pkin(2) + pkin(1);
t225 = -pkin(8) - pkin(7);
t230 = t222 * t209 + t224 * t225;
t229 = t224 * t209 - t222 * t225;
t228 = g(1) * t224 + g(2) * t222;
t227 = pkin(3) * t214 + qJ(4) * t213;
t200 = t214 * t238 + t233;
t202 = t214 * t234 - t237;
t226 = g(1) * t202 + g(2) * t200 + t211 * t241;
t220 = -pkin(9) - qJ(4);
t203 = t214 * t233 + t238;
t201 = t214 * t237 - t234;
t199 = -g(3) * t214 + t228 * t213;
t198 = -g(1) * t203 - g(2) * t201 - t212 * t241;
t1 = [0, -t228, g(1) * t222 - g(2) * t224, 0, 0, 0, 0, 0, -g(3) * t221 - t228 * t223, -g(3) * t223 + t228 * t221, 0, 0, 0, 0, 0, -t228 * t214 - t241, t199, -g(1) * (t214 * t231 + t236) - g(2) * (t214 * t235 - t232) - t219 * t241, -g(1) * (-t214 * t232 + t235) - g(2) * (-t214 * t236 - t231) + t218 * t241, -t199, -g(1) * (t227 * t224 + t229) - g(2) * (t227 * t222 + t230) - g(3) * (t213 * pkin(3) - t214 * qJ(4) + t240) 0, 0, 0, 0, 0, t198, t226, t198, -t199, -t226, -g(1) * (pkin(4) * t236 + t203 * pkin(5) + t202 * qJ(6) + t224 * t239 + t229) - g(2) * (-pkin(4) * t232 + t201 * pkin(5) + t200 * qJ(6) + t222 * t239 + t230) - g(3) * (t214 * t220 + t240) + (-g(3) * (pkin(5) * t212 + qJ(6) * t211 + t207) + t228 * t220) * t213;];
U_reg  = t1;
