% Calculate minimal parameter regressor of potential energy for
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:23:09
% EndTime: 2019-03-09 18:23:09
% DurationCPUTime: 0.07s
% Computational Cost: add. (86->39), mult. (96->53), div. (0->0), fcn. (101->10), ass. (0->28)
t216 = qJ(2) + qJ(3);
t214 = cos(t216);
t234 = g(3) * t214;
t215 = qJ(5) + qJ(6);
t211 = sin(t215);
t219 = sin(qJ(1));
t233 = t219 * t211;
t213 = cos(t215);
t232 = t219 * t213;
t217 = sin(qJ(5));
t231 = t219 * t217;
t220 = cos(qJ(5));
t230 = t219 * t220;
t222 = cos(qJ(1));
t229 = t222 * t211;
t228 = t222 * t213;
t227 = t222 * t217;
t226 = t222 * t220;
t225 = g(1) * t222 + g(2) * t219;
t212 = sin(t216);
t221 = cos(qJ(2));
t224 = t221 * pkin(2) + pkin(3) * t214 + qJ(4) * t212 + pkin(1);
t223 = -pkin(8) - pkin(7);
t218 = sin(qJ(2));
t209 = g(1) * t219 - g(2) * t222;
t208 = g(3) * t212 + t225 * t214;
t207 = t225 * t212 - t234;
t1 = [0, -t225, t209, 0, 0, 0, 0, 0, -g(3) * t218 - t225 * t221, -g(3) * t221 + t225 * t218, 0, 0, 0, 0, 0, -t208, t207, -t209, t208, -t207, -g(3) * (t218 * pkin(2) + t212 * pkin(3) - t214 * qJ(4) + pkin(6)) + (-g(1) * t224 - g(2) * t223) * t222 + (g(1) * t223 - g(2) * t224) * t219, 0, 0, 0, 0, 0, -g(1) * (t212 * t227 + t230) - g(2) * (t212 * t231 - t226) + t217 * t234, -g(1) * (t212 * t226 - t231) - g(2) * (t212 * t230 + t227) + t220 * t234, 0, 0, 0, 0, 0, -g(1) * (t212 * t229 + t232) - g(2) * (t212 * t233 - t228) + t211 * t234, -g(1) * (t212 * t228 - t233) - g(2) * (t212 * t232 + t229) + t213 * t234;];
U_reg  = t1;
