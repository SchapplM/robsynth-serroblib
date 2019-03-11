% Calculate minimal parameter regressor of potential energy for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:51:37
% EndTime: 2019-03-09 02:51:37
% DurationCPUTime: 0.11s
% Computational Cost: add. (125->54), mult. (137->71), div. (0->0), fcn. (136->10), ass. (0->33)
t215 = pkin(9) + qJ(3);
t209 = sin(t215);
t219 = cos(pkin(9));
t239 = t219 * pkin(2) + qJ(4) * t209 + pkin(1);
t211 = cos(t215);
t238 = g(3) * t211;
t221 = sin(qJ(1));
t236 = t211 * t221;
t222 = cos(qJ(1));
t235 = t211 * t222;
t214 = pkin(10) + qJ(6);
t208 = sin(t214);
t234 = t221 * t208;
t210 = cos(t214);
t233 = t221 * t210;
t216 = sin(pkin(10));
t232 = t221 * t216;
t218 = cos(pkin(10));
t231 = t221 * t218;
t230 = t222 * t208;
t229 = t222 * t210;
t228 = t222 * t216;
t227 = t222 * t218;
t226 = pkin(3) * t235 + t239 * t222;
t220 = -pkin(7) - qJ(2);
t225 = pkin(3) * t236 + t222 * t220 + t239 * t221;
t224 = g(1) * t222 + g(2) * t221;
t217 = sin(pkin(9));
t223 = t217 * pkin(2) + t209 * pkin(3) - t211 * qJ(4) + pkin(6);
t202 = g(1) * t221 - g(2) * t222;
t197 = g(3) * t209 + t224 * t211;
t196 = t224 * t209 - t238;
t1 = [0, -t224, t202, -g(3) * t217 - t224 * t219, -g(3) * t219 + t224 * t217, -t202, -g(1) * (t222 * pkin(1) + t221 * qJ(2)) - g(2) * (t221 * pkin(1) - t222 * qJ(2)) - g(3) * pkin(6), 0, 0, 0, 0, 0, -t197, t196, -t202, t197, -t196, -g(1) * (-t221 * t220 + t226) - g(2) * t225 - g(3) * t223, -g(1) * (t209 * t228 + t231) - g(2) * (t209 * t232 - t227) + t216 * t238, -g(1) * (t209 * t227 - t232) - g(2) * (t209 * t231 + t228) + t218 * t238, -t197, -g(1) * (qJ(5) * t235 + (pkin(4) - t220) * t221 + t226) - g(2) * (-t222 * pkin(4) + qJ(5) * t236 + t225) - g(3) * (t209 * qJ(5) + t223) 0, 0, 0, 0, 0, -g(1) * (t209 * t230 + t233) - g(2) * (t209 * t234 - t229) + t208 * t238, -g(1) * (t209 * t229 - t234) - g(2) * (t209 * t233 + t230) + t210 * t238;];
U_reg  = t1;
