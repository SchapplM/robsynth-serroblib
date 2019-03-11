% Calculate minimal parameter regressor of potential energy for
% S6RPRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:16:18
% EndTime: 2019-03-09 03:16:19
% DurationCPUTime: 0.10s
% Computational Cost: add. (169->62), mult. (170->81), div. (0->0), fcn. (177->10), ass. (0->40)
t215 = pkin(9) + qJ(3);
t210 = sin(t215);
t239 = g(3) * t210;
t217 = sin(pkin(9));
t238 = t217 * pkin(2) + pkin(6);
t218 = cos(pkin(10));
t206 = t218 * pkin(4) + pkin(3);
t212 = cos(t215);
t237 = t206 * t212;
t214 = pkin(10) + qJ(5);
t209 = sin(t214);
t222 = sin(qJ(1));
t236 = t222 * t209;
t211 = cos(t214);
t235 = t222 * t211;
t216 = sin(pkin(10));
t234 = t222 * t216;
t233 = t222 * t218;
t223 = cos(qJ(1));
t232 = t223 * t209;
t231 = t223 * t211;
t230 = t223 * t216;
t229 = t223 * t218;
t219 = cos(pkin(9));
t207 = t219 * pkin(2) + pkin(1);
t221 = -pkin(7) - qJ(2);
t228 = t222 * t207 + t223 * t221;
t227 = t223 * t207 - t222 * t221;
t226 = g(1) * t223 + g(2) * t222;
t225 = pkin(3) * t212 + qJ(4) * t210;
t197 = t212 * t236 + t231;
t199 = t212 * t232 - t235;
t224 = g(1) * t199 + g(2) * t197 + t209 * t239;
t220 = -pkin(8) - qJ(4);
t203 = g(1) * t222 - g(2) * t223;
t200 = t212 * t231 + t236;
t198 = t212 * t235 - t232;
t196 = -g(3) * t212 + t226 * t210;
t195 = -g(1) * t200 - g(2) * t198 - t211 * t239;
t1 = [0, -t226, t203, -g(3) * t217 - t226 * t219, -g(3) * t219 + t226 * t217, -t203, -g(1) * (t223 * pkin(1) + t222 * qJ(2)) - g(2) * (t222 * pkin(1) - t223 * qJ(2)) - g(3) * pkin(6), 0, 0, 0, 0, 0, -t226 * t212 - t239, t196, -g(1) * (t212 * t229 + t234) - g(2) * (t212 * t233 - t230) - t218 * t239, -g(1) * (-t212 * t230 + t233) - g(2) * (-t212 * t234 - t229) + t216 * t239, -t196, -g(1) * (t225 * t223 + t227) - g(2) * (t225 * t222 + t228) - g(3) * (t210 * pkin(3) - t212 * qJ(4) + t238) 0, 0, 0, 0, 0, t195, t224, t195, -t196, -t224, -g(1) * (pkin(4) * t234 + t200 * pkin(5) + t199 * qJ(6) + t223 * t237 + t227) - g(2) * (-pkin(4) * t230 + t198 * pkin(5) + t197 * qJ(6) + t222 * t237 + t228) - g(3) * (t212 * t220 + t238) + (-g(3) * (pkin(5) * t211 + qJ(6) * t209 + t206) + t226 * t220) * t210;];
U_reg  = t1;
