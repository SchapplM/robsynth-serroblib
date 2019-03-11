% Calculate minimal parameter regressor of potential energy for
% S6PRPRRR3
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
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:34:11
% EndTime: 2019-03-08 20:34:11
% DurationCPUTime: 0.14s
% Computational Cost: add. (154->63), mult. (231->108), div. (0->0), fcn. (287->14), ass. (0->31)
t222 = sin(pkin(11));
t223 = sin(pkin(6));
t237 = t222 * t223;
t225 = cos(pkin(11));
t236 = t223 * t225;
t228 = sin(qJ(2));
t235 = t223 * t228;
t230 = cos(qJ(2));
t234 = t223 * t230;
t226 = cos(pkin(6));
t233 = t226 * t228;
t232 = t226 * t230;
t220 = pkin(12) + qJ(4);
t210 = t222 * t228 - t225 * t232;
t212 = t222 * t232 + t225 * t228;
t231 = -g(1) * t212 - g(2) * t210 + g(3) * t234;
t229 = cos(qJ(6));
t227 = sin(qJ(6));
t224 = cos(pkin(12));
t221 = sin(pkin(12));
t219 = qJ(5) + t220;
t218 = cos(t220);
t217 = sin(t220);
t216 = cos(t219);
t215 = sin(t219);
t213 = -t222 * t233 + t225 * t230;
t211 = t222 * t230 + t225 * t233;
t209 = t226 * t215 + t216 * t235;
t208 = t213 * t216 + t215 * t237;
t207 = t211 * t216 - t215 * t236;
t1 = [-g(3) * qJ(1), 0, -g(1) * t213 - g(2) * t211 - g(3) * t235, -t231, -g(1) * (t213 * t224 + t221 * t237) - g(2) * (t211 * t224 - t221 * t236) - g(3) * (t226 * t221 + t224 * t235) -g(1) * (-t213 * t221 + t224 * t237) - g(2) * (-t211 * t221 - t224 * t236) - g(3) * (-t221 * t235 + t226 * t224) t231, -g(1) * (t225 * pkin(1) + t213 * pkin(2) + pkin(7) * t237 + t212 * qJ(3)) - g(2) * (t222 * pkin(1) + t211 * pkin(2) - pkin(7) * t236 + t210 * qJ(3)) - g(3) * (t226 * pkin(7) + qJ(1) + (pkin(2) * t228 - qJ(3) * t230) * t223) 0, 0, 0, 0, 0, -g(1) * (t213 * t218 + t217 * t237) - g(2) * (t211 * t218 - t217 * t236) - g(3) * (t226 * t217 + t218 * t235) -g(1) * (-t213 * t217 + t218 * t237) - g(2) * (-t211 * t217 - t218 * t236) - g(3) * (-t217 * t235 + t226 * t218) 0, 0, 0, 0, 0, -g(1) * t208 - g(2) * t207 - g(3) * t209, -g(1) * (-t213 * t215 + t216 * t237) - g(2) * (-t211 * t215 - t216 * t236) - g(3) * (-t215 * t235 + t226 * t216) 0, 0, 0, 0, 0, -g(1) * (t208 * t229 + t212 * t227) - g(2) * (t207 * t229 + t210 * t227) - g(3) * (t209 * t229 - t227 * t234) -g(1) * (-t208 * t227 + t212 * t229) - g(2) * (-t207 * t227 + t210 * t229) - g(3) * (-t209 * t227 - t229 * t234);];
U_reg  = t1;
