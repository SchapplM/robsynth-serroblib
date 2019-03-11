% Calculate minimal parameter regressor of potential energy for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:43:56
% EndTime: 2019-03-09 08:43:56
% DurationCPUTime: 0.15s
% Computational Cost: add. (127->63), mult. (191->77), div. (0->0), fcn. (198->8), ass. (0->41)
t214 = sin(qJ(2));
t239 = qJ(3) * t214 + pkin(1);
t216 = cos(qJ(2));
t238 = g(3) * t216;
t237 = pkin(2) * t214 + pkin(6);
t210 = pkin(9) + qJ(5);
t204 = sin(t210);
t215 = sin(qJ(1));
t235 = t215 * t204;
t205 = cos(t210);
t234 = t215 * t205;
t211 = sin(pkin(9));
t233 = t215 * t211;
t212 = cos(pkin(9));
t232 = t215 * t212;
t231 = t215 * t216;
t217 = cos(qJ(1));
t230 = t216 * t217;
t229 = t217 * t204;
t228 = t217 * t205;
t227 = t217 * t211;
t226 = t217 * t212;
t225 = t214 * t233;
t224 = t214 * t227;
t223 = pkin(2) * t231 + t215 * t239;
t222 = pkin(2) * t230 + t215 * pkin(7) + t217 * t239;
t221 = -qJ(3) * t216 + t237;
t220 = g(1) * t217 + g(2) * t215;
t219 = -t217 * pkin(7) + t223;
t191 = -t214 * t228 + t235;
t193 = t214 * t234 + t229;
t218 = g(1) * t191 - g(2) * t193 + t205 * t238;
t213 = -pkin(8) - qJ(4);
t203 = pkin(4) * t212 + pkin(3);
t198 = g(1) * t215 - g(2) * t217;
t196 = g(3) * t214 + t216 * t220;
t195 = t214 * t220 - t238;
t194 = t214 * t235 - t228;
t192 = t214 * t229 + t234;
t190 = -g(1) * t192 - g(2) * t194 + t204 * t238;
t1 = [0, -t220, t198, 0, 0, 0, 0, 0, -t196, t195, -t198, t196, -t195, -g(1) * t222 - g(2) * t219 - g(3) * t221, -g(1) * (t224 + t232) - g(2) * (t225 - t226) + t211 * t238, -g(1) * (t214 * t226 - t233) - g(2) * (t214 * t232 + t227) + t212 * t238, -t196, -g(1) * (pkin(3) * t215 + qJ(4) * t230 + t222) - g(2) * (qJ(4) * t231 + (-pkin(3) - pkin(7)) * t217 + t223) - g(3) * (qJ(4) * t214 + t221) 0, 0, 0, 0, 0, t190, t218, t190, -t196, -t218, -g(1) * (pkin(4) * t224 + t192 * pkin(5) + t191 * qJ(6) + t215 * t203 + t222) - g(2) * (pkin(4) * t225 + t194 * pkin(5) - t193 * qJ(6) - t217 * t203 + t219) - g(3) * (-t214 * t213 + t237) + (-g(3) * (-pkin(4) * t211 - pkin(5) * t204 + qJ(6) * t205 - qJ(3)) + t220 * t213) * t216;];
U_reg  = t1;
