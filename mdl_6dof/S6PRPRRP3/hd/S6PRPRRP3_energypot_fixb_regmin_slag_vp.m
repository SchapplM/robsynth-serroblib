% Calculate minimal parameter regressor of potential energy for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:07:44
% EndTime: 2019-03-08 20:07:44
% DurationCPUTime: 0.17s
% Computational Cost: add. (180->78), mult. (306->119), div. (0->0), fcn. (370->12), ass. (0->44)
t212 = sin(pkin(11));
t238 = pkin(3) * t212;
t213 = sin(pkin(10));
t216 = cos(pkin(10));
t221 = sin(qJ(2));
t217 = cos(pkin(6));
t223 = cos(qJ(2));
t229 = t217 * t223;
t196 = t213 * t221 - t216 * t229;
t220 = sin(qJ(5));
t237 = t196 * t220;
t198 = t213 * t229 + t216 * t221;
t236 = t198 * t220;
t214 = sin(pkin(6));
t235 = t213 * t214;
t234 = t214 * t216;
t233 = t214 * t221;
t232 = t214 * t223;
t231 = t217 * t212;
t230 = t217 * t221;
t228 = t220 * t223;
t227 = pkin(7) * t217 + qJ(1);
t226 = pkin(1) * t216 + pkin(7) * t235;
t197 = t213 * t223 + t216 * t230;
t211 = pkin(11) + qJ(4);
t206 = sin(t211);
t207 = cos(t211);
t190 = t197 * t206 + t207 * t234;
t199 = -t213 * t230 + t216 * t223;
t192 = t199 * t206 - t207 * t235;
t194 = t206 * t233 - t207 * t217;
t225 = g(1) * t192 + g(2) * t190 + g(3) * t194;
t224 = -g(1) * t198 - g(2) * t196 + g(3) * t232;
t222 = cos(qJ(5));
t219 = -pkin(8) - qJ(3);
t218 = -qJ(6) - pkin(9);
t215 = cos(pkin(11));
t208 = t213 * pkin(1);
t205 = pkin(5) * t222 + pkin(4);
t204 = pkin(3) * t215 + pkin(2);
t195 = t206 * t217 + t207 * t233;
t193 = t199 * t207 + t206 * t235;
t191 = t197 * t207 - t206 * t234;
t1 = [-g(3) * qJ(1), 0, -g(1) * t199 - g(2) * t197 - g(3) * t233, -t224, -g(1) * (t199 * t215 + t212 * t235) - g(2) * (t197 * t215 - t212 * t234) - g(3) * (t215 * t233 + t231) -g(1) * (-t199 * t212 + t215 * t235) - g(2) * (-t197 * t212 - t215 * t234) - g(3) * (-t212 * t233 + t215 * t217) t224, -g(1) * (pkin(2) * t199 + qJ(3) * t198 + t226) - g(2) * (pkin(2) * t197 - pkin(7) * t234 + qJ(3) * t196 + t208) - g(3) * ((pkin(2) * t221 - qJ(3) * t223) * t214 + t227) 0, 0, 0, 0, 0, -g(1) * t193 - g(2) * t191 - g(3) * t195, t225, 0, 0, 0, 0, 0, -g(1) * (t193 * t222 + t236) - g(2) * (t191 * t222 + t237) - g(3) * (t195 * t222 - t214 * t228) -g(1) * (-t193 * t220 + t198 * t222) - g(2) * (-t191 * t220 + t196 * t222) - g(3) * (-t195 * t220 - t222 * t232) -t225, -g(1) * (pkin(5) * t236 - t192 * t218 + t193 * t205 - t198 * t219 + t199 * t204 + t226) - g(2) * (pkin(5) * t237 - t190 * t218 + t191 * t205 - t196 * t219 + t197 * t204 + t208) - g(3) * (pkin(3) * t231 - t194 * t218 + t195 * t205 + t227) + (-g(1) * t213 * t238 - g(3) * (-pkin(5) * t228 + t204 * t221 + t219 * t223) - g(2) * (-pkin(7) - t238) * t216) * t214;];
U_reg  = t1;
