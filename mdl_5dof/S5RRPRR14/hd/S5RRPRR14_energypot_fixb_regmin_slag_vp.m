% Calculate minimal parameter regressor of potential energy for
% S5RRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR14_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR14_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:38:37
% EndTime: 2019-12-31 20:38:37
% DurationCPUTime: 0.13s
% Computational Cost: add. (105->53), mult. (198->91), div. (0->0), fcn. (245->12), ass. (0->30)
t215 = sin(pkin(5));
t219 = sin(qJ(2));
t232 = t215 * t219;
t220 = sin(qJ(1));
t231 = t215 * t220;
t222 = cos(qJ(2));
t230 = t215 * t222;
t223 = cos(qJ(1));
t229 = t215 * t223;
t228 = t220 * t219;
t227 = t220 * t222;
t226 = t223 * t219;
t225 = t223 * t222;
t217 = cos(pkin(5));
t206 = -t217 * t225 + t228;
t208 = t217 * t227 + t226;
t224 = -g(1) * t208 - g(2) * t206 + g(3) * t230;
t221 = cos(qJ(5));
t218 = sin(qJ(5));
t216 = cos(pkin(10));
t214 = sin(pkin(10));
t213 = pkin(10) + qJ(4);
t212 = cos(t213);
t211 = sin(t213);
t209 = -t217 * t228 + t225;
t207 = t217 * t226 + t227;
t205 = t217 * t211 + t212 * t232;
t204 = t209 * t212 + t211 * t231;
t203 = t207 * t212 - t211 * t229;
t1 = [0, -g(1) * t223 - g(2) * t220, g(1) * t220 - g(2) * t223, 0, 0, 0, 0, 0, -g(1) * t209 - g(2) * t207 - g(3) * t232, -t224, -g(1) * (t209 * t216 + t214 * t231) - g(2) * (t207 * t216 - t214 * t229) - g(3) * (t217 * t214 + t216 * t232), -g(1) * (-t209 * t214 + t216 * t231) - g(2) * (-t207 * t214 - t216 * t229) - g(3) * (-t214 * t232 + t217 * t216), t224, -g(1) * (t223 * pkin(1) + t209 * pkin(2) + pkin(7) * t231 + t208 * qJ(3)) - g(2) * (t220 * pkin(1) + t207 * pkin(2) - pkin(7) * t229 + t206 * qJ(3)) - g(3) * (t217 * pkin(7) + pkin(6) + (pkin(2) * t219 - qJ(3) * t222) * t215), 0, 0, 0, 0, 0, -g(1) * t204 - g(2) * t203 - g(3) * t205, -g(1) * (-t209 * t211 + t212 * t231) - g(2) * (-t207 * t211 - t212 * t229) - g(3) * (-t211 * t232 + t217 * t212), 0, 0, 0, 0, 0, -g(1) * (t204 * t221 + t208 * t218) - g(2) * (t203 * t221 + t206 * t218) - g(3) * (t205 * t221 - t218 * t230), -g(1) * (-t204 * t218 + t208 * t221) - g(2) * (-t203 * t218 + t206 * t221) - g(3) * (-t205 * t218 - t221 * t230);];
U_reg = t1;
