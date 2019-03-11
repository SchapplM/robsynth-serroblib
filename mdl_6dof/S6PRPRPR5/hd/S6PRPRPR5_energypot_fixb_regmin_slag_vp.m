% Calculate minimal parameter regressor of potential energy for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:45:02
% EndTime: 2019-03-08 19:45:02
% DurationCPUTime: 0.16s
% Computational Cost: add. (186->74), mult. (318->114), div. (0->0), fcn. (387->12), ass. (0->40)
t212 = sin(pkin(11));
t234 = pkin(3) * t212;
t213 = sin(pkin(10));
t214 = sin(pkin(6));
t233 = t213 * t214;
t216 = cos(pkin(10));
t232 = t214 * t216;
t220 = sin(qJ(2));
t231 = t214 * t220;
t222 = cos(qJ(2));
t230 = t214 * t222;
t217 = cos(pkin(6));
t229 = t217 * t212;
t228 = t217 * t220;
t227 = t217 * t222;
t226 = t217 * pkin(7) + qJ(1);
t225 = t216 * pkin(1) + pkin(7) * t233;
t198 = t213 * t222 + t216 * t228;
t211 = pkin(11) + qJ(4);
t206 = sin(t211);
t207 = cos(t211);
t191 = t198 * t206 + t207 * t232;
t200 = -t213 * t228 + t216 * t222;
t193 = t200 * t206 - t207 * t233;
t195 = t206 * t231 - t217 * t207;
t224 = g(1) * t193 + g(2) * t191 + g(3) * t195;
t192 = t198 * t207 - t206 * t232;
t194 = t200 * t207 + t206 * t233;
t196 = t217 * t206 + t207 * t231;
t223 = g(1) * t194 + g(2) * t192 + g(3) * t196;
t197 = t213 * t220 - t216 * t227;
t199 = t213 * t227 + t216 * t220;
t190 = -g(1) * t199 - g(2) * t197 + g(3) * t230;
t221 = cos(qJ(6));
t219 = sin(qJ(6));
t218 = -pkin(8) - qJ(3);
t215 = cos(pkin(11));
t208 = t213 * pkin(1);
t205 = t215 * pkin(3) + pkin(2);
t1 = [-g(3) * qJ(1), 0, -g(1) * t200 - g(2) * t198 - g(3) * t231, -t190, -g(1) * (t200 * t215 + t212 * t233) - g(2) * (t198 * t215 - t212 * t232) - g(3) * (t215 * t231 + t229) -g(1) * (-t200 * t212 + t215 * t233) - g(2) * (-t198 * t212 - t215 * t232) - g(3) * (-t212 * t231 + t217 * t215) t190, -g(1) * (t200 * pkin(2) + t199 * qJ(3) + t225) - g(2) * (t198 * pkin(2) - pkin(7) * t232 + t197 * qJ(3) + t208) - g(3) * ((pkin(2) * t220 - qJ(3) * t222) * t214 + t226) 0, 0, 0, 0, 0, -t223, t224, t190, t223, -t224, -g(1) * (t194 * pkin(4) + t193 * qJ(5) - t199 * t218 + t200 * t205 + t225) - g(2) * (t192 * pkin(4) + t191 * qJ(5) - t197 * t218 + t198 * t205 + t208) - g(3) * (pkin(3) * t229 + t196 * pkin(4) + t195 * qJ(5) + t226) + (-g(1) * t213 * t234 - g(3) * (t205 * t220 + t218 * t222) - g(2) * (-pkin(7) - t234) * t216) * t214, 0, 0, 0, 0, 0, -g(1) * (t193 * t219 + t199 * t221) - g(2) * (t191 * t219 + t197 * t221) - g(3) * (t195 * t219 - t221 * t230) -g(1) * (t193 * t221 - t199 * t219) - g(2) * (t191 * t221 - t197 * t219) - g(3) * (t195 * t221 + t219 * t230);];
U_reg  = t1;
