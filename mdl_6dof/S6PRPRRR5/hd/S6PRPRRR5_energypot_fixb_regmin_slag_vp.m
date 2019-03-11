% Calculate minimal parameter regressor of potential energy for
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:43:35
% EndTime: 2019-03-08 20:43:35
% DurationCPUTime: 0.13s
% Computational Cost: add. (108->53), mult. (205->92), div. (0->0), fcn. (253->12), ass. (0->31)
t199 = sin(pkin(11));
t200 = sin(pkin(6));
t218 = t199 * t200;
t201 = cos(pkin(11));
t217 = t200 * t201;
t204 = sin(qJ(4));
t216 = t200 * t204;
t205 = sin(qJ(2));
t215 = t200 * t205;
t207 = cos(qJ(4));
t214 = t200 * t207;
t208 = cos(qJ(2));
t213 = t200 * t208;
t202 = cos(pkin(6));
t212 = t202 * t205;
t211 = t202 * t208;
t191 = t199 * t205 - t201 * t211;
t193 = t199 * t211 + t201 * t205;
t210 = -g(1) * t193 - g(2) * t191 + g(3) * t213;
t192 = t199 * t208 + t201 * t212;
t194 = -t199 * t212 + t201 * t208;
t209 = g(1) * t194 + g(2) * t192 + g(3) * t215;
t206 = cos(qJ(6));
t203 = sin(qJ(6));
t198 = qJ(4) + qJ(5);
t197 = cos(t198);
t196 = sin(t198);
t190 = -t196 * t213 + t197 * t202;
t189 = t191 * t196 - t197 * t217;
t188 = t193 * t196 + t197 * t218;
t1 = [-g(3) * qJ(1), 0, -t209, -t210, t209, t210, -g(1) * (pkin(1) * t201 + pkin(2) * t194 + pkin(7) * t218 + qJ(3) * t193) - g(2) * (pkin(1) * t199 + pkin(2) * t192 - pkin(7) * t217 + qJ(3) * t191) - g(3) * (t202 * pkin(7) + qJ(1) + (pkin(2) * t205 - qJ(3) * t208) * t200) 0, 0, 0, 0, 0, -g(1) * (t193 * t204 + t199 * t214) - g(2) * (t191 * t204 - t201 * t214) - g(3) * (t202 * t207 - t204 * t213) -g(1) * (t193 * t207 - t199 * t216) - g(2) * (t191 * t207 + t201 * t216) - g(3) * (-t202 * t204 - t207 * t213) 0, 0, 0, 0, 0, -g(1) * t188 - g(2) * t189 - g(3) * t190, -g(1) * (t193 * t197 - t196 * t218) - g(2) * (t191 * t197 + t196 * t217) - g(3) * (-t196 * t202 - t197 * t213) 0, 0, 0, 0, 0, -g(1) * (t188 * t206 + t194 * t203) - g(2) * (t189 * t206 + t192 * t203) - g(3) * (t190 * t206 + t203 * t215) -g(1) * (-t188 * t203 + t194 * t206) - g(2) * (-t189 * t203 + t192 * t206) - g(3) * (-t190 * t203 + t206 * t215);];
U_reg  = t1;
