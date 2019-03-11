% Calculate minimal parameter regressor of potential energy for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:35:19
% EndTime: 2019-03-09 03:35:19
% DurationCPUTime: 0.05s
% Computational Cost: add. (80->29), mult. (63->43), div. (0->0), fcn. (61->10), ass. (0->23)
t192 = qJ(3) + pkin(11) + qJ(5);
t187 = sin(t192);
t208 = g(3) * t187;
t207 = qJ(2) + pkin(6);
t193 = qJ(1) + pkin(10);
t190 = sin(t193);
t195 = sin(qJ(6));
t206 = t190 * t195;
t198 = cos(qJ(6));
t205 = t190 * t198;
t191 = cos(t193);
t204 = t191 * t195;
t203 = t191 * t198;
t202 = g(1) * t191 + g(2) * t190;
t197 = sin(qJ(1));
t200 = cos(qJ(1));
t201 = -g(1) * t200 - g(2) * t197;
t199 = cos(qJ(3));
t196 = sin(qJ(3));
t194 = -qJ(4) - pkin(7);
t189 = t199 * pkin(3) + pkin(2);
t188 = cos(t192);
t1 = [0, t201, g(1) * t197 - g(2) * t200, t201 * pkin(1) - g(3) * t207, 0, 0, 0, 0, 0, -g(3) * t196 - t202 * t199, -g(3) * t199 + t202 * t196, -g(1) * t190 + g(2) * t191, -g(1) * (t200 * pkin(1) + t191 * t189 - t190 * t194) - g(2) * (t197 * pkin(1) + t190 * t189 + t191 * t194) - g(3) * (t196 * pkin(3) + t207) 0, 0, 0, 0, 0, -t202 * t188 - t208, -g(3) * t188 + t202 * t187, 0, 0, 0, 0, 0, -g(1) * (t188 * t203 + t206) - g(2) * (t188 * t205 - t204) - t198 * t208, -g(1) * (-t188 * t204 + t205) - g(2) * (-t188 * t206 - t203) + t195 * t208;];
U_reg  = t1;
