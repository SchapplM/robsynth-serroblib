% Calculate minimal parameter regressor of potential energy for
% S6RPRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:08:44
% EndTime: 2019-03-09 06:08:44
% DurationCPUTime: 0.08s
% Computational Cost: add. (108->42), mult. (98->54), div. (0->0), fcn. (96->10), ass. (0->25)
t198 = pkin(10) + qJ(3);
t196 = qJ(4) + t198;
t191 = sin(t196);
t213 = g(3) * t191;
t202 = sin(qJ(5));
t203 = sin(qJ(1));
t212 = t203 * t202;
t204 = cos(qJ(5));
t211 = t203 * t204;
t205 = cos(qJ(1));
t210 = t205 * t202;
t209 = t205 * t204;
t208 = pkin(5) * t202 + pkin(7) + pkin(8) + qJ(2);
t207 = g(1) * t205 + g(2) * t203;
t192 = cos(t196);
t193 = t204 * pkin(5) + pkin(4);
t195 = cos(t198);
t200 = cos(pkin(10));
t201 = -qJ(6) - pkin(9);
t206 = -t200 * pkin(2) - pkin(3) * t195 + t191 * t201 - t192 * t193 - pkin(1);
t199 = sin(pkin(10));
t194 = sin(t198);
t190 = g(1) * t203 - g(2) * t205;
t188 = -g(3) * t192 + t207 * t191;
t1 = [0, -t207, t190, -g(3) * t199 - t207 * t200, -g(3) * t200 + t207 * t199, -t190, -g(1) * (t205 * pkin(1) + t203 * qJ(2)) - g(2) * (t203 * pkin(1) - t205 * qJ(2)) - g(3) * pkin(6), 0, 0, 0, 0, 0, -g(3) * t194 - t207 * t195, -g(3) * t195 + t207 * t194, 0, 0, 0, 0, 0, -t207 * t192 - t213, t188, 0, 0, 0, 0, 0, -g(1) * (t192 * t209 + t212) - g(2) * (t192 * t211 - t210) - t204 * t213, -g(1) * (-t192 * t210 + t211) - g(2) * (-t192 * t212 - t209) + t202 * t213, -t188, -g(3) * (t199 * pkin(2) + pkin(3) * t194 + t191 * t193 + t192 * t201 + pkin(6)) + (g(1) * t206 + g(2) * t208) * t205 + (-g(1) * t208 + g(2) * t206) * t203;];
U_reg  = t1;
