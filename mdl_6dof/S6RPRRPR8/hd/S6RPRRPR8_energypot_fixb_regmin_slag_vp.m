% Calculate minimal parameter regressor of potential energy for
% S6RPRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:24:50
% EndTime: 2019-03-09 05:24:50
% DurationCPUTime: 0.07s
% Computational Cost: add. (71->42), mult. (96->57), div. (0->0), fcn. (98->8), ass. (0->28)
t194 = cos(qJ(3));
t207 = g(3) * t194;
t185 = qJ(4) + pkin(10) + qJ(6);
t182 = sin(t185);
t192 = sin(qJ(1));
t206 = t192 * t182;
t183 = cos(t185);
t205 = t192 * t183;
t190 = sin(qJ(4));
t204 = t192 * t190;
t193 = cos(qJ(4));
t203 = t192 * t193;
t195 = cos(qJ(1));
t202 = t195 * t182;
t201 = t195 * t183;
t200 = t195 * t190;
t199 = t195 * t193;
t198 = pkin(4) * t190 + pkin(7);
t197 = g(1) * (t195 * pkin(1) + t192 * qJ(2));
t180 = g(1) * t192 - g(2) * t195;
t184 = t193 * pkin(4) + pkin(3);
t189 = -qJ(5) - pkin(8);
t191 = sin(qJ(3));
t196 = t184 * t191 + t189 * t194;
t187 = t192 * pkin(1);
t181 = g(1) * t195 + g(2) * t192;
t179 = -g(3) * t191 + t180 * t194;
t1 = [0, -t181, t180, t181, -t180, -t197 - g(2) * (-t195 * qJ(2) + t187) - g(3) * pkin(6), 0, 0, 0, 0, 0, -t180 * t191 - t207, -t179, 0, 0, 0, 0, 0, -g(1) * (t191 * t203 + t200) - g(2) * (-t191 * t199 + t204) - t193 * t207, -g(1) * (-t191 * t204 + t199) - g(2) * (t191 * t200 + t203) + t190 * t207, t179, -t197 - g(2) * t187 - g(3) * (t194 * t184 - t191 * t189 + pkin(2) + pkin(6)) + (-g(1) * t196 - g(2) * t198) * t192 + (-g(1) * t198 - g(2) * (-qJ(2) - t196)) * t195, 0, 0, 0, 0, 0, -g(1) * (t191 * t205 + t202) - g(2) * (-t191 * t201 + t206) - t183 * t207, -g(1) * (-t191 * t206 + t201) - g(2) * (t191 * t202 + t205) + t182 * t207;];
U_reg  = t1;
