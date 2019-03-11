% Calculate minimal parameter regressor of potential energy for
% S6RPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:55:16
% EndTime: 2019-03-09 06:55:16
% DurationCPUTime: 0.05s
% Computational Cost: add. (75->22), mult. (59->34), div. (0->0), fcn. (60->12), ass. (0->23)
t191 = qJ(3) + qJ(4);
t189 = qJ(5) + t191;
t183 = sin(t189);
t204 = g(3) * t183;
t190 = qJ(1) + pkin(11);
t185 = sin(t190);
t192 = sin(qJ(6));
t203 = t185 * t192;
t195 = cos(qJ(6));
t202 = t185 * t195;
t186 = cos(t190);
t201 = t186 * t192;
t200 = t186 * t195;
t199 = g(1) * t186 + g(2) * t185;
t194 = sin(qJ(1));
t197 = cos(qJ(1));
t198 = -g(1) * t197 - g(2) * t194;
t196 = cos(qJ(3));
t193 = sin(qJ(3));
t188 = cos(t191);
t187 = sin(t191);
t184 = cos(t189);
t1 = [0, t198, g(1) * t194 - g(2) * t197, -g(3) * (qJ(2) + pkin(6)) + t198 * pkin(1), 0, 0, 0, 0, 0, -g(3) * t193 - t199 * t196, -g(3) * t196 + t199 * t193, 0, 0, 0, 0, 0, -g(3) * t187 - t199 * t188, -g(3) * t188 + t199 * t187, 0, 0, 0, 0, 0, -t199 * t184 - t204, -g(3) * t184 + t199 * t183, 0, 0, 0, 0, 0, -g(1) * (t184 * t200 + t203) - g(2) * (t184 * t202 - t201) - t195 * t204, -g(1) * (-t184 * t201 + t202) - g(2) * (-t184 * t203 - t200) + t192 * t204;];
U_reg  = t1;
