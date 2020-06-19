% Calculate potential energy for
% S2RR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% m [3x1]
%   mass of all robot links (including the base)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S2RR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_energypot_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_energypot_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR3_energypot_fixb_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR3_energypot_fixb_slag_vp1: rSges has to be [3x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:21
% EndTime: 2020-06-19 09:14:21
% DurationCPUTime: 0.05s
% Computational Cost: add. (21->18), mult. (22->22), div. (0->0), fcn. (10->4), ass. (0->6)
t11 = cos(qJ(1));
t10 = sin(qJ(1));
t9 = qJ(1) + qJ(2);
t8 = cos(t9);
t7 = sin(t9);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t11 * rSges(2,1) - t10 * rSges(2,2)) + g(2) * (t10 * rSges(2,1) + t11 * rSges(2,2)) + g(3) * (pkin(2) + rSges(2,3))) - m(3) * (g(1) * (t8 * rSges(3,1) - t7 * rSges(3,2) + t11 * pkin(1)) + g(2) * (t7 * rSges(3,1) + t8 * rSges(3,2) + t10 * pkin(1)) + g(3) * (pkin(3) + pkin(2) + rSges(3,3)));
U = t1;
