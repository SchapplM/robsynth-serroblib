% Calculate potential energy for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRP6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP6_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:10
% EndTime: 2019-12-31 16:30:10
% DurationCPUTime: 0.27s
% Computational Cost: add. (79->34), mult. (143->29), div. (0->0), fcn. (130->6), ass. (0->15)
t61 = sin(qJ(3));
t63 = cos(qJ(3));
t77 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3);
t78 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1);
t88 = t77 * t61 + t78 * t63 - mrSges(3,1);
t85 = mrSges(3,2) - mrSges(5,2) - mrSges(4,3);
t81 = -m(4) - m(5);
t84 = -m(3) + t81;
t83 = t78 * t61 - t77 * t63 + mrSges(2,2) - mrSges(3,3);
t62 = sin(qJ(2));
t64 = cos(qJ(2));
t82 = -mrSges(2,1) + (t81 * pkin(5) + t85) * t62 + (t81 * pkin(2) + t88) * t64;
t60 = cos(pkin(6));
t59 = sin(pkin(6));
t1 = (-mrSges(1,3) - mrSges(2,3) + t81 * (t62 * pkin(2) - t64 * pkin(5) + qJ(1)) + (-m(2) - m(3)) * qJ(1) - t85 * t64 + t88 * t62) * g(3) + (-mrSges(1,2) + t84 * (t59 * pkin(1) - t60 * pkin(4)) - t83 * t60 + t82 * t59) * g(2) + (-mrSges(1,1) + t84 * (t60 * pkin(1) + t59 * pkin(4)) + t83 * t59 + t82 * t60) * g(1);
U = t1;
