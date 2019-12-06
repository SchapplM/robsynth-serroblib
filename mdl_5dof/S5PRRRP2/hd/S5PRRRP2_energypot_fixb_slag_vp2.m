% Calculate potential energy for
% S5PRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRP2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP2_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP2_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:34
% EndTime: 2019-12-05 16:41:34
% DurationCPUTime: 0.23s
% Computational Cost: add. (139->47), mult. (99->38), div. (0->0), fcn. (68->8), ass. (0->21)
t64 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t63 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t62 = -m(5) - m(6);
t51 = sin(qJ(4));
t52 = cos(qJ(4));
t61 = -t63 * t51 + t64 * t52 - mrSges(4,1);
t60 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2);
t49 = sin(pkin(8));
t45 = t49 * pkin(1);
t50 = cos(pkin(8));
t46 = t50 * pkin(1);
t59 = pkin(5) + qJ(1);
t48 = pkin(8) + qJ(2);
t42 = sin(t48);
t58 = pkin(2) * t42 + t45;
t43 = cos(t48);
t57 = pkin(2) * t43 + t46;
t44 = qJ(3) + t48;
t41 = cos(t44);
t40 = sin(t44);
t1 = (-m(2) * qJ(1) - m(3) * t59 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - mrSges(4,3) + t63 * t52 + t64 * t51 + (-m(4) + t62) * (pkin(6) + t59)) * g(3) + (-m(3) * t45 - m(4) * t58 - t49 * mrSges(2,1) - t42 * mrSges(3,1) - t50 * mrSges(2,2) - t43 * mrSges(3,2) - mrSges(1,2) + t62 * (t40 * pkin(3) - t41 * pkin(7) + t58) - t60 * t41 + t61 * t40) * g(2) + (-m(3) * t46 - m(4) * t57 - t50 * mrSges(2,1) - t43 * mrSges(3,1) + t49 * mrSges(2,2) + t42 * mrSges(3,2) - mrSges(1,1) + t62 * (t41 * pkin(3) + t40 * pkin(7) + t57) + t61 * t41 + t60 * t40) * g(1);
U = t1;
