% Calculate potential energy for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% m [6x1]
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
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:48:39
% EndTime: 2022-01-20 11:48:39
% DurationCPUTime: 0.22s
% Computational Cost: add. (129->45), mult. (104->31), div. (0->0), fcn. (73->8), ass. (0->20)
t61 = -m(5) - m(6);
t65 = mrSges(5,2) + mrSges(6,2);
t64 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t62 = -m(3) - m(4);
t63 = -mrSges(2,1) + pkin(1) * (t61 + t62);
t46 = qJ(3) + qJ(4);
t37 = sin(t46);
t39 = cos(t46);
t48 = sin(qJ(3));
t50 = cos(qJ(3));
t60 = -m(4) * pkin(2) - t50 * mrSges(4,1) + t48 * mrSges(4,2) - mrSges(3,1) + t61 * (t50 * pkin(3) + pkin(2)) + t64 * t39 + t65 * t37;
t52 = -pkin(8) - pkin(7);
t58 = m(4) * pkin(7) - m(5) * t52 + m(6) * (qJ(5) - t52) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t53 = pkin(6) + pkin(5);
t51 = cos(qJ(1));
t49 = sin(qJ(1));
t47 = qJ(1) + qJ(2);
t40 = cos(t47);
t38 = sin(t47);
t1 = (-m(2) * pkin(5) - t48 * mrSges(4,1) - t50 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t61 * (t48 * pkin(3) + t53) + t62 * t53 - t65 * t39 + t64 * t37) * g(3) + (-t51 * mrSges(2,2) + t60 * t38 + t58 * t40 + t49 * t63 - mrSges(1,2)) * g(2) + (t49 * mrSges(2,2) - t58 * t38 + t60 * t40 + t51 * t63 - mrSges(1,1)) * g(1);
U = t1;
