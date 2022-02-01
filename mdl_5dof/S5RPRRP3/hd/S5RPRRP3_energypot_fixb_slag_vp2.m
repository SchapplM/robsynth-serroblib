% Calculate potential energy for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP3_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:29:29
% EndTime: 2022-01-23 09:29:30
% DurationCPUTime: 0.23s
% Computational Cost: add. (129->45), mult. (104->31), div. (0->0), fcn. (73->8), ass. (0->20)
t61 = -m(5) - m(6);
t65 = mrSges(5,2) + mrSges(6,2);
t64 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t62 = -m(3) - m(4);
t63 = -mrSges(2,1) + pkin(1) * (t61 + t62);
t47 = qJ(3) + qJ(4);
t39 = sin(t47);
t40 = cos(t47);
t49 = sin(qJ(3));
t51 = cos(qJ(3));
t60 = -m(4) * pkin(2) - t51 * mrSges(4,1) + t49 * mrSges(4,2) - mrSges(3,1) + t61 * (t51 * pkin(3) + pkin(2)) + t64 * t40 + t65 * t39;
t53 = -pkin(7) - pkin(6);
t58 = m(4) * pkin(6) - m(5) * t53 - m(6) * (-qJ(5) + t53) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t48 = qJ(2) + pkin(5);
t52 = cos(qJ(1));
t50 = sin(qJ(1));
t46 = qJ(1) + pkin(8);
t38 = cos(t46);
t37 = sin(t46);
t1 = (-m(2) * pkin(5) - t49 * mrSges(4,1) - t51 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t61 * (t49 * pkin(3) + t48) + t62 * t48 - t65 * t40 + t64 * t39) * g(3) + (-t52 * mrSges(2,2) + t60 * t37 + t58 * t38 + t63 * t50 - mrSges(1,2)) * g(2) + (t50 * mrSges(2,2) - t58 * t37 + t60 * t38 + t52 * t63 - mrSges(1,1)) * g(1);
U = t1;
