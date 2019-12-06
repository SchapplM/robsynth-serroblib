% Calculate potential energy for
% S5PRRRP3
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRP3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP3_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:27
% EndTime: 2019-12-05 16:43:27
% DurationCPUTime: 0.24s
% Computational Cost: add. (129->45), mult. (104->31), div. (0->0), fcn. (73->8), ass. (0->20)
t61 = -m(5) - m(6);
t65 = mrSges(5,2) + mrSges(6,2);
t64 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t62 = -m(3) - m(4);
t63 = -mrSges(2,1) + pkin(1) * (t61 + t62);
t47 = qJ(3) + qJ(4);
t39 = sin(t47);
t40 = cos(t47);
t51 = sin(qJ(3));
t52 = cos(qJ(3));
t60 = -m(4) * pkin(2) - t52 * mrSges(4,1) + t51 * mrSges(4,2) - mrSges(3,1) + t61 * (t52 * pkin(3) + pkin(2)) + t64 * t40 + t65 * t39;
t53 = -pkin(7) - pkin(6);
t58 = m(4) * pkin(6) - m(5) * t53 - m(6) * (-qJ(5) + t53) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t50 = pkin(5) + qJ(1);
t49 = cos(pkin(8));
t48 = sin(pkin(8));
t46 = pkin(8) + qJ(2);
t38 = cos(t46);
t37 = sin(t46);
t1 = (-m(2) * qJ(1) - t51 * mrSges(4,1) - t52 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t61 * (t51 * pkin(3) + t50) + t62 * t50 - t65 * t40 + t64 * t39) * g(3) + (-t49 * mrSges(2,2) + t60 * t37 + t58 * t38 + t48 * t63 - mrSges(1,2)) * g(2) + (t48 * mrSges(2,2) - t58 * t37 + t60 * t38 + t49 * t63 - mrSges(1,1)) * g(1);
U = t1;
