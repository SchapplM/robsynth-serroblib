% Calculate potential energy for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:26
% EndTime: 2020-01-03 12:03:27
% DurationCPUTime: 0.35s
% Computational Cost: add. (135->51), mult. (104->38), div. (0->0), fcn. (73->10), ass. (0->22)
t63 = -m(3) - m(4);
t64 = pkin(1) * (m(5) + m(6) - t63) + mrSges(2,1);
t47 = pkin(9) + qJ(4);
t40 = qJ(5) + t47;
t35 = sin(t40);
t36 = cos(t40);
t50 = cos(pkin(9));
t37 = t50 * pkin(3) + pkin(2);
t38 = sin(t47);
t39 = cos(t47);
t49 = sin(pkin(9));
t61 = mrSges(3,1) + m(4) * pkin(2) + t50 * mrSges(4,1) - t49 * mrSges(4,2) + m(5) * t37 + t39 * mrSges(5,1) - t38 * mrSges(5,2) + m(6) * (pkin(4) * t39 + t37) + t36 * mrSges(6,1) - t35 * mrSges(6,2);
t51 = -pkin(7) - qJ(3);
t60 = m(4) * qJ(3) - m(5) * t51 - m(6) * (-pkin(8) + t51) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t54 = pkin(6) + pkin(5);
t58 = t49 * pkin(3) + t54;
t53 = cos(qJ(1));
t52 = sin(qJ(1));
t48 = qJ(1) + qJ(2);
t42 = cos(t48);
t41 = sin(t48);
t1 = (-t52 * mrSges(2,2) + t60 * t41 + t61 * t42 + t64 * t53 - mrSges(1,3)) * g(3) + (-t53 * mrSges(2,2) - t61 * t41 + t60 * t42 - t64 * t52 - mrSges(1,2)) * g(2) + (-mrSges(1,1) - m(2) * pkin(5) - mrSges(2,3) - mrSges(3,3) - t49 * mrSges(4,1) - t50 * mrSges(4,2) - m(5) * t58 - t38 * mrSges(5,1) - t39 * mrSges(5,2) - m(6) * (pkin(4) * t38 + t58) - t35 * mrSges(6,1) - t36 * mrSges(6,2) + t63 * t54) * g(1);
U = t1;
