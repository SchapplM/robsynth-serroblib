% Calculate potential energy for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPPR6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_energypot_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR6_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:35
% EndTime: 2019-12-31 16:40:35
% DurationCPUTime: 0.24s
% Computational Cost: add. (70->41), mult. (120->38), div. (0->0), fcn. (101->6), ass. (0->18)
t69 = mrSges(3,2) - mrSges(4,3);
t45 = sin(pkin(6));
t46 = cos(pkin(6));
t68 = pkin(2) * t46 + qJ(3) * t45;
t67 = -m(5) * pkin(3) - mrSges(3,1) - mrSges(4,1);
t65 = -m(4) - m(5);
t48 = sin(qJ(1));
t64 = t68 * t48;
t47 = sin(qJ(4));
t49 = cos(qJ(4));
t52 = t45 * t47 + t46 * t49;
t53 = t45 * t49 - t46 * t47;
t63 = -mrSges(5,1) * t52 - mrSges(5,2) * t53 + t69 * t45 + t67 * t46 - mrSges(2,1);
t62 = mrSges(4,2) + mrSges(3,3) - mrSges(2,2) - mrSges(5,3);
t43 = t48 * pkin(1);
t50 = cos(qJ(1));
t58 = -t50 * qJ(2) + t43;
t1 = (-t53 * mrSges(5,1) + t52 * mrSges(5,2) - mrSges(1,3) - mrSges(2,3) + t65 * (t45 * pkin(2) - t46 * qJ(3) + pkin(4)) - t69 * t46 + t67 * t45 + (-m(2) - m(3)) * pkin(4)) * g(3) + (-mrSges(1,2) - m(3) * t58 - m(4) * (t58 + t64) - m(5) * (t43 + t64) + (-m(5) * (pkin(5) - qJ(2)) + t62) * t50 + t63 * t48) * g(2) + (-mrSges(1,1) + (m(5) * pkin(5) - t62) * t48 + (-m(3) + t65) * (t50 * pkin(1) + t48 * qJ(2)) + (t65 * t68 + t63) * t50) * g(1);
U = t1;
