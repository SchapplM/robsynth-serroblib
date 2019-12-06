% Calculate kinetic energy for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR6_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 19:00:00
% EndTime: 2019-12-05 19:00:01
% DurationCPUTime: 0.25s
% Computational Cost: add. (412->71), mult. (579->120), div. (0->0), fcn. (346->8), ass. (0->33)
t84 = qJD(1) + qJD(2);
t99 = t84 / 0.2e1;
t88 = sin(qJ(2));
t97 = qJD(1) * pkin(1);
t80 = t84 * pkin(7) + t88 * t97;
t98 = t80 * mrSges(4,3);
t87 = sin(qJ(3));
t95 = pkin(8) * t84 + t80;
t74 = qJD(3) * pkin(3) - t95 * t87;
t91 = cos(qJ(3));
t75 = t95 * t91;
t86 = sin(qJ(4));
t90 = cos(qJ(4));
t67 = t86 * t74 + t90 * t75;
t83 = qJD(3) + qJD(4);
t92 = cos(qJ(2));
t96 = t92 * t97;
t66 = t90 * t74 - t86 * t75;
t78 = -t96 + (-pkin(3) * t91 - pkin(2)) * t84;
t89 = cos(qJ(5));
t85 = sin(qJ(5));
t82 = qJD(5) + t83;
t81 = -t84 * pkin(2) - t96;
t77 = (t86 * t91 + t87 * t90) * t84;
t76 = (-t86 * t87 + t90 * t91) * t84;
t70 = -t76 * pkin(4) + t78;
t69 = t85 * t76 + t89 * t77;
t68 = t89 * t76 - t85 * t77;
t65 = t76 * pkin(9) + t67;
t64 = t83 * pkin(4) - t77 * pkin(9) + t66;
t63 = t85 * t64 + t89 * t65;
t62 = t89 * t64 - t85 * t65;
t1 = m(4) * (t81 ^ 2 + (t87 ^ 2 + t91 ^ 2) * t80 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t78 ^ 2) / 0.2e1 + m(6) * (t62 ^ 2 + t63 ^ 2 + t70 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t88 ^ 2 + t92 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t66 * mrSges(5,1) - t67 * mrSges(5,2) + Ifges(5,3) * t83 / 0.2e1) * t83 + (t62 * mrSges(6,1) - t63 * mrSges(6,2) + Ifges(6,3) * t82 / 0.2e1) * t82 + (Ifges(4,3) * qJD(3) / 0.2e1 + (-t87 * mrSges(4,1) - t91 * mrSges(4,2)) * t80) * qJD(3) + (t78 * mrSges(5,2) - t66 * mrSges(5,3) + Ifges(5,5) * t83 + Ifges(5,1) * t77 / 0.2e1) * t77 + (t70 * mrSges(6,2) - t62 * mrSges(6,3) + Ifges(6,5) * t82 + Ifges(6,1) * t69 / 0.2e1) * t69 + (-t78 * mrSges(5,1) + t67 * mrSges(5,3) + Ifges(5,4) * t77 + Ifges(5,6) * t83 + Ifges(5,2) * t76 / 0.2e1) * t76 + (-t70 * mrSges(6,1) + t63 * mrSges(6,3) + Ifges(6,4) * t69 + Ifges(6,6) * t82 + Ifges(6,2) * t68 / 0.2e1) * t68 + (Ifges(3,3) * t99 + (mrSges(3,1) * t92 - mrSges(3,2) * t88) * t97 + (-t81 * mrSges(4,1) + Ifges(4,6) * qJD(3) + (Ifges(4,2) * t99 + t98) * t91) * t91 + (Ifges(4,4) * t91 * t84 + t81 * mrSges(4,2) + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t99 + t98) * t87) * t87) * t84;
T = t1;
