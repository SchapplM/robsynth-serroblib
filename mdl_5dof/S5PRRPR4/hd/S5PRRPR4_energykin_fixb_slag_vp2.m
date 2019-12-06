% Calculate kinetic energy for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR4_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR4_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR4_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:21:23
% EndTime: 2019-12-05 16:21:23
% DurationCPUTime: 0.27s
% Computational Cost: add. (231->70), mult. (539->116), div. (0->0), fcn. (346->8), ass. (0->30)
t87 = sin(qJ(2));
t80 = qJD(2) * pkin(6) + qJD(1) * t87;
t95 = t80 * mrSges(4,3);
t94 = qJD(2) / 0.2e1;
t86 = sin(qJ(3));
t92 = qJ(4) * qJD(2) + t80;
t74 = qJD(3) * pkin(3) - t92 * t86;
t89 = cos(qJ(3));
t75 = t92 * t89;
t83 = sin(pkin(9));
t84 = cos(pkin(9));
t67 = t83 * t74 + t84 * t75;
t90 = cos(qJ(2));
t93 = qJD(1) * t90;
t66 = t84 * t74 - t75 * t83;
t78 = -t93 + qJD(4) + (-pkin(3) * t89 - pkin(2)) * qJD(2);
t88 = cos(qJ(5));
t85 = sin(qJ(5));
t82 = qJD(3) + qJD(5);
t81 = -qJD(2) * pkin(2) - t93;
t77 = (t83 * t89 + t84 * t86) * qJD(2);
t76 = (-t83 * t86 + t84 * t89) * qJD(2);
t70 = -pkin(4) * t76 + t78;
t69 = t76 * t85 + t77 * t88;
t68 = t76 * t88 - t77 * t85;
t65 = pkin(7) * t76 + t67;
t64 = qJD(3) * pkin(4) - pkin(7) * t77 + t66;
t63 = t64 * t85 + t65 * t88;
t62 = t64 * t88 - t65 * t85;
t1 = m(4) * (t81 ^ 2 + (t86 ^ 2 + t89 ^ 2) * t80 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t78 ^ 2) / 0.2e1 + m(6) * (t62 ^ 2 + t63 ^ 2 + t70 ^ 2) / 0.2e1 + (m(3) * (t87 ^ 2 + t90 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t62 * mrSges(6,1) - t63 * mrSges(6,2) + Ifges(6,3) * t82 / 0.2e1) * t82 + (t78 * mrSges(5,2) - t66 * mrSges(5,3) + Ifges(5,1) * t77 / 0.2e1) * t77 + (-t78 * mrSges(5,1) + t67 * mrSges(5,3) + Ifges(5,4) * t77 + Ifges(5,2) * t76 / 0.2e1) * t76 + (t70 * mrSges(6,2) - t62 * mrSges(6,3) + Ifges(6,5) * t82 + Ifges(6,1) * t69 / 0.2e1) * t69 + (-t70 * mrSges(6,1) + t63 * mrSges(6,3) + Ifges(6,4) * t69 + Ifges(6,6) * t82 + Ifges(6,2) * t68 / 0.2e1) * t68 + (Ifges(3,3) * t94 + (t90 * mrSges(3,1) - t87 * mrSges(3,2)) * qJD(1) + (-t81 * mrSges(4,1) + (Ifges(4,2) * t94 + t95) * t89) * t89 + (Ifges(4,4) * qJD(2) * t89 + t81 * mrSges(4,2) + (Ifges(4,1) * t94 + t95) * t86) * t86) * qJD(2) + (t66 * mrSges(5,1) - t67 * mrSges(5,2) + Ifges(5,5) * t77 + Ifges(5,6) * t76 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(3) + (-t86 * mrSges(4,1) - t89 * mrSges(4,2)) * t80 + (Ifges(4,5) * t86 + Ifges(4,6) * t89) * qJD(2)) * qJD(3);
T = t1;
