% Calculate kinetic energy for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR7_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR7_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR7_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:00
% EndTime: 2019-12-31 16:36:01
% DurationCPUTime: 0.14s
% Computational Cost: add. (104->50), mult. (249->88), div. (0->0), fcn. (141->8), ass. (0->26)
t79 = sin(qJ(2));
t75 = sin(pkin(4));
t88 = qJD(1) * t75;
t70 = qJD(2) * pkin(6) + t79 * t88;
t78 = sin(qJ(3));
t81 = cos(qJ(3));
t76 = cos(pkin(4));
t87 = qJD(1) * t76;
t65 = t81 * t70 + t78 * t87;
t86 = qJD(2) * t78;
t85 = t81 * qJD(2);
t82 = cos(qJ(2));
t84 = t82 * t88;
t64 = -t78 * t70 + t81 * t87;
t80 = cos(qJ(4));
t77 = sin(qJ(4));
t73 = qJD(4) - t85;
t71 = -qJD(2) * pkin(2) - t84;
t69 = t77 * qJD(3) + t80 * t86;
t68 = t80 * qJD(3) - t77 * t86;
t66 = -t84 + (-pkin(3) * t81 - pkin(7) * t78 - pkin(2)) * qJD(2);
t63 = qJD(3) * pkin(7) + t65;
t62 = -qJD(3) * pkin(3) - t64;
t61 = t80 * t63 + t77 * t66;
t60 = -t77 * t63 + t80 * t66;
t1 = m(5) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + m(4) * (t64 ^ 2 + t65 ^ 2 + t71 ^ 2) / 0.2e1 + (t60 * mrSges(5,1) - t61 * mrSges(5,2) + Ifges(5,3) * t73 / 0.2e1) * t73 + (t64 * mrSges(4,1) - t65 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t62 * mrSges(5,2) - t60 * mrSges(5,3) + Ifges(5,5) * t73 + Ifges(5,1) * t69 / 0.2e1) * t69 + (m(3) * (t76 ^ 2 + (t79 ^ 2 + t82 ^ 2) * t75 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (-t62 * mrSges(5,1) + t61 * mrSges(5,3) + Ifges(5,4) * t69 + Ifges(5,6) * t73 + Ifges(5,2) * t68 / 0.2e1) * t68 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t82 - mrSges(3,2) * t79) * t88 + (-t71 * mrSges(4,1) + t65 * mrSges(4,3) + Ifges(4,6) * qJD(3) + Ifges(4,2) * t85 / 0.2e1) * t81 + (t71 * mrSges(4,2) - t64 * mrSges(4,3) + Ifges(4,5) * qJD(3) + (Ifges(4,4) * t81 + Ifges(4,1) * t78 / 0.2e1) * qJD(2)) * t78) * qJD(2);
T = t1;
