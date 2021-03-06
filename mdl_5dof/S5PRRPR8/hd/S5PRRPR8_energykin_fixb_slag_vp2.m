% Calculate kinetic energy for
% S5PRRPR8
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR8_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR8_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:15
% EndTime: 2019-12-31 17:42:15
% DurationCPUTime: 0.12s
% Computational Cost: add. (146->42), mult. (248->73), div. (0->0), fcn. (140->8), ass. (0->21)
t83 = cos(qJ(2));
t73 = qJD(2) * pkin(2) + qJD(1) * t83;
t79 = sin(qJ(3));
t82 = cos(qJ(3));
t80 = sin(qJ(2));
t85 = qJD(1) * t80;
t70 = t82 * t73 - t79 * t85;
t75 = qJD(2) + qJD(3);
t68 = pkin(3) * t75 + t70;
t71 = t73 * t79 + t82 * t85;
t76 = sin(pkin(9));
t77 = cos(pkin(9));
t66 = t76 * t68 + t77 * t71;
t65 = t68 * t77 - t71 * t76;
t81 = cos(qJ(5));
t78 = sin(qJ(5));
t64 = pkin(7) * t75 + t66;
t63 = -pkin(4) * t75 - t65;
t62 = qJD(4) * t78 + t64 * t81;
t61 = qJD(4) * t81 - t64 * t78;
t1 = m(4) * (t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (qJD(4) ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + m(6) * (t61 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + (m(2) / 0.2e1 + m(3) * (t80 ^ 2 + t83 ^ 2) / 0.2e1) * qJD(1) ^ 2 + (t61 * mrSges(6,1) - t62 * mrSges(6,2) + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (Ifges(3,3) * qJD(2) / 0.2e1 + (t83 * mrSges(3,1) - t80 * mrSges(3,2)) * qJD(1)) * qJD(2) + (t70 * mrSges(4,1) + t65 * mrSges(5,1) - t71 * mrSges(4,2) - t66 * mrSges(5,2) + (-t63 * mrSges(6,1) + t62 * mrSges(6,3) + Ifges(6,6) * qJD(5)) * t81 + (t63 * mrSges(6,2) - t61 * mrSges(6,3) + Ifges(6,5) * qJD(5)) * t78 + (Ifges(5,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(6,2) * t81 ^ 2 / 0.2e1 + (Ifges(6,4) * t81 + Ifges(6,1) * t78 / 0.2e1) * t78) * t75) * t75;
T = t1;
