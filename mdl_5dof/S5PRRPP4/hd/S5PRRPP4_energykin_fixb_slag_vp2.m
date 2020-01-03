% Calculate kinetic energy for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPP4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP4_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP4_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:43
% EndTime: 2019-12-31 17:40:44
% DurationCPUTime: 0.16s
% Computational Cost: add. (98->60), mult. (213->79), div. (0->0), fcn. (74->2), ass. (0->17)
t73 = pkin(3) + pkin(4);
t65 = sin(qJ(3));
t66 = cos(qJ(3));
t72 = qJD(2) * pkin(6);
t60 = t65 * qJD(1) + t66 * t72;
t71 = qJ(5) * qJD(2);
t58 = qJD(3) * qJ(4) + t60;
t70 = qJ(4) * t65 + pkin(2);
t59 = t66 * qJD(1) - t65 * t72;
t69 = qJD(4) - t59;
t67 = qJD(2) ^ 2;
t57 = (-pkin(3) * t66 - t70) * qJD(2);
t56 = -qJD(3) * pkin(3) + t69;
t55 = -t66 * t71 + t58;
t54 = qJD(5) + (t73 * t66 + t70) * qJD(2);
t53 = -t73 * qJD(3) - t65 * t71 + t69;
t1 = m(6) * (t53 ^ 2 + t54 ^ 2 + t55 ^ 2) / 0.2e1 + m(4) * (t67 * pkin(2) ^ 2 + t59 ^ 2 + t60 ^ 2) / 0.2e1 + m(5) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + t67 * Ifges(3,3) / 0.2e1 + (m(3) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + ((-t57 * mrSges(5,1) + t54 * mrSges(6,1) + t58 * mrSges(5,2) + t60 * mrSges(4,3) - t55 * mrSges(6,3) + (pkin(2) * mrSges(4,1) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t66) * qJD(2)) * t66 + (t56 * mrSges(5,2) + t54 * mrSges(6,2) - t59 * mrSges(4,3) - t57 * mrSges(5,3) - t53 * mrSges(6,3) + (-pkin(2) * mrSges(4,2) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t65 + (Ifges(4,4) - Ifges(6,4) - Ifges(5,5)) * t66) * qJD(2)) * t65) * qJD(2) + (t59 * mrSges(4,1) - t56 * mrSges(5,1) - t53 * mrSges(6,1) - t60 * mrSges(4,2) + t55 * mrSges(6,2) + t58 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(3) + ((Ifges(4,6) - Ifges(5,6) + Ifges(6,6)) * t66 + (Ifges(5,4) + Ifges(4,5) - Ifges(6,5)) * t65) * qJD(2)) * qJD(3);
T = t1;
