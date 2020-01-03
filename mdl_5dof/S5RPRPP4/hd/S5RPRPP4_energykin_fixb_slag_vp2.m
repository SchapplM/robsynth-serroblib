% Calculate kinetic energy for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
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
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPP4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP4_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP4_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:14
% EndTime: 2019-12-31 18:14:15
% DurationCPUTime: 0.21s
% Computational Cost: add. (179->66), mult. (371->96), div. (0->0), fcn. (180->4), ass. (0->23)
t63 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t75 = t63 * mrSges(4,3);
t74 = qJD(1) / 0.2e1;
t71 = cos(qJ(3));
t73 = -qJ(4) * qJD(1) + t63;
t57 = qJD(3) * pkin(3) + t73 * t71;
t70 = sin(qJ(3));
t58 = t73 * t70;
t68 = sin(pkin(7));
t69 = cos(pkin(7));
t53 = t68 * t57 + t69 * t58;
t67 = qJD(1) * qJ(2);
t61 = t70 * qJD(1) * pkin(3) + qJD(4) + t67;
t52 = t69 * t57 - t68 * t58;
t72 = qJD(1) ^ 2;
t66 = t72 * qJ(2) ^ 2;
t65 = -qJD(1) * pkin(1) + qJD(2);
t60 = (-t68 * t70 + t69 * t71) * qJD(1);
t59 = (t68 * t71 + t69 * t70) * qJD(1);
t54 = t59 * pkin(4) - t60 * qJ(5) + t61;
t51 = qJD(3) * qJ(5) + t53;
t50 = -qJD(3) * pkin(4) + qJD(5) - t52;
t1 = m(5) * (t52 ^ 2 + t53 ^ 2 + t61 ^ 2) / 0.2e1 + m(3) * (t65 ^ 2 + t66) / 0.2e1 + m(4) * (t66 + (t70 ^ 2 + t71 ^ 2) * t63 ^ 2) / 0.2e1 + m(6) * (t50 ^ 2 + t51 ^ 2 + t54 ^ 2) / 0.2e1 + (Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + qJ(2) * mrSges(3,3)) * t72 + (t61 * mrSges(5,2) + t50 * mrSges(6,2) - t52 * mrSges(5,3) - t54 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t60) * t60 + (t61 * mrSges(5,1) + t54 * mrSges(6,1) - t51 * mrSges(6,2) - t53 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t59 + (-Ifges(5,4) + Ifges(6,5)) * t60) * t59 + (t65 * mrSges(3,2) + (mrSges(4,2) * t67 + (Ifges(4,1) * t74 - t75) * t71) * t71 + ((qJ(2) * mrSges(4,1) - Ifges(4,4) * t71) * qJD(1) + (Ifges(4,2) * t74 - t75) * t70) * t70) * qJD(1) + (t52 * mrSges(5,1) - t50 * mrSges(6,1) - t53 * mrSges(5,2) + t51 * mrSges(6,3) + (Ifges(4,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(3) + (t71 * mrSges(4,1) - t70 * mrSges(4,2)) * t63 + (Ifges(6,4) + Ifges(5,5)) * t60 + (-Ifges(5,6) + Ifges(6,6)) * t59 + (Ifges(4,5) * t71 - Ifges(4,6) * t70) * qJD(1)) * qJD(3);
T = t1;
