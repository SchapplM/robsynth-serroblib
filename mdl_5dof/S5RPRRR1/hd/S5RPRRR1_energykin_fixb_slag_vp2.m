% Calculate kinetic energy for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(1,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_energykin_fixb_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR1_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR1_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:24:50
% EndTime: 2019-07-18 13:24:51
% DurationCPUTime: 0.27s
% Computational Cost: add. (139->60), mult. (322->100), div. (0->0), fcn. (192->6), ass. (0->28)
t57 = sin(qJ(4));
t60 = cos(qJ(4));
t61 = cos(qJ(3));
t67 = qJD(1) * qJ(2);
t65 = t61 * t67;
t48 = -t60 * qJD(2) + t57 * t65;
t73 = t48 ^ 2;
t63 = qJD(1) ^ 2;
t71 = t63 * qJ(2) ^ 2;
t70 = qJ(2) * mrSges(4,3);
t58 = sin(qJ(3));
t69 = t58 * qJD(1);
t68 = t61 * qJD(1);
t66 = t58 * t67;
t51 = t60 * qJD(3) - t57 * t69;
t62 = qJD(2) ^ 2;
t59 = cos(qJ(5));
t56 = sin(qJ(5));
t54 = t58 ^ 2 * t71;
t53 = qJD(4) - t68;
t52 = t57 * qJD(3) + t60 * t69;
t50 = t57 * qJD(2) + t60 * t65;
t47 = qJD(5) - t51;
t46 = t59 * t50 + t56 * t66;
t45 = -t56 * t50 + t59 * t66;
t44 = t59 * t52 + t56 * t53;
t43 = -t56 * t52 + t59 * t53;
t1 = m(3) * (t62 + t71) / 0.2e1 + m(4) * (t61 ^ 2 * t71 + t54 + t62) / 0.2e1 + m(5) * (t50 ^ 2 + t54 + t73) / 0.2e1 + Ifges(4,3) * qJD(3) ^ 2 / 0.2e1 + m(6) * (t45 ^ 2 + t46 ^ 2 + t73) / 0.2e1 + t63 * qJ(2) * mrSges(3,3) + (t48 * mrSges(5,3) + Ifges(5,1) * t52 / 0.2e1) * t52 + (t45 * mrSges(6,1) - t46 * mrSges(6,2) + Ifges(6,3) * t47 / 0.2e1) * t47 + (-t48 * mrSges(5,1) - t50 * mrSges(5,2) + Ifges(5,5) * t52 + Ifges(5,3) * t53 / 0.2e1) * t53 + (t48 * mrSges(6,2) - t45 * mrSges(6,3) + Ifges(6,5) * t47 + Ifges(6,1) * t44 / 0.2e1) * t44 + (t50 * mrSges(5,3) + Ifges(5,4) * t52 + Ifges(5,6) * t53 + Ifges(5,2) * t51 / 0.2e1) * t51 + (-t48 * mrSges(6,1) + t46 * mrSges(6,3) + Ifges(6,4) * t44 + Ifges(6,6) * t47 + Ifges(6,2) * t43 / 0.2e1) * t43 + (-qJD(2) * mrSges(3,1) + (-qJD(2) * mrSges(4,1) + (Ifges(4,2) / 0.2e1 + t70) * t68 + (-qJ(2) * mrSges(4,2) + Ifges(4,6)) * qJD(3)) * t61 + (qJD(2) * mrSges(4,2) + Ifges(4,5) * qJD(3) + (-qJD(3) * mrSges(4,1) - t51 * mrSges(5,1) + t52 * mrSges(5,2)) * qJ(2) + (Ifges(4,4) * t61 + (Ifges(4,1) / 0.2e1 + t70) * t58) * qJD(1)) * t58) * qJD(1) + (Ifges(3,2) + Ifges(2,3)) * t63 / 0.2e1;
T  = t1;
