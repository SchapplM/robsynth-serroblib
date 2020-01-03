% Calculate kinetic energy for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR7_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR7_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:28
% EndTime: 2019-12-31 17:59:28
% DurationCPUTime: 0.23s
% Computational Cost: add. (136->58), mult. (271->91), div. (0->0), fcn. (114->6), ass. (0->26)
t72 = sin(pkin(8));
t81 = -pkin(1) * t72 - qJ(3);
t68 = t81 * qJD(1);
t85 = t68 ^ 2;
t84 = m(3) / 0.2e1;
t73 = cos(pkin(8));
t82 = -pkin(1) * t73 - pkin(2);
t64 = qJD(3) + (-pkin(6) + t82) * qJD(1);
t75 = sin(qJ(4));
t77 = cos(qJ(4));
t61 = t77 * qJD(2) + t75 * t64;
t83 = qJD(1) * t77;
t60 = -t75 * qJD(2) + t77 * t64;
t78 = qJD(2) ^ 2;
t76 = cos(qJ(5));
t74 = sin(qJ(5));
t70 = t75 * qJD(1) + qJD(5);
t67 = t82 * qJD(1) + qJD(3);
t66 = t74 * qJD(4) + t76 * t83;
t65 = t76 * qJD(4) - t74 * t83;
t62 = (pkin(4) * t75 - pkin(7) * t77 - t81) * qJD(1);
t59 = qJD(4) * pkin(7) + t61;
t58 = -qJD(4) * pkin(4) - t60;
t57 = t76 * t59 + t74 * t62;
t56 = -t74 * t59 + t76 * t62;
t1 = m(6) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + m(5) * (t60 ^ 2 + t61 ^ 2 + t85) / 0.2e1 + t78 * t84 + m(4) * (t67 ^ 2 + t78 + t85) / 0.2e1 + (t56 * mrSges(6,1) - t57 * mrSges(6,2) + Ifges(6,3) * t70 / 0.2e1) * t70 + (t60 * mrSges(5,1) - t61 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t58 * mrSges(6,2) - t56 * mrSges(6,3) + Ifges(6,5) * t70 + Ifges(6,1) * t66 / 0.2e1) * t66 + (-t58 * mrSges(6,1) + t57 * mrSges(6,3) + Ifges(6,4) * t66 + Ifges(6,6) * t70 + Ifges(6,2) * t65 / 0.2e1) * t65 + (t67 * mrSges(4,2) - t68 * mrSges(4,3) + (-t68 * mrSges(5,2) - t60 * mrSges(5,3) + Ifges(5,5) * qJD(4) + Ifges(5,1) * t83 / 0.2e1) * t77 + (-t68 * mrSges(5,1) - t61 * mrSges(5,3) - Ifges(5,6) * qJD(4)) * t75 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + (t73 * mrSges(3,1) - t72 * mrSges(3,2) + (t72 ^ 2 + t73 ^ 2) * t84 * pkin(1)) * pkin(1) + (-Ifges(5,4) * t77 + Ifges(5,2) * t75 / 0.2e1) * t75) * qJD(1)) * qJD(1);
T = t1;
