% Calculate kinetic energy for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR9_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR9_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR9_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:38
% EndTime: 2019-12-31 18:23:38
% DurationCPUTime: 0.28s
% Computational Cost: add. (160->70), mult. (340->104), div. (0->0), fcn. (150->6), ass. (0->29)
t88 = m(3) / 0.2e1;
t87 = -pkin(3) - pkin(7);
t74 = sin(pkin(8));
t70 = (pkin(1) * t74 + pkin(6)) * qJD(1);
t77 = sin(qJ(3));
t79 = cos(qJ(3));
t65 = t77 * qJD(2) + t79 * t70;
t86 = qJD(1) * t79;
t85 = t77 * qJD(1);
t75 = cos(pkin(8));
t84 = -pkin(1) * t75 - pkin(2);
t64 = t79 * qJD(2) - t77 * t70;
t62 = -qJD(3) * qJ(4) - t65;
t83 = qJD(4) - t64;
t82 = -qJ(4) * t77 + t84;
t78 = cos(qJ(5));
t76 = sin(qJ(5));
t72 = qJD(5) + t85;
t71 = t84 * qJD(1);
t69 = t78 * qJD(3) - t76 * t86;
t68 = -t76 * qJD(3) - t78 * t86;
t63 = (-pkin(3) * t79 + t82) * qJD(1);
t61 = -qJD(3) * pkin(3) + t83;
t60 = (t79 * t87 + t82) * qJD(1);
t59 = pkin(4) * t86 - t62;
t58 = pkin(4) * t85 + qJD(3) * t87 + t83;
t57 = t76 * t58 + t78 * t60;
t56 = t78 * t58 - t76 * t60;
t1 = qJD(2) ^ 2 * t88 + m(4) * (t64 ^ 2 + t65 ^ 2 + t71 ^ 2) / 0.2e1 + m(6) * (t56 ^ 2 + t57 ^ 2 + t59 ^ 2) / 0.2e1 + m(5) * (t61 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + (t56 * mrSges(6,1) - t57 * mrSges(6,2) + Ifges(6,3) * t72 / 0.2e1) * t72 + (t59 * mrSges(6,2) - t56 * mrSges(6,3) + Ifges(6,5) * t72 + Ifges(6,1) * t69 / 0.2e1) * t69 + (-t59 * mrSges(6,1) + t57 * mrSges(6,3) + Ifges(6,4) * t69 + Ifges(6,6) * t72 + Ifges(6,2) * t68 / 0.2e1) * t68 + (t64 * mrSges(4,1) - t65 * mrSges(4,2) + t61 * mrSges(5,2) - t62 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(3)) * qJD(3) + ((-t71 * mrSges(4,1) - t62 * mrSges(5,1) + t63 * mrSges(5,2) + t65 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t86) * t79 + (t61 * mrSges(5,1) + t71 * mrSges(4,2) - t64 * mrSges(4,3) - t63 * mrSges(5,3)) * t77 + ((-Ifges(5,5) + Ifges(4,6)) * t79 + (-Ifges(5,4) + Ifges(4,5)) * t77) * qJD(3) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t75 * mrSges(3,1) - t74 * mrSges(3,2) + (t74 ^ 2 + t75 ^ 2) * t88 * pkin(1)) * pkin(1) + ((Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t77 + (Ifges(4,4) + Ifges(5,6)) * t79) * t77) * qJD(1)) * qJD(1);
T = t1;
