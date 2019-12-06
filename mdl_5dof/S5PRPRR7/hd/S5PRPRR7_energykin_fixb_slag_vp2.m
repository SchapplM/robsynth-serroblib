% Calculate kinetic energy for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR7_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR7_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:39
% EndTime: 2019-12-05 15:59:39
% DurationCPUTime: 0.19s
% Computational Cost: add. (128->54), mult. (247->90), div. (0->0), fcn. (114->6), ass. (0->24)
t72 = sin(qJ(2));
t80 = t72 * qJD(1);
t67 = qJD(2) * qJ(3) + t80;
t83 = t67 ^ 2;
t75 = cos(qJ(2));
t78 = -t75 * qJD(1) + qJD(3);
t65 = (-pkin(2) - pkin(6)) * qJD(2) + t78;
t82 = t65 * mrSges(5,3);
t81 = qJD(2) / 0.2e1;
t79 = -pkin(7) * qJD(2) + t65;
t74 = cos(qJ(4));
t73 = cos(qJ(5));
t71 = sin(qJ(4));
t70 = sin(qJ(5));
t69 = qJD(4) + qJD(5);
t66 = -qJD(2) * pkin(2) + t78;
t63 = t80 + (pkin(4) * t71 + qJ(3)) * qJD(2);
t62 = (-t70 * t71 + t73 * t74) * qJD(2);
t61 = (-t70 * t74 - t71 * t73) * qJD(2);
t60 = t79 * t71;
t59 = qJD(4) * pkin(4) + t79 * t74;
t58 = t70 * t59 + t73 * t60;
t57 = t73 * t59 - t70 * t60;
t1 = m(5) * (t83 + (t71 ^ 2 + t74 ^ 2) * t65 ^ 2) / 0.2e1 + m(6) * (t57 ^ 2 + t58 ^ 2 + t63 ^ 2) / 0.2e1 + m(4) * (t66 ^ 2 + t83) / 0.2e1 + (m(3) * (t72 ^ 2 + t75 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t57 * mrSges(6,1) - t58 * mrSges(6,2) + Ifges(6,3) * t69 / 0.2e1) * t69 + (Ifges(5,3) * qJD(4) / 0.2e1 + (t74 * mrSges(5,1) - t71 * mrSges(5,2)) * t65) * qJD(4) + (t63 * mrSges(6,2) - t57 * mrSges(6,3) + Ifges(6,5) * t69 + Ifges(6,1) * t62 / 0.2e1) * t62 + (-t63 * mrSges(6,1) + t58 * mrSges(6,3) + Ifges(6,4) * t62 + Ifges(6,6) * t69 + Ifges(6,2) * t61 / 0.2e1) * t61 + (t66 * mrSges(4,2) + t67 * mrSges(4,3) + (t75 * mrSges(3,1) - t72 * mrSges(3,2)) * qJD(1) + (t67 * mrSges(5,2) + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t81 - t82) * t74) * t74 + (t67 * mrSges(5,1) - Ifges(5,6) * qJD(4) + (Ifges(5,2) * t81 - t82) * t71) * t71 + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1 - Ifges(5,4) * t74 * t71) * qJD(2)) * qJD(2);
T = t1;
