% Calculate kinetic energy for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR8_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR8_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR8_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:02
% EndTime: 2019-12-31 18:01:02
% DurationCPUTime: 0.11s
% Computational Cost: add. (160->46), mult. (248->70), div. (0->0), fcn. (88->6), ass. (0->24)
t81 = m(3) / 0.2e1;
t70 = -qJD(1) + qJD(4);
t80 = t70 / 0.2e1;
t68 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t72 = cos(pkin(8));
t67 = t72 * t68;
t71 = sin(pkin(8));
t64 = t67 + (-qJ(2) * t71 - pkin(3)) * qJD(1);
t79 = qJ(2) * qJD(1);
t66 = t71 * t68 + t72 * t79;
t74 = sin(qJ(4));
t76 = cos(qJ(4));
t61 = t74 * t64 + t76 * t66;
t60 = t76 * t64 - t74 * t66;
t77 = qJD(3) ^ 2;
t75 = cos(qJ(5));
t73 = sin(qJ(5));
t69 = -qJD(1) * pkin(1) + qJD(2);
t65 = -t71 * t79 + t67;
t59 = t70 * pkin(7) + t61;
t58 = -t70 * pkin(4) - t60;
t57 = t73 * qJD(3) + t75 * t59;
t56 = t75 * qJD(3) - t73 * t59;
t1 = m(6) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + t69 ^ 2 * t81 + m(4) * (t65 ^ 2 + t66 ^ 2 + t77) / 0.2e1 + m(5) * (t60 ^ 2 + t61 ^ 2 + t77) / 0.2e1 + (t56 * mrSges(6,1) - t57 * mrSges(6,2) + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (-t61 * mrSges(5,2) + t60 * mrSges(5,1) + Ifges(5,3) * t80 + (Ifges(6,2) * t75 * t80 - t58 * mrSges(6,1) + t57 * mrSges(6,3) + Ifges(6,6) * qJD(5)) * t75 + (t58 * mrSges(6,2) - t56 * mrSges(6,3) + Ifges(6,5) * qJD(5) + (Ifges(6,4) * t75 + Ifges(6,1) * t73 / 0.2e1) * t70) * t73) * t70 + (-t69 * mrSges(3,1) - t65 * mrSges(4,1) + t66 * mrSges(4,2) + (Ifges(2,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + (qJ(2) * t81 + mrSges(3,3)) * qJ(2)) * qJD(1)) * qJD(1);
T = t1;
