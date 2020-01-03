% Calculate kinetic energy for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR4_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR4_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:04
% EndTime: 2019-12-31 17:45:05
% DurationCPUTime: 0.16s
% Computational Cost: add. (129->52), mult. (270->86), div. (0->0), fcn. (124->6), ass. (0->24)
t85 = m(3) / 0.2e1;
t75 = sin(pkin(7));
t70 = (-pkin(1) * t75 - qJ(3)) * qJD(1);
t77 = cos(pkin(7));
t83 = -pkin(1) * t77 - pkin(2);
t67 = qJD(3) + (-qJ(4) + t83) * qJD(1);
t74 = sin(pkin(8));
t76 = cos(pkin(8));
t61 = t76 * qJD(2) + t74 * t67;
t84 = qJD(1) * t74;
t68 = qJD(4) - t70;
t60 = -t74 * qJD(2) + t76 * t67;
t80 = qJD(2) ^ 2;
t79 = cos(qJ(5));
t78 = sin(qJ(5));
t69 = t83 * qJD(1) + qJD(3);
t66 = (-t74 * t78 + t76 * t79) * qJD(1);
t65 = (-t74 * t79 - t76 * t78) * qJD(1);
t64 = pkin(4) * t84 + t68;
t59 = -pkin(6) * t84 + t61;
t58 = -t76 * qJD(1) * pkin(6) + t60;
t57 = t78 * t58 + t79 * t59;
t56 = t79 * t58 - t78 * t59;
t1 = m(5) * (t60 ^ 2 + t61 ^ 2 + t68 ^ 2) / 0.2e1 + m(6) * (t56 ^ 2 + t57 ^ 2 + t64 ^ 2) / 0.2e1 + t80 * t85 + m(4) * (t69 ^ 2 + t70 ^ 2 + t80) / 0.2e1 + (t64 * mrSges(6,2) - t56 * mrSges(6,3) + Ifges(6,1) * t66 / 0.2e1) * t66 + (-t64 * mrSges(6,1) + t57 * mrSges(6,3) + Ifges(6,4) * t66 + Ifges(6,2) * t65 / 0.2e1) * t65 + (t56 * mrSges(6,1) - t57 * mrSges(6,2) + Ifges(6,5) * t66 + Ifges(6,6) * t65 + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (t69 * mrSges(4,2) - t70 * mrSges(4,3) + t68 * (mrSges(5,1) * t74 + mrSges(5,2) * t76) + (-t60 * t76 - t61 * t74) * mrSges(5,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + (t77 * mrSges(3,1) - t75 * mrSges(3,2) + (t75 ^ 2 + t77 ^ 2) * t85 * pkin(1)) * pkin(1) + Ifges(5,1) * t76 ^ 2 / 0.2e1 + (-Ifges(5,4) * t76 + Ifges(5,2) * t74 / 0.2e1) * t74) * qJD(1)) * qJD(1);
T = t1;
