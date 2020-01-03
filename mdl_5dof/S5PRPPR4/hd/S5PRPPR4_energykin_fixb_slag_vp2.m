% Calculate kinetic energy for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPPR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR4_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR4_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR4_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:43
% EndTime: 2019-12-31 17:36:43
% DurationCPUTime: 0.13s
% Computational Cost: add. (104->53), mult. (244->81), div. (0->0), fcn. (126->4), ass. (0->21)
t75 = sin(pkin(8));
t76 = cos(pkin(8));
t81 = qJ(3) * qJD(2);
t70 = t75 * qJD(1) + t76 * t81;
t82 = qJD(2) * t76;
t80 = qJ(4) * t75 + pkin(2);
t69 = qJD(1) * t76 - t75 * t81;
t67 = qJD(4) - t69;
t78 = cos(qJ(5));
t77 = sin(qJ(5));
t73 = -qJD(2) * pkin(2) + qJD(3);
t68 = t70 ^ 2;
t66 = (t75 * t78 - t76 * t77) * qJD(2);
t65 = (-t75 * t77 - t76 * t78) * qJD(2);
t64 = qJD(3) + (-pkin(3) * t76 - t80) * qJD(2);
t63 = -pkin(6) * t82 + t70;
t62 = -pkin(6) * qJD(2) * t75 + t67;
t61 = -qJD(3) + ((pkin(3) + pkin(4)) * t76 + t80) * qJD(2);
t60 = t62 * t77 + t63 * t78;
t59 = t62 * t78 - t63 * t77;
t1 = m(4) * (t69 ^ 2 + t73 ^ 2 + t68) / 0.2e1 + m(5) * (t64 ^ 2 + t67 ^ 2 + t68) / 0.2e1 + m(6) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + (m(3) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t61 * mrSges(6,2) - t59 * mrSges(6,3) + Ifges(6,1) * t66 / 0.2e1) * t66 + (-t61 * mrSges(6,1) + t60 * mrSges(6,3) + Ifges(6,4) * t66 + Ifges(6,2) * t65 / 0.2e1) * t65 + (t59 * mrSges(6,1) - t60 * mrSges(6,2) + Ifges(6,5) * t66 + Ifges(6,6) * t65 + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (Ifges(3,3) * qJD(2) / 0.2e1 + (-t73 * mrSges(4,1) - t64 * mrSges(5,1) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t82 + (mrSges(5,2) + mrSges(4,3)) * t70) * t76 + (t73 * mrSges(4,2) + t67 * mrSges(5,2) - t69 * mrSges(4,3) - t64 * mrSges(5,3) + ((Ifges(4,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t75 + (Ifges(4,4) - Ifges(5,5)) * t76) * qJD(2)) * t75) * qJD(2);
T = t1;
