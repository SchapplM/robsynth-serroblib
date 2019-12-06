% Calculate kinetic energy for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRP3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP3_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP3_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:27
% EndTime: 2019-12-05 16:43:28
% DurationCPUTime: 0.21s
% Computational Cost: add. (166->64), mult. (367->96), div. (0->0), fcn. (208->4), ass. (0->22)
t78 = cos(qJ(3));
t73 = t78 * qJD(1);
t76 = sin(qJ(3));
t82 = qJD(2) * t76;
t64 = qJD(3) * pkin(3) + t73 + (-pkin(7) - pkin(6)) * t82;
t81 = qJD(2) * t78;
t69 = pkin(6) * t81 + t76 * qJD(1);
t65 = pkin(7) * t81 + t69;
t75 = sin(qJ(4));
t77 = cos(qJ(4));
t59 = t75 * t64 + t77 * t65;
t58 = t77 * t64 - t75 * t65;
t70 = (-pkin(3) * t78 - pkin(2)) * qJD(2);
t79 = qJD(2) ^ 2;
t74 = qJD(3) + qJD(4);
t68 = -pkin(6) * t82 + t73;
t67 = (t75 * t78 + t76 * t77) * qJD(2);
t66 = (-t75 * t76 + t77 * t78) * qJD(2);
t60 = -t66 * pkin(4) + qJD(5) + t70;
t57 = t66 * qJ(5) + t59;
t56 = t74 * pkin(4) - t67 * qJ(5) + t58;
t1 = m(5) * (t58 ^ 2 + t59 ^ 2 + t70 ^ 2) / 0.2e1 + m(6) * (t56 ^ 2 + t57 ^ 2 + t60 ^ 2) / 0.2e1 + m(4) * (t79 * pkin(2) ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + t79 * Ifges(3,3) / 0.2e1 + (m(3) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t68 * mrSges(4,1) - t69 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + ((t69 * mrSges(4,3) + Ifges(4,6) * qJD(3) + (pkin(2) * mrSges(4,1) + Ifges(4,2) * t78 / 0.2e1) * qJD(2)) * t78 + (-t68 * mrSges(4,3) + Ifges(4,5) * qJD(3) + (-pkin(2) * mrSges(4,2) + Ifges(4,4) * t78 + Ifges(4,1) * t76 / 0.2e1) * qJD(2)) * t76) * qJD(2) + (t58 * mrSges(5,1) + t56 * mrSges(6,1) - t59 * mrSges(5,2) - t57 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t74) * t74 + (t70 * mrSges(5,2) + t60 * mrSges(6,2) - t58 * mrSges(5,3) - t56 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t67 + (Ifges(5,5) + Ifges(6,5)) * t74) * t67 + (-t70 * mrSges(5,1) - t60 * mrSges(6,1) + t59 * mrSges(5,3) + t57 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t66 + (Ifges(5,6) + Ifges(6,6)) * t74 + (Ifges(5,4) + Ifges(6,4)) * t67) * t66;
T = t1;
