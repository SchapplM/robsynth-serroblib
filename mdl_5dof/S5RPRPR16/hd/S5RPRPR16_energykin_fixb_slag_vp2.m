% Calculate kinetic energy for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR16_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR16_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR16_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:40
% EndTime: 2019-12-31 18:38:40
% DurationCPUTime: 0.21s
% Computational Cost: add. (153->70), mult. (291->101), div. (0->0), fcn. (110->4), ass. (0->27)
t63 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t79 = t63 * mrSges(4,3);
t68 = qJD(1) * qJ(2);
t70 = sin(qJ(3));
t77 = qJD(1) * t70;
t78 = pkin(3) * t77 + t68;
t72 = cos(qJ(3));
t76 = t72 * qJD(1);
t75 = qJD(3) * qJ(4);
t74 = pkin(4) * qJD(1) - t63;
t73 = qJD(1) ^ 2;
t71 = cos(qJ(5));
t69 = sin(qJ(5));
t67 = t73 * qJ(2) ^ 2;
t66 = -pkin(1) * qJD(1) + qJD(2);
t64 = qJD(5) + t76;
t61 = qJD(3) * t71 + t69 * t77;
t60 = -qJD(3) * t69 + t71 * t77;
t59 = -t63 * t70 - t75;
t58 = -qJ(4) * t76 + t78;
t57 = -qJD(3) * pkin(3) - t72 * t63 + qJD(4);
t56 = -t70 * t74 + t75;
t55 = (pkin(7) * t70 - qJ(4) * t72) * qJD(1) + t78;
t54 = qJD(4) + t74 * t72 + (-pkin(3) - pkin(7)) * qJD(3);
t53 = t54 * t69 + t55 * t71;
t52 = t54 * t71 - t55 * t69;
t1 = m(3) * (t66 ^ 2 + t67) / 0.2e1 + m(4) * (t67 + (t70 ^ 2 + t72 ^ 2) * t63 ^ 2) / 0.2e1 + m(5) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + m(6) * (t52 ^ 2 + t53 ^ 2 + t56 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + qJ(2) * mrSges(3,3)) * t73 + (t52 * mrSges(6,1) - t53 * mrSges(6,2) + Ifges(6,3) * t64 / 0.2e1) * t64 + (t56 * mrSges(6,2) - t52 * mrSges(6,3) + Ifges(6,5) * t64 + Ifges(6,1) * t61 / 0.2e1) * t61 + (-t56 * mrSges(6,1) + t53 * mrSges(6,3) + Ifges(6,4) * t61 + Ifges(6,6) * t64 + Ifges(6,2) * t60 / 0.2e1) * t60 + (t57 * mrSges(5,2) - t59 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * qJD(3) + (mrSges(4,1) * t72 - mrSges(4,2) * t70) * t63) * qJD(3) + (t66 * mrSges(3,2) + (mrSges(4,2) * t68 + t57 * mrSges(5,1) - t58 * mrSges(5,3) + (-t79 + (Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1) * qJD(1)) * t72 + (-Ifges(5,4) + Ifges(4,5)) * qJD(3)) * t72 + (t59 * mrSges(5,1) - t58 * mrSges(5,2) + (qJ(2) * mrSges(4,1) + (-Ifges(4,4) - Ifges(5,6)) * t72) * qJD(1) + (-t79 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(1)) * t70 + (Ifges(5,5) - Ifges(4,6)) * qJD(3)) * t70) * qJD(1);
T = t1;
