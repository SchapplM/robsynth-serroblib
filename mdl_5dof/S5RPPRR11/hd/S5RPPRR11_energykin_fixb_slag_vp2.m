% Calculate kinetic energy for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR11_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR11_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR11_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:29
% EndTime: 2019-12-31 18:05:30
% DurationCPUTime: 0.16s
% Computational Cost: add. (122->55), mult. (222->86), div. (0->0), fcn. (80->4), ass. (0->23)
t72 = -pkin(1) - qJ(3);
t60 = -t72 * qJD(1) - qJD(2);
t76 = t60 ^ 2;
t75 = m(3) / 0.2e1;
t63 = qJD(1) * qJ(2) + qJD(3);
t59 = -qJD(1) * pkin(6) + t63;
t74 = t59 * mrSges(5,3);
t73 = qJD(1) / 0.2e1;
t69 = cos(qJ(4));
t71 = qJD(1) * t69;
t68 = cos(qJ(5));
t67 = sin(qJ(4));
t66 = sin(qJ(5));
t64 = -qJD(1) * pkin(1) + qJD(2);
t62 = t67 * qJD(1) + qJD(5);
t57 = t66 * qJD(4) + t68 * t71;
t56 = t68 * qJD(4) - t66 * t71;
t55 = -qJD(4) * pkin(4) - t69 * t59;
t54 = qJD(4) * pkin(7) + t67 * t59;
t53 = -qJD(2) + (pkin(4) * t67 - pkin(7) * t69 - t72) * qJD(1);
t52 = t66 * t53 + t68 * t54;
t51 = t68 * t53 - t66 * t54;
t1 = t64 ^ 2 * t75 + m(6) * (t51 ^ 2 + t52 ^ 2 + t55 ^ 2) / 0.2e1 + m(4) * (t63 ^ 2 + t76) / 0.2e1 + m(5) * (t76 + (t67 ^ 2 + t69 ^ 2) * t59 ^ 2) / 0.2e1 + (t51 * mrSges(6,1) - t52 * mrSges(6,2) + Ifges(6,3) * t62 / 0.2e1) * t62 + (Ifges(5,3) * qJD(4) / 0.2e1 + (t69 * mrSges(5,1) - t67 * mrSges(5,2)) * t59) * qJD(4) + (t55 * mrSges(6,2) - t51 * mrSges(6,3) + Ifges(6,5) * t62 + Ifges(6,1) * t57 / 0.2e1) * t57 + (-t55 * mrSges(6,1) + t52 * mrSges(6,3) + Ifges(6,4) * t57 + Ifges(6,6) * t62 + Ifges(6,2) * t56 / 0.2e1) * t56 + (t64 * mrSges(3,2) + t63 * mrSges(4,2) + t60 * mrSges(4,3) + (t60 * mrSges(5,2) + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t73 - t74) * t69) * t69 + (-Ifges(5,4) * t71 + t60 * mrSges(5,1) - Ifges(5,6) * qJD(4) + (Ifges(5,2) * t73 - t74) * t67) * t67 + (Ifges(4,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + (qJ(2) * t75 + mrSges(3,3)) * qJ(2)) * qJD(1)) * qJD(1);
T = t1;
