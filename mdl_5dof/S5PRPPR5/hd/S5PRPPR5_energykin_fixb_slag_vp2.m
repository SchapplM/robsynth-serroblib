% Calculate kinetic energy for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPPR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR5_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR5_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:37:58
% EndTime: 2019-12-31 17:37:58
% DurationCPUTime: 0.10s
% Computational Cost: add. (100->42), mult. (176->68), div. (0->0), fcn. (68->6), ass. (0->17)
t74 = cos(qJ(2));
t77 = -t74 * qJD(1) + qJD(3);
t65 = (-pkin(2) - pkin(3)) * qJD(2) + t77;
t72 = sin(qJ(2));
t68 = qJD(2) * qJ(3) + t72 * qJD(1);
t69 = sin(pkin(8));
t70 = cos(pkin(8));
t63 = t69 * t65 + t70 * t68;
t62 = t70 * t65 - t69 * t68;
t73 = cos(qJ(5));
t71 = sin(qJ(5));
t67 = -qJD(2) * pkin(2) + t77;
t61 = -qJD(2) * pkin(6) + t63;
t60 = qJD(2) * pkin(4) - t62;
t59 = t71 * qJD(4) + t73 * t61;
t58 = t73 * qJD(4) - t71 * t61;
t1 = m(6) * (t58 ^ 2 + t59 ^ 2 + t60 ^ 2) / 0.2e1 + m(4) * (t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(5) * (qJD(4) ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + (m(3) * (t72 ^ 2 + t74 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t58 * mrSges(6,1) - t59 * mrSges(6,2) + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (-t67 * mrSges(4,1) - t62 * mrSges(5,1) + t63 * mrSges(5,2) + t68 * mrSges(4,3) + (t74 * mrSges(3,1) - t72 * mrSges(3,2)) * qJD(1) + (t60 * mrSges(6,1) - t59 * mrSges(6,3) - Ifges(6,6) * qJD(5)) * t73 + (-t60 * mrSges(6,2) + t58 * mrSges(6,3) - Ifges(6,5) * qJD(5)) * t71 + (Ifges(3,3) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + t73 ^ 2 * Ifges(6,2) / 0.2e1 + (Ifges(6,4) * t73 + Ifges(6,1) * t71 / 0.2e1) * t71) * qJD(2)) * qJD(2);
T = t1;
