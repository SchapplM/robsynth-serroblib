% Calculate kinetic energy for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR5_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR5_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR5_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:36
% EndTime: 2019-12-31 17:35:36
% DurationCPUTime: 0.08s
% Computational Cost: add. (83->37), mult. (142->63), div. (0->0), fcn. (64->6), ass. (0->19)
t65 = qJD(3) + qJD(4);
t75 = t65 / 0.2e1;
t71 = cos(qJ(3));
t63 = qJD(3) * pkin(3) + t71 * qJD(2);
t67 = sin(qJ(4));
t70 = cos(qJ(4));
t68 = sin(qJ(3));
t74 = qJD(2) * t68;
t61 = t67 * t63 + t70 * t74;
t60 = t70 * t63 - t67 * t74;
t73 = qJD(1) ^ 2;
t72 = qJD(2) ^ 2;
t69 = cos(qJ(5));
t66 = sin(qJ(5));
t59 = t65 * pkin(7) + t61;
t58 = -t65 * pkin(4) - t60;
t57 = -t66 * qJD(1) + t69 * t59;
t56 = -t69 * qJD(1) - t66 * t59;
t1 = m(2) * t73 / 0.2e1 + m(6) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + m(5) * (t60 ^ 2 + t61 ^ 2 + t73) / 0.2e1 + m(3) * (t72 + t73) / 0.2e1 + m(4) * (t73 + (t68 ^ 2 + t71 ^ 2) * t72) / 0.2e1 + (t56 * mrSges(6,1) - t57 * mrSges(6,2) + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (Ifges(4,3) * qJD(3) / 0.2e1 + (t71 * mrSges(4,1) - t68 * mrSges(4,2)) * qJD(2)) * qJD(3) + (-t61 * mrSges(5,2) + t60 * mrSges(5,1) + Ifges(5,3) * t75 + (Ifges(6,2) * t69 * t75 - t58 * mrSges(6,1) + t57 * mrSges(6,3) + Ifges(6,6) * qJD(5)) * t69 + (t58 * mrSges(6,2) - t56 * mrSges(6,3) + Ifges(6,5) * qJD(5) + (Ifges(6,4) * t69 + Ifges(6,1) * t66 / 0.2e1) * t65) * t66) * t65;
T = t1;
