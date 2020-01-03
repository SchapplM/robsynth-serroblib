% Calculate kinetic energy for
% S5RRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP2_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP2_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP2_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:26
% EndTime: 2019-12-31 20:51:26
% DurationCPUTime: 0.18s
% Computational Cost: add. (178->66), mult. (255->89), div. (0->0), fcn. (84->4), ass. (0->21)
t79 = pkin(3) + pkin(4);
t65 = qJD(1) + qJD(2);
t68 = sin(qJ(2));
t77 = pkin(1) * qJD(1);
t63 = t65 * pkin(7) + t68 * t77;
t78 = t63 * mrSges(4,3);
t69 = cos(qJ(3));
t59 = qJD(3) * qJ(4) + t69 * t63;
t76 = qJ(5) * t65;
t67 = sin(qJ(3));
t75 = t67 * t63 + qJD(4);
t70 = cos(qJ(2));
t74 = t70 * t77;
t73 = qJ(4) * t67 + pkin(2);
t64 = -t65 * pkin(2) - t74;
t58 = -qJD(3) * pkin(3) + t75;
t57 = -t69 * t76 + t59;
t56 = -t74 + (-pkin(3) * t69 - t73) * t65;
t55 = -t79 * qJD(3) - t67 * t76 + t75;
t54 = t74 + qJD(5) + (t79 * t69 + t73) * t65;
t1 = m(6) * (t54 ^ 2 + t55 ^ 2 + t57 ^ 2) / 0.2e1 + m(5) * (t56 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + m(4) * (t64 ^ 2 + (t67 ^ 2 + t69 ^ 2) * t63 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t68 ^ 2 + t70 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (-t58 * mrSges(5,1) - t55 * mrSges(6,1) + t57 * mrSges(6,2) + t59 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * qJD(3) + (-t67 * mrSges(4,1) - t69 * mrSges(4,2)) * t63) * qJD(3) + (Ifges(3,3) * t65 / 0.2e1 + (mrSges(3,1) * t70 - mrSges(3,2) * t68) * t77 + (-t64 * mrSges(4,1) - t56 * mrSges(5,1) + t54 * mrSges(6,1) + t59 * mrSges(5,2) - t57 * mrSges(6,3) + (t78 + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t65) * t69 + (Ifges(4,6) - Ifges(5,6) + Ifges(6,6)) * qJD(3)) * t69 + (-t55 * mrSges(6,3) + t58 * mrSges(5,2) + t54 * mrSges(6,2) - t56 * mrSges(5,3) + t64 * mrSges(4,2) + (t78 + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t65) * t67 + (Ifges(4,4) - Ifges(6,4) - Ifges(5,5)) * t65 * t69 + (Ifges(5,4) + Ifges(4,5) - Ifges(6,5)) * qJD(3)) * t67) * t65;
T = t1;
