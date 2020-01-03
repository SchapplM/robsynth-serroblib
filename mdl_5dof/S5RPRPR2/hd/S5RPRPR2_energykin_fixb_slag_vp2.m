% Calculate kinetic energy for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR2_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR2_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:33:35
% EndTime: 2020-01-03 11:33:35
% DurationCPUTime: 0.17s
% Computational Cost: add. (202->55), mult. (354->92), div. (0->0), fcn. (184->8), ass. (0->30)
t97 = m(3) / 0.2e1;
t82 = qJD(1) + qJD(3);
t96 = pkin(7) * t82;
t86 = cos(pkin(8));
t77 = (pkin(1) * t86 + pkin(2)) * qJD(1);
t88 = sin(qJ(3));
t90 = cos(qJ(3));
t84 = sin(pkin(8));
t95 = pkin(1) * qJD(1) * t84;
t73 = t88 * t77 + t90 * t95;
t71 = qJ(4) * t82 + t73;
t83 = sin(pkin(9));
t85 = cos(pkin(9));
t67 = t83 * qJD(2) + t85 * t71;
t72 = t77 * t90 - t88 * t95;
t94 = qJD(4) - t72;
t91 = qJD(2) ^ 2;
t89 = cos(qJ(5));
t87 = sin(qJ(5));
t81 = t85 * qJD(2);
t75 = (t83 * t89 + t85 * t87) * t82;
t74 = (-t83 * t87 + t85 * t89) * t82;
t70 = -pkin(3) * t82 + t94;
t68 = (-pkin(4) * t85 - pkin(3)) * t82 + t94;
t66 = -t71 * t83 + t81;
t65 = t85 * t96 + t67;
t64 = t81 + (-t71 - t96) * t83;
t63 = t64 * t87 + t65 * t89;
t62 = t64 * t89 - t65 * t87;
t1 = m(5) * (t66 ^ 2 + t67 ^ 2 + t70 ^ 2) / 0.2e1 + m(6) * (t62 ^ 2 + t63 ^ 2 + t68 ^ 2) / 0.2e1 + m(4) * (t72 ^ 2 + t73 ^ 2 + t91) / 0.2e1 + t91 * t97 + (t68 * mrSges(6,2) - t62 * mrSges(6,3) + Ifges(6,1) * t75 / 0.2e1) * t75 + (-t68 * mrSges(6,1) + t63 * mrSges(6,3) + Ifges(6,4) * t75 + Ifges(6,2) * t74 / 0.2e1) * t74 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t86 * mrSges(3,1) - t84 * mrSges(3,2) + (t84 ^ 2 + t86 ^ 2) * t97 * pkin(1)) * pkin(1)) * qJD(1) ^ 2 + (t62 * mrSges(6,1) - t63 * mrSges(6,2) + Ifges(6,5) * t75 + Ifges(6,6) * t74 + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (-t73 * mrSges(4,2) + t72 * mrSges(4,1) + t70 * (-mrSges(5,1) * t85 + mrSges(5,2) * t83) + (Ifges(5,2) * t85 ^ 2 / 0.2e1 + Ifges(4,3) / 0.2e1 + (Ifges(5,4) * t85 + Ifges(5,1) * t83 / 0.2e1) * t83) * t82 + (-t66 * t83 + t67 * t85) * mrSges(5,3)) * t82;
T = t1;
