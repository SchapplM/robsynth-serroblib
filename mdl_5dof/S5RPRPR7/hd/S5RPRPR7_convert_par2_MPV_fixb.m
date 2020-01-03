% Return the minimum parameter vector for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
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
% MPV [20x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RPRPR7_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR7_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR7_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t104 = 2 * pkin(7) * mrSges(6,3) + Ifges(6,2);
t97 = pkin(7) ^ 2;
t83 = m(6) * t97 + Ifges(5,1) + t104;
t98 = (pkin(4) ^ 2);
t87 = t98 * m(6) + Ifges(5,2);
t106 = t83 - t87;
t93 = sin(pkin(9));
t95 = cos(pkin(9));
t105 = t93 * t95;
t90 = t93 ^ 2;
t91 = t95 ^ 2;
t103 = t91 - t90;
t101 = pkin(7) * m(6) + mrSges(6,3);
t84 = t101 * pkin(4) + Ifges(5,4);
t102 = t84 * t105;
t85 = mrSges(5,2) - t101;
t88 = m(6) * pkin(4) + mrSges(5,1);
t100 = -t93 * t85 + t95 * t88;
t86 = -pkin(6) * m(4) + mrSges(3,2) - mrSges(4,3);
t89 = m(4) * pkin(2) + mrSges(3,1);
t94 = sin(pkin(8));
t96 = cos(pkin(8));
t99 = -t94 * t86 + t96 * t89;
t1 = [Ifges(2,3) + Ifges(3,3) + Ifges(4,2) + t90 * t83 + 0.2e1 * t102 + t91 * t87 + (2 * pkin(6) * mrSges(4,3)) + ((pkin(2) ^ 2 + pkin(6) ^ 2) * m(4)) + 0.2e1 * t99 * pkin(1); mrSges(2,1) + t99; t96 * t86 + t94 * t89 + mrSges(2,2); m(3) + m(4); t106 * t103 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t102; t103 * t84 + t106 * t105 + Ifges(4,4); t95 * Ifges(5,5) - t93 * Ifges(5,6) + Ifges(4,5); t93 * Ifges(5,5) + t95 * Ifges(5,6) + Ifges(4,6); Ifges(4,3) + Ifges(5,3) + ((t97 + t98) * m(6)) + 0.2e1 * t100 * pkin(3) + t104; mrSges(4,1) + t100; t95 * t85 + t93 * t88 + mrSges(4,2); mrSges(5,3); m(5) + m(6); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
