% Return the minimum parameter vector for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% MPV [23x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRRPP4_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_convert_par2_MPV_fixb: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP4_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP4_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t86 = sin(pkin(8));
t84 = t86 ^ 2;
t87 = cos(pkin(8));
t85 = t87 ^ 2;
t102 = t85 - t84;
t101 = t86 * t87;
t89 = Ifges(5,2) + Ifges(6,3);
t92 = Ifges(5,1) + Ifges(6,1);
t100 = t89 - t92;
t91 = Ifges(5,4) - Ifges(6,5);
t99 = t91 * t101;
t98 = -pkin(7) * m(4) - mrSges(4,3);
t97 = (mrSges(3,3) - t98);
t96 = t87 * mrSges(5,1) - t86 * mrSges(5,2);
t95 = pkin(2) ^ 2;
t94 = pkin(7) ^ 2;
t93 = (m(3) + m(4));
t90 = Ifges(5,5) + Ifges(6,4);
t88 = Ifges(5,6) - Ifges(6,6);
t83 = t94 + t95;
t1 = [Ifges(2,3) + Ifges(3,2) + Ifges(4,2) + t84 * t92 + 0.2e1 * t99 + t85 * t89 + (2 * pkin(7) * mrSges(4,3)) + t83 * m(4) + (2 * pkin(6) * t97) + ((pkin(1) ^ 2 + pkin(6) ^ 2) * t93); pkin(1) * t93 + mrSges(2,1); -pkin(6) * t93 + mrSges(2,2) - t97; Ifges(3,1) - Ifges(3,2) + (-t83 + t94) * m(4); Ifges(3,4); t98 * pkin(2) + Ifges(3,5); Ifges(3,6); t95 * m(4) + Ifges(3,3); m(4) * pkin(2) + mrSges(3,1); mrSges(3,2); -t102 * t100 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t99; -t100 * t101 + t102 * t91 + Ifges(4,4); -t86 * t88 + t87 * t90 + Ifges(4,5); t86 * t90 + t87 * t88 + Ifges(4,6); 0.2e1 * pkin(3) * t96 + Ifges(6,2) + Ifges(4,3) + Ifges(5,3); mrSges(4,1) + t96; t86 * mrSges(5,1) + t87 * mrSges(5,2) + mrSges(4,2); mrSges(5,3); m(5); mrSges(6,1); mrSges(6,2); mrSges(6,3); m(6);];
MPV = t1;
