% Return the minimum parameter vector for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% MPV [22x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRRPR3_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR3_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t88 = (pkin(8) ^ 2);
t89 = (pkin(4) ^ 2);
t94 = 2 * pkin(8) * mrSges(6,3) + Ifges(6,2);
t78 = Ifges(5,2) + (t88 + t89) * m(6) + t94;
t79 = m(6) * t88 + Ifges(5,1) + t94;
t96 = -t78 + t79;
t85 = sin(pkin(9));
t86 = cos(pkin(9));
t95 = t85 * t86;
t82 = t85 ^ 2;
t83 = t86 ^ 2;
t93 = t83 - t82;
t92 = Ifges(5,4) * t95;
t91 = -pkin(8) * m(6) - mrSges(6,3);
t81 = m(6) * pkin(4) + mrSges(5,1);
t90 = -t85 * mrSges(5,2) + t86 * t81;
t87 = m(3) + m(4);
t80 = t91 * pkin(4) + Ifges(5,5);
t1 = [pkin(1) ^ 2 * t87 + Ifges(2,3); pkin(1) * t87 + mrSges(2,1); mrSges(2,2); Ifges(3,3) + Ifges(4,2) + t82 * t79 + 0.2e1 * t92 + t83 * t78 + (2 * pkin(7) * mrSges(4,3)) + ((pkin(2) ^ 2 + pkin(7) ^ 2) * m(4)); m(4) * pkin(2) + mrSges(3,1); -pkin(7) * m(4) + mrSges(3,2) - mrSges(4,3); t96 * t93 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t92; t93 * Ifges(5,4) + t96 * t95 + Ifges(4,4); -t85 * Ifges(5,6) + t86 * t80 + Ifges(4,5); t86 * Ifges(5,6) + t85 * t80 + Ifges(4,6); (t89 * m(6)) + 0.2e1 * pkin(3) * t90 + Ifges(4,3) + Ifges(5,3); mrSges(4,1) + t90; t86 * mrSges(5,2) + t85 * t81 + mrSges(4,2); mrSges(5,3) - t91; m(5) + m(6); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
