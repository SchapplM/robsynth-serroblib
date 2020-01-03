% Return the minimum parameter vector for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRPRP5_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_convert_par2_MPV_fixb: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP5_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP5_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP5_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t93 = (pkin(7) ^ 2);
t94 = (pkin(3) ^ 2);
t100 = (-Ifges(5,2) - Ifges(6,3));
t97 = 2 * pkin(7) * mrSges(5,3) - t100;
t84 = Ifges(4,2) + (t93 + t94) * m(5) + t97;
t85 = t93 * m(5) + Ifges(4,1) + t97;
t102 = -t84 + t85;
t91 = sin(pkin(8));
t92 = cos(pkin(8));
t101 = t91 * t92;
t88 = t91 ^ 2;
t89 = t92 ^ 2;
t99 = t89 - t88;
t98 = Ifges(4,4) * t101;
t96 = -pkin(7) * m(5) - mrSges(5,3);
t87 = m(5) * pkin(3) + mrSges(4,1);
t95 = -t91 * mrSges(4,2) + t92 * t87;
t86 = t96 * pkin(3) + Ifges(4,5);
t1 = [Ifges(2,3) + Ifges(3,2) + t88 * t85 + 0.2e1 * t98 + t89 * t84 + (2 * pkin(6) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(6) ^ 2) * m(3)); m(3) * pkin(1) + mrSges(2,1); -pkin(6) * m(3) + mrSges(2,2) - mrSges(3,3); t102 * t99 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t98; t99 * Ifges(4,4) + t102 * t101 + Ifges(3,4); -t91 * Ifges(4,6) + t92 * t86 + Ifges(3,5); t92 * Ifges(4,6) + t91 * t86 + Ifges(3,6); (t94 * m(5)) + 0.2e1 * pkin(2) * t95 + Ifges(3,3) + Ifges(4,3); mrSges(3,1) + t95; t92 * mrSges(4,2) + t91 * t87 + mrSges(3,2); mrSges(4,3) - t96; m(4) + m(5); Ifges(5,1) + Ifges(6,1) + t100; Ifges(5,4) - Ifges(6,5); Ifges(5,5) + Ifges(6,4); Ifges(5,6) - Ifges(6,6); Ifges(5,3) + Ifges(6,2); mrSges(5,1); mrSges(5,2); mrSges(6,1); mrSges(6,2); mrSges(6,3); m(6);];
MPV = t1;
