% Return the minimum parameter vector for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% MPV [27x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RPRRR13_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_convert_par2_MPV_fixb: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR13_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR13_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t87 = (m(5) + m(6));
t85 = m(4) + t87;
t96 = (pkin(6) * t85);
t88 = (pkin(8) ^ 2);
t90 = (pkin(4) ^ 2);
t95 = (Ifges(5,2) + (t88 + t90) * m(6));
t94 = -pkin(8) * m(6) - mrSges(6,3);
t83 = (mrSges(5,3) - t94);
t93 = pkin(7) * t87 + t83;
t92 = 2 * pkin(8) * mrSges(6,3) + 2 * pkin(7) * t83 + Ifges(6,2) + t95;
t91 = (pkin(3) ^ 2);
t89 = pkin(7) ^ 2;
t1 = [t91 * t87 + Ifges(3,1) + Ifges(4,2) + Ifges(2,3) + (2 * mrSges(4,3) + t96) * pkin(6); mrSges(2,1); mrSges(2,2); mrSges(3,2) - mrSges(4,3) - t96; mrSges(3,3); m(3) + t85; Ifges(4,1) - Ifges(4,2) + (t89 - t91) * t87 + t92; t93 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t89 + t91) * t87 + t92; pkin(3) * t87 + mrSges(4,1); mrSges(4,2) - t93; m(6) * t88 + Ifges(5,1) - t95; Ifges(5,4); t94 * pkin(4) + Ifges(5,5); Ifges(5,6); t90 * m(6) + Ifges(5,3); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
