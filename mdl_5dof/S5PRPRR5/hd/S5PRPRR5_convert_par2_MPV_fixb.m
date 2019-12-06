% Return the minimum parameter vector for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5PRPRR5_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR5_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR5_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR5_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t94 = -pkin(7) * m(6) - mrSges(6,3);
t86 = (m(5) + m(6));
t88 = (pkin(7) ^ 2);
t90 = (pkin(4) ^ 2);
t93 = (Ifges(5,2) + (t88 + t90) * m(6));
t92 = (mrSges(5,3) - t94);
t91 = 2 * pkin(7) * mrSges(6,3) + 2 * pkin(6) * t92 + Ifges(6,2) + t93;
t89 = pkin(6) ^ 2;
t85 = cos(pkin(9));
t84 = sin(pkin(9));
t1 = [m(2) + m(3); Ifges(3,3) + t85 ^ 2 * (Ifges(4,2) + (pkin(3) ^ 2 + t89) * t86 + t91) + (0.2e1 * t85 * Ifges(4,4) + (t89 * t86 + Ifges(4,1) + t91) * t84) * t84; mrSges(3,1); mrSges(3,2); pkin(3) * t86 + mrSges(4,1); mrSges(4,2); pkin(6) * t86 + mrSges(4,3) + t92; m(4) + t86; m(6) * t88 + Ifges(5,1) - t93; Ifges(5,4); t94 * pkin(4) + Ifges(5,5); Ifges(5,6); t90 * m(6) + Ifges(5,3); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
