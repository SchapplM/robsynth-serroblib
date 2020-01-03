% Return the minimum parameter vector for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RPRPR14_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_convert_par2_MPV_fixb: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR14_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR14_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR14_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t88 = (pkin(7) ^ 2);
t94 = 2 * pkin(7) * mrSges(6,3) + Ifges(6,2);
t78 = m(6) * t88 + Ifges(5,1) + t94;
t89 = (pkin(4) ^ 2);
t81 = t89 * m(6) + Ifges(5,2);
t97 = t78 - t81;
t96 = (pkin(6) * m(4));
t86 = sin(pkin(8));
t87 = cos(pkin(8));
t95 = t86 * t87;
t83 = t86 ^ 2;
t84 = t87 ^ 2;
t93 = t84 - t83;
t91 = pkin(7) * m(6) + mrSges(6,3);
t79 = t91 * pkin(4) + Ifges(5,4);
t92 = t79 * t95;
t80 = mrSges(5,2) - t91;
t82 = m(6) * pkin(4) + mrSges(5,1);
t90 = -t86 * t80 + t87 * t82;
t1 = [0.2e1 * t92 + t83 * t78 + t84 * t81 + Ifges(3,1) + Ifges(4,2) + Ifges(2,3) + ((2 * mrSges(4,3) + t96) * pkin(6)); mrSges(2,1); mrSges(2,2); mrSges(3,2) - mrSges(4,3) - t96; mrSges(3,3); m(3) + m(4); t97 * t93 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t92; t93 * t79 + t97 * t95 + Ifges(4,4); t87 * Ifges(5,5) - t86 * Ifges(5,6) + Ifges(4,5); t86 * Ifges(5,5) + t87 * Ifges(5,6) + Ifges(4,6); Ifges(4,3) + Ifges(5,3) + ((t88 + t89) * m(6)) + 0.2e1 * t90 * pkin(3) + t94; mrSges(4,1) + t90; t87 * t80 + t86 * t82 + mrSges(4,2); mrSges(5,3); m(5) + m(6); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
