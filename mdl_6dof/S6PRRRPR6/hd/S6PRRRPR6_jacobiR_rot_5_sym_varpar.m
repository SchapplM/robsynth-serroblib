% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:56
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR6_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:56:40
% EndTime: 2019-02-22 09:56:40
% DurationCPUTime: 0.07s
% Computational Cost: add. (51->24), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->29)
t135 = sin(pkin(6));
t139 = sin(qJ(3));
t151 = t135 * t139;
t142 = cos(qJ(3));
t150 = t135 * t142;
t143 = cos(qJ(2));
t149 = t135 * t143;
t137 = cos(pkin(6));
t140 = sin(qJ(2));
t148 = t137 * t140;
t147 = t137 * t143;
t138 = sin(qJ(4));
t146 = t138 * t142;
t141 = cos(qJ(4));
t145 = t141 * t142;
t144 = t142 * t143;
t136 = cos(pkin(11));
t134 = sin(pkin(11));
t132 = t137 * t139 + t140 * t150;
t131 = t137 * t142 - t140 * t151;
t130 = -t134 * t148 + t136 * t143;
t129 = t134 * t147 + t136 * t140;
t128 = t134 * t143 + t136 * t148;
t127 = t134 * t140 - t136 * t147;
t126 = t130 * t142 + t134 * t151;
t125 = -t130 * t139 + t134 * t150;
t124 = t128 * t142 - t136 * t151;
t123 = -t128 * t139 - t136 * t150;
t1 = [0, -t129 * t145 + t130 * t138, t125 * t141, -t126 * t138 + t129 * t141, 0, 0; 0, -t127 * t145 + t128 * t138, t123 * t141, -t124 * t138 + t127 * t141, 0, 0; 0 (t138 * t140 + t141 * t144) * t135, t131 * t141, -t132 * t138 - t141 * t149, 0, 0; 0, -t129 * t139, t126, 0, 0, 0; 0, -t127 * t139, t124, 0, 0, 0; 0, t139 * t149, t132, 0, 0, 0; 0, -t129 * t146 - t130 * t141, t125 * t138, t126 * t141 + t129 * t138, 0, 0; 0, -t127 * t146 - t128 * t141, t123 * t138, t124 * t141 + t127 * t138, 0, 0; 0 (t138 * t144 - t140 * t141) * t135, t131 * t138, t132 * t141 - t138 * t149, 0, 0;];
JR_rot  = t1;
