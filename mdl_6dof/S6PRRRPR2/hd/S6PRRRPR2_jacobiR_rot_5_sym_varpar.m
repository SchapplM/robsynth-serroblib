% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:54
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR2_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:54:38
% EndTime: 2019-02-22 09:54:38
% DurationCPUTime: 0.07s
% Computational Cost: add. (97->23), mult. (158->52), div. (0->0), fcn. (231->10), ass. (0->36)
t134 = sin(pkin(11));
t137 = cos(pkin(11));
t140 = cos(qJ(2));
t138 = cos(pkin(6));
t139 = sin(qJ(2));
t142 = t138 * t139;
t127 = t134 * t140 + t137 * t142;
t132 = qJ(3) + qJ(4);
t130 = sin(t132);
t131 = cos(t132);
t135 = sin(pkin(6));
t144 = t135 * t137;
t119 = -t127 * t130 - t131 * t144;
t133 = sin(pkin(12));
t151 = t119 * t133;
t129 = -t134 * t142 + t137 * t140;
t145 = t134 * t135;
t121 = -t129 * t130 + t131 * t145;
t150 = t121 * t133;
t143 = t135 * t139;
t124 = -t130 * t143 + t138 * t131;
t149 = t124 * t133;
t148 = t131 * t133;
t136 = cos(pkin(12));
t147 = t131 * t136;
t146 = t131 * t140;
t141 = t138 * t140;
t128 = -t134 * t141 - t137 * t139;
t126 = -t134 * t139 + t137 * t141;
t125 = t138 * t130 + t131 * t143;
t123 = t124 * t136;
t122 = t129 * t131 + t130 * t145;
t120 = t127 * t131 - t130 * t144;
t118 = t121 * t136;
t117 = t119 * t136;
t1 = [0, t128 * t147 + t129 * t133, t118, t118, 0, 0; 0, t126 * t147 + t127 * t133, t117, t117, 0, 0; 0 (t133 * t139 + t136 * t146) * t135, t123, t123, 0, 0; 0, -t128 * t148 + t129 * t136, -t150, -t150, 0, 0; 0, -t126 * t148 + t127 * t136, -t151, -t151, 0, 0; 0 (-t133 * t146 + t136 * t139) * t135, -t149, -t149, 0, 0; 0, t128 * t130, t122, t122, 0, 0; 0, t126 * t130, t120, t120, 0, 0; 0, t135 * t140 * t130, t125, t125, 0, 0;];
JR_rot  = t1;
