% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:59
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRP6_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:59:24
% EndTime: 2019-02-22 11:59:25
% DurationCPUTime: 0.15s
% Computational Cost: add. (125->31), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->37)
t143 = cos(pkin(6));
t145 = sin(qJ(2));
t149 = cos(qJ(1));
t151 = t149 * t145;
t146 = sin(qJ(1));
t148 = cos(qJ(2));
t153 = t146 * t148;
t134 = t143 * t151 + t153;
t141 = qJ(3) + pkin(11);
t139 = sin(t141);
t140 = cos(t141);
t142 = sin(pkin(6));
t156 = t142 * t149;
t128 = -t134 * t140 + t139 * t156;
t150 = t149 * t148;
t154 = t146 * t145;
t133 = -t143 * t150 + t154;
t144 = sin(qJ(5));
t147 = cos(qJ(5));
t164 = t128 * t144 + t133 * t147;
t163 = t128 * t147 - t133 * t144;
t160 = t140 * t144;
t159 = t140 * t147;
t158 = t142 * t145;
t157 = t142 * t146;
t155 = t144 * t148;
t152 = t147 * t148;
t126 = -t134 * t139 - t140 * t156;
t136 = -t143 * t154 + t150;
t135 = t143 * t153 + t151;
t132 = t143 * t139 + t140 * t158;
t131 = -t139 * t158 + t143 * t140;
t130 = t136 * t140 + t139 * t157;
t129 = t136 * t139 - t140 * t157;
t125 = t130 * t147 + t135 * t144;
t124 = -t130 * t144 + t135 * t147;
t1 = [t163, -t135 * t159 + t136 * t144, -t129 * t147, 0, t124, 0; t125, -t133 * t159 + t134 * t144, t126 * t147, 0, t164, 0; 0 (t140 * t152 + t144 * t145) * t142, t131 * t147, 0, -t132 * t144 - t142 * t152, 0; -t164, t135 * t160 + t136 * t147, t129 * t144, 0, -t125, 0; t124, t133 * t160 + t134 * t147, -t126 * t144, 0, t163, 0; 0 (-t140 * t155 + t145 * t147) * t142, -t131 * t144, 0, -t132 * t147 + t142 * t155, 0; t126, -t135 * t139, t130, 0, 0, 0; t129, -t133 * t139, -t128, 0, 0, 0; 0, t142 * t148 * t139, t132, 0, 0, 0;];
JR_rot  = t1;
