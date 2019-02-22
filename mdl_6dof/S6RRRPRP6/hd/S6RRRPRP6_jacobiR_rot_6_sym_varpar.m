% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function JR_rot = S6RRRPRP6_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:59:15
% EndTime: 2019-02-22 11:59:15
% DurationCPUTime: 0.16s
% Computational Cost: add. (125->31), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->37)
t145 = cos(pkin(6));
t147 = sin(qJ(2));
t151 = cos(qJ(1));
t153 = t151 * t147;
t148 = sin(qJ(1));
t150 = cos(qJ(2));
t155 = t148 * t150;
t136 = t145 * t153 + t155;
t143 = qJ(3) + pkin(11);
t141 = sin(t143);
t142 = cos(t143);
t144 = sin(pkin(6));
t158 = t144 * t151;
t130 = -t136 * t142 + t141 * t158;
t152 = t151 * t150;
t156 = t148 * t147;
t135 = -t145 * t152 + t156;
t146 = sin(qJ(5));
t149 = cos(qJ(5));
t166 = t130 * t146 + t135 * t149;
t165 = t130 * t149 - t135 * t146;
t162 = t142 * t146;
t161 = t142 * t149;
t160 = t144 * t147;
t159 = t144 * t148;
t157 = t146 * t150;
t154 = t149 * t150;
t128 = -t136 * t141 - t142 * t158;
t138 = -t145 * t156 + t152;
t137 = t145 * t155 + t153;
t134 = t145 * t141 + t142 * t160;
t133 = -t141 * t160 + t145 * t142;
t132 = t138 * t142 + t141 * t159;
t131 = t138 * t141 - t142 * t159;
t127 = t132 * t149 + t137 * t146;
t126 = -t132 * t146 + t137 * t149;
t1 = [t165, -t137 * t161 + t138 * t146, -t131 * t149, 0, t126, 0; t127, -t135 * t161 + t136 * t146, t128 * t149, 0, t166, 0; 0 (t142 * t154 + t146 * t147) * t144, t133 * t149, 0, -t134 * t146 - t144 * t154, 0; -t166, t137 * t162 + t138 * t149, t131 * t146, 0, -t127, 0; t126, t135 * t162 + t136 * t149, -t128 * t146, 0, t165, 0; 0 (-t142 * t157 + t147 * t149) * t144, -t133 * t146, 0, -t134 * t149 + t144 * t157, 0; t128, -t137 * t141, t132, 0, 0, 0; t131, -t135 * t141, -t130, 0, 0, 0; 0, t144 * t150 * t141, t134, 0, 0, 0;];
JR_rot  = t1;
