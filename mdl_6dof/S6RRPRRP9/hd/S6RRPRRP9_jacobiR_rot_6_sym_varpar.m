% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:36
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRP9_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:36:03
% EndTime: 2019-02-22 11:36:03
% DurationCPUTime: 0.15s
% Computational Cost: add. (125->31), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->37)
t144 = cos(pkin(6));
t146 = sin(qJ(2));
t150 = cos(qJ(1));
t152 = t150 * t146;
t147 = sin(qJ(1));
t149 = cos(qJ(2));
t154 = t147 * t149;
t135 = t144 * t152 + t154;
t142 = pkin(11) + qJ(4);
t140 = sin(t142);
t141 = cos(t142);
t143 = sin(pkin(6));
t157 = t143 * t150;
t129 = -t135 * t141 + t140 * t157;
t151 = t150 * t149;
t155 = t147 * t146;
t134 = -t144 * t151 + t155;
t145 = sin(qJ(5));
t148 = cos(qJ(5));
t165 = t129 * t145 + t134 * t148;
t164 = t129 * t148 - t134 * t145;
t161 = t141 * t145;
t160 = t141 * t148;
t159 = t143 * t146;
t158 = t143 * t147;
t156 = t145 * t149;
t153 = t148 * t149;
t127 = -t135 * t140 - t141 * t157;
t137 = -t144 * t155 + t151;
t136 = t144 * t154 + t152;
t133 = t144 * t140 + t141 * t159;
t132 = -t140 * t159 + t144 * t141;
t131 = t137 * t141 + t140 * t158;
t130 = t137 * t140 - t141 * t158;
t126 = t131 * t148 + t136 * t145;
t125 = -t131 * t145 + t136 * t148;
t1 = [t164, -t136 * t160 + t137 * t145, 0, -t130 * t148, t125, 0; t126, -t134 * t160 + t135 * t145, 0, t127 * t148, t165, 0; 0 (t141 * t153 + t145 * t146) * t143, 0, t132 * t148, -t133 * t145 - t143 * t153, 0; -t165, t136 * t161 + t137 * t148, 0, t130 * t145, -t126, 0; t125, t134 * t161 + t135 * t148, 0, -t127 * t145, t164, 0; 0 (-t141 * t156 + t146 * t148) * t143, 0, -t132 * t145, -t133 * t148 + t143 * t156, 0; t127, -t136 * t140, 0, t131, 0, 0; t130, -t134 * t140, 0, -t129, 0, 0; 0, t143 * t149 * t140, 0, t133, 0, 0;];
JR_rot  = t1;
