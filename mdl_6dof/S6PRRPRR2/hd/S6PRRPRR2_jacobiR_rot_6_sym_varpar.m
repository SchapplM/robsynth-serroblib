% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:48
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRR2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:48:37
% EndTime: 2019-02-22 09:48:37
% DurationCPUTime: 0.08s
% Computational Cost: add. (167->29), mult. (219->64), div. (0->0), fcn. (320->10), ass. (0->38)
t145 = qJ(3) + pkin(12);
t142 = cos(t145);
t146 = qJ(5) + qJ(6);
t143 = sin(t146);
t161 = t142 * t143;
t144 = cos(t146);
t160 = t142 * t144;
t152 = cos(qJ(2));
t159 = t142 * t152;
t147 = sin(pkin(11));
t148 = sin(pkin(6));
t158 = t147 * t148;
t149 = cos(pkin(11));
t157 = t148 * t149;
t151 = sin(qJ(2));
t156 = t148 * t151;
t155 = t148 * t152;
t150 = cos(pkin(6));
t154 = t150 * t151;
t153 = t150 * t152;
t141 = sin(t145);
t139 = -t147 * t154 + t149 * t152;
t138 = t147 * t153 + t149 * t151;
t137 = t147 * t152 + t149 * t154;
t136 = t147 * t151 - t149 * t153;
t135 = t150 * t141 + t142 * t156;
t134 = -t141 * t156 + t150 * t142;
t133 = t139 * t142 + t141 * t158;
t132 = -t139 * t141 + t142 * t158;
t131 = t137 * t142 - t141 * t157;
t130 = -t137 * t141 - t142 * t157;
t129 = -t135 * t144 + t143 * t155;
t128 = -t135 * t143 - t144 * t155;
t127 = -t133 * t144 - t138 * t143;
t126 = -t133 * t143 + t138 * t144;
t125 = -t131 * t144 - t136 * t143;
t124 = -t131 * t143 + t136 * t144;
t1 = [0, -t138 * t160 + t139 * t143, t132 * t144, 0, t126, t126; 0, -t136 * t160 + t137 * t143, t130 * t144, 0, t124, t124; 0 (t143 * t151 + t144 * t159) * t148, t134 * t144, 0, t128, t128; 0, t138 * t161 + t139 * t144, -t132 * t143, 0, t127, t127; 0, t136 * t161 + t137 * t144, -t130 * t143, 0, t125, t125; 0 (-t143 * t159 + t144 * t151) * t148, -t134 * t143, 0, t129, t129; 0, -t138 * t141, t133, 0, 0, 0; 0, -t136 * t141, t131, 0, 0, 0; 0, t141 * t155, t135, 0, 0, 0;];
JR_rot  = t1;
